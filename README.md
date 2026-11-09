# methylspan

Download DNA methylation data from the [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/) and test whether it predicts cancer survival beyond clinical data.

The project began as a 2023 MPH practicum on pancreatic cancer survival. That paper and slides are in [docs/2023-practicum/](docs/2023-practicum/), and the code used for them is tagged `v2023-practicum`. The current code fixes problems found in that analysis; see [What changed since 2023](#what-changed-since-2023).

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e '.[dev,notebooks]'
nbstripout --install --attributes .gitattributes   # strip notebook outputs on commit
```

Requires Python 3.10+.

## Workflow

Each step reads the previous step's output, so you only download once.

```bash
# 1. See what GDC has for a site, then select a cohort.
#    Defaults: primary tumors only, one file per patient. Prints counts by project/platform/sample type.
methylspan files --site pancreas --platform 450K --out data/paad_450k_files.csv

# 2. Download (asks before starting; resumes if interrupted)
methylspan download data/paad_450k_files.csv --raw-dir data/raw

# 3. Build a probes x patients float32 matrix, cached as parquet
methylspan matrix data/paad_450k_files.csv --raw-dir data/raw --out data/paad_450k_betas.parquet

# 4. Survival table (time, event, age, stage, ...) for the same patients
methylspan clinical data/paad_450k_files.csv --out data/paad_clinical.csv
```

GDC periodically reprocesses files and gives them new IDs. To reuse an old download, pass `--resolve-superseded` to `methylspan matrix`. It looks up each old folder's current ID through the GDC `/history` endpoint:

```bash
methylspan matrix data/paad_450k_files.csv \
  --raw-dir path/to/gdc_download_YYYYMMDD --resolve-superseded \
  --out data/paad_450k_betas.parquet
```

Then model in Python. [notebooks/pdac_survival.ipynb](notebooks/pdac_survival.ipynb) has the full example:

```python
from methylspan import matrix, survival
import pandas as pd

clinical = pd.read_csv("data/paad_clinical.csv", index_col="case_id")
betas = matrix.filter_probes(matrix.load_matrix("data/paad_450k_betas.parquet"))
X, y = survival.assemble(clinical, betas.T)

models = {
    "clinical": survival.make_model(use_methylation=False),
    "clinical + coxnet": survival.tune(survival.make_model(n_top=5000, estimator=survival.coxnet())),
}
survival.summarize(survival.compare_models(models, X, y, n_splits=5, n_repeats=3))
```

## Results: TCGA-PAAD

The cohort is 183 primary tumors on the 450K array, one per patient, with 99 deaths. After dropping probes missing in more than 10% of samples, 393,592 probes remain. Models were compared with 5-fold cross-validation repeated 3 times (15 folds), using the same folds for every model.

| Model | Harrell's C | Uno's C | Mean time-dep. AUC | Integrated Brier score |
|---|---|---|---|---|
| Clinical only (age + stage) | 0.547 ± 0.074 | 0.543 | 0.548 | 0.190 |
| Methylation, 15 PCs (2023 approach) | **0.610** ± 0.064 | 0.605 | 0.621 | 0.191 |
| Clinical + 15 PCs | 0.593 ± 0.075 | 0.589 | 0.604 | 0.198 |
| Clinical + elastic-net Cox (top 5k probes) | 0.585 ± 0.063 | 0.583 | 0.611 | 0.190 |

What this shows:

- **The 2023 C-index of about 0.66 falls to 0.61** once the cohort is cleaned and preprocessing is fit inside each fold.
- **Methylation does add some ranking ability.** The 15-PC model beats the clinical model on 13 of 15 matched folds, by +0.06 C on average. That is a small effect, and fold-to-fold variation (±0.07) is about as large.
- **Clinical data is a weak baseline here.** 151 of the 183 patients are stage II, so stage barely separates them.
- **None of the models predicts survival time better.** The integrated Brier scores are essentially identical (about 0.19), so methylation helps rank patients by risk but doesn't improve predicted survival curves.
- **The cohort includes 8 neuroendocrine tumors.** Restricting to ductal adenocarcinoma is a sensible next run.

## Modules

| Module | What it does |
|---|---|
| `gdc` | Paginated `/files` and `/cases` queries, streamed batch downloads with safe tar extraction, and superseded-ID lookup |
| `cohort` | `describe()` composition table; `select_cohort()` filters by sample type, platform and project, keeping one file per patient |
| `matrix` | Reads SeSAMe beta files into a probes x patients float32 matrix; parquet save/load; missingness filter |
| `clinical` | Survival table using the primary diagnosis and the latest follow-up across all records |
| `mapping` | Illumina 450K/EPIC manifest loader, probe-to-gene mapping by region (promoter by default), gene-level summaries |
| `survival` | Cox pipelines with clinical and/or methylation features, in-fold preprocessing, nested alpha tuning, and Harrell/Uno C, time-dependent AUC and integrated Brier score |

Probe-to-gene mapping needs an Illumina manifest CSV (450K or EPIC) from Illumina's support site. Data files are not tracked by git.

## What changed since 2023

The 2023 analysis queried every pancreas methylation file and pooled them. GDC now returns 565 such files from 445 patients. They include adjacent normal tissue, lab-grown organoid models (HCMI), 450K and EPIC arrays, and 108 patients with more than one file. So the models saw duplicated patients and large tissue/platform differences that have nothing to do with prognosis. Other fixes:

- Survival time comes from the diagnosis flagged as primary; the old code used the first one listed. Alive patients use their latest follow-up from any record.
- Probe selection, scaling and PCA are fit on training folds only.
- The old "median days error" metric included censored patients and scored failed (infinite) predictions as perfect. It is replaced by standard censoring-aware metrics.
- Gene mapping splits multi-gene annotations, pairs each gene with its region, and defaults to promoter probes. It no longer averages promoter and gene-body probes together.
- GDC queries paginate (the old code stopped at 1,000 results). A leftover `[:10]` test slice limited every "full" fetch to 10 files. Downloads stream to disk in resumable batches.
- The package installs with `pip`, has tests (`pytest`), and no longer relies on `sys.path` hacks, `os.chdir` or hard-coded paths.

## Tests

```bash
pytest
```

Tests use synthetic data and a fake HTTP session, so they run offline.
