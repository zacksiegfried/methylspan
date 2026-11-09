"""Command line entry point: `methylspan files|download|matrix|clinical`."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from . import clinical, cohort, gdc, matrix


def _add_cohort_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--sample-type", action="append", dest="sample_types",
                   help="Keep this sample type (repeatable). Default: 'Primary Tumor'. Use 'all' to keep everything.")
    p.add_argument("--platform", action="append", dest="platforms", choices=["450K", "EPIC", "EPICv2", "27K"],
                   help="Keep this platform (repeatable). Default: all, with a warning if mixed.")
    p.add_argument("--keep-duplicates", action="store_true",
                   help="Keep several files per patient (not suitable for survival models).")


def cmd_files(args: argparse.Namespace) -> None:
    files = gdc.query_methylation_files(primary_site=args.site, project_ids=args.project)
    print("All methylation beta files:\n")
    print(cohort.describe(files).to_string(), "\n")
    sample_types = None if args.sample_types == ["all"] else (args.sample_types or ["Primary Tumor"])
    selected = cohort.select_cohort(files, sample_types=sample_types, platforms=args.platforms,
                                    one_per_case=not args.keep_duplicates)
    print("Selected cohort:\n")
    print(cohort.describe(selected).to_string())
    print(f"\n{len(selected)} files, {selected['file_size'].sum() / 1e9:.2f} GB")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")


def cmd_download(args: argparse.Namespace) -> None:
    files = pd.read_csv(args.files)
    print(f"{len(files)} files, {files['file_size'].sum() / 1e9:.2f} GB -> {args.raw_dir}")
    if not args.yes and input("Download? [y/N] ").strip().lower() != "y":
        return
    gdc.download_files(files, args.raw_dir)


def cmd_matrix(args: argparse.Namespace) -> None:
    files = pd.read_csv(args.files)
    aliases = None
    if args.resolve_superseded:
        old_ids = matrix.unmatched_dirs(files, args.raw_dir)
        print(f"Looking up current versions of {len(old_ids)} unmatched folders")
        current = gdc.resolve_superseded(old_ids)
        aliases = {new: old for old, new in current.items()}
    m = matrix.build_beta_matrix(files, args.raw_dir, join=args.join, aliases=aliases)
    matrix.save_matrix(m, args.out)
    print(f"Wrote {m.shape[0]} probes x {m.shape[1]} samples to {args.out}")


def cmd_clinical(args: argparse.Namespace) -> None:
    files = pd.read_csv(args.files)
    table = clinical.survival_table(gdc.query_cases(files["case_id"].dropna().unique()))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out)
    usable = clinical.usable_for_survival(table)
    print(f"Wrote {len(table)} patients to {args.out} ({len(usable)} usable, {int(usable['event'].sum())} deaths)")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="methylspan", description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("files", help="List GDC methylation files and select a cohort")
    p.add_argument("--site", help="GDC primary site, e.g. pancreas")
    p.add_argument("--project", action="append", help="GDC project id, e.g. TCGA-PAAD (repeatable)")
    _add_cohort_args(p)
    p.add_argument("--out", default="data/files.csv")
    p.set_defaults(func=cmd_files)

    p = sub.add_parser("download", help="Download the files listed in a cohort CSV")
    p.add_argument("files")
    p.add_argument("--raw-dir", default="data/raw")
    p.add_argument("-y", "--yes", action="store_true", help="Skip the size confirmation")
    p.set_defaults(func=cmd_download)

    p = sub.add_parser("matrix", help="Build a probes x patients parquet matrix from downloaded files")
    p.add_argument("files")
    p.add_argument("--raw-dir", default="data/raw")
    p.add_argument("--join", choices=["inner", "outer"], default="inner")
    p.add_argument("--resolve-superseded", action="store_true",
                   help="Match folders named by old GDC file ids (e.g. a 2023 download) to current ids")
    p.add_argument("--out", default="data/betas.parquet")
    p.set_defaults(func=cmd_matrix)

    p = sub.add_parser("clinical", help="Fetch a survival table for the patients in a cohort CSV")
    p.add_argument("files")
    p.add_argument("--out", default="data/clinical.csv")
    p.set_defaults(func=cmd_clinical)

    args = parser.parse_args(argv)
    if args.command == "files" and not (args.site or args.project):
        parser.error("files needs --site and/or --project")
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s %(message)s")
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
