# methylspan

## NCI GDC Data Portal Tools

### Accessing methylation data (No permanent download method)

A complete methylation beta array containing CpG information for all patients of a specified cancer primary site can be created using the MethylDataFetch.getMethylBetaArrays function

A new directory 'data' is created where NCI GDC data files are downloaded, read into python, and then deleted

```
import sys
sys.path.insert(0, '/Users/user_name/Documents/methylspan')

from MethylDataFetch import getMethylBetaArrays

data = getMethylBetaArrays('thymus')
```
Returns a dictionary in the form ```{'data': df, 'primary_site': 'thymus'}```

The dataframe contains methylation information for all patients with Thymus as primary cancer site. Rows contains Ilumina CpG index. Columns are named the unique patient 'file_uuid' and contain methylation beta values. Saves a csv file under methylspan/analysis/comp.csv with the same information. Conditional formatting can be used to quickly highlight any genes with significant differences in methylation content.

### Accessing methylation data (Permanent download method)

MethylDataFetch.py contains functions related to downloading methylation data from NCI GDC data portal. 
Run the script and input primary site name to download all methylation beta array files from cases with corresponding primary site into the methylspan/beta_arrays directory.

**_Script will display download size and ask if you wish to proceed before downloading_**

## Methylation Analysis

AnalysisFunctions.py contains complete analysis functions that are designed to go from input primary site(s) to interpretable output data.

AnalysisFunctions.visualGeneMethylContent(primary_site_list)
Takes a list object of primary site names and returns a dataframe with indivdual gene IDs as rows. Columns contain primary site name with mean methylation content data.

## Survival analysis on cancer patients using methylation data

### Mapping Ilumina CpG IDs to genomic positions

MethylMapping.py contains functions for mapping CpG ids to genes.


### modeling 

notebooks is a directory of unfinished jupyter notebooks relating to modeling the methylation and meta data
