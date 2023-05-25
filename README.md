# methylspan

## NCI GDC Data Portal Tools

### Accessing methylation data

A complete methylation beta array containing CpG information for all patients of a specified cancer primary site can be created using the MethylDataFetch.getMethylBetaArrays function

The following will return a pandas dataframe for all patients with Thymus as primary cancer site. The first column 'cpg' contains Ilumina CpG index followed by columns named after patient 'file_uuid' containing methylation beta values 

```
import sys
sys.path.insert(0, '/Users/user_name/Documents/methylspan')

from MethylDataFetch import getMethylBetaArrays

frame = getMethylBetaArrays('thymus')
```

## Survival analysis on cancer patients using methylation data

### Methylation data management

MethylDataFetch.py and MSM.py must be placed in the same directory. Upon downloading methylation data the newly made beta_arrays directory must be kept in the same directory as these two scripts.

MethylDataFetch.py contains functions related to downloading methylation data from NCI GDC data portal. 
Run the script and input primary site name to download all methylation beta array files from cases with corresponding primary site.

**_Script will display download size and ask if you wish to proceed before downloading_**

MSM.py contains functions to import and format methylation beta arrays and meta data into pandas data frames. These are mainly used in the modeling notebooks to import data.


### Mapping Ilumina CpG IDs to genomic positions

MethylMapping.py contains functions for mapping CpG ids to genes.


### modeling 

notebooks is a directory of unfinished jupyter notebooks relating to modeling the methylation and meta data
