# methylspan

## Survival analysis on cancer patients using methylation data

### methylation data management

MethylDataFetch.py and MSM.py must be placed in the same directory. Upon downloading methylation data the newly made beta_arrays directory must be kept in the same directory as these two scripts.

MethylDataFetch.py contains functions related to downloading methylation data from NCI GDC data portal. 
Run the script and input primary site name to download all methylation beta array files from cases with corresponding primary site.

**_Script will display download size and ask if you wish to proceed before downloading_**

MSM.py contains functions to import and format methylation beta arrays and meta data into pandas data frames. These are mainly used in the modeling notebooks to import data.


### modeling 

notebooks is a directory of unfinished jupyter notebooks relating to modeling the methylation and meta data
