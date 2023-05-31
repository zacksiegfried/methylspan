import csv
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

### Could add platform specific mapping either by buidling other functions and subsetting sample data set
# or by tagging sample data set with platform type (not sure if this makes programatic sense)
def methylMap450(methylation_beta_array, drop_missing = 'Y'):
    """Transforms column names from IlmnID to UCSC_RefGene_Name using 450 manifest,
        drops CpG sites with no mapping information by default"""

    # sets files directory as working directory
    os.chdir(script_dir)
    df_450 = csv.DictReader(open('mapping_files/humanmethylation450.csv'))

    mapping_dict = {}
    for row in df_450:
        mapping_dict[row['IlmnID']] = row['UCSC_RefGene_Name']
    
    methylation_beta_array.rename(columns=mapping_dict, inplace=True)

    if drop_missing == 'Y':
        methylation_beta_array = methylation_beta_array.drop([""], axis=1)

    return(methylation_beta_array)


def averageGene(mapped_sample_data):
    """Averages mapped methylation beta array beta values for all similar genes (not scientific)"""

    df = mapped_sample_data.groupby(by=mapped_sample_data.columns, axis=1).mean()

    return(df)


def fullGene(mapped_sample_data, gene):
    """Filters mapped methylation beta array by a specificed gene and identifys indivdual CpG sites"""

    df = mapped_sample_data[gene]
    cols = []
    count = 1
    for column in df.columns:
        if column == gene:
            cols.append(f'{gene}_{count}')
            count+=1
            continue
        cols.append(column)
    df.columns = cols

    return(df)