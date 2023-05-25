import MethylDataFetch
import os
import random
import pandas as pd
from functools import reduce

"""Contains functions related to formatting metadata and methylation beta arrays"""


def pullMethylMetaData(primary_site):
    """Returns list of dictionaries with metadata matching primary site cases that are dead with days to death data"""

    full_file_info = MethylDataFetch.getMethylMetaData(str(primary_site))

    ### could build into pull "Dead", pull "Alive" functions if useful
    #if file["cases"][0]["demographic"]["vital_status"] == str(vital_status):
        #if file["cases"][0]["demographic"]["days_to_death"] != None:
    
    print("Number of methylation files for all cases for primary site " + str(primary_site) + ": " + str(len(full_file_info)))
    return(full_file_info)


def methylDataFormat(all_files, sample_size = None):
    """Returns methylation array as a df. Index column is file_is, column headers are CpG site #"""

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Uses all the data by default, allows reduction of cases to speed up output
    if sample_size == None:
        sample_set = all_files
    else:
        sample_set = random.sample(all_files, sample_size)

    # find sub directory name
    for file in os.listdir(str(script_dir + "/beta_arrays")):
        d = os.path.join(str(script_dir + "/beta_arrays"), file)
        if os.path.isdir(d):
            sub_folder = d


    # makes list of dataframes that match the file_names in sample_set (SLOW)
    sample_df_list = []
    # iterates over list of metadata dictonaires that have been randomly selected (=sample_no = 100)
    for entry in sample_set:
        # iterates over the files inside the sub directory of beta_arrays (should all be directories with generated names)
        for file in os.listdir(str(sub_folder)):
            d = os.path.join(str(sub_folder), file)
            if os.path.isdir(d): # if file is directory
                # for file in sub sub directory, if file name == file name in training set metadata add to df list
                for data in os.listdir(str(d)):
                    if str(data) == entry['file_name']:
                        file_path = str(d + "/" + data)
                        df = pd.read_csv(file_path, sep="\t", header=None)
                        df.columns = ['file_id', entry['file_id']]
                        sample_df_list.append(df)
        


    # merging methylation datasets, dropping cases with lots of missing values and CpG sites with no reads
    sample_df = reduce(lambda x, y: pd.merge(x, y, on = 'file_id'), sample_df_list)
    sample_df = sample_df.transpose()
    sample_df.columns = sample_df.iloc[0]
    sample_df = sample_df.iloc[1:]
    sample_df = sample_df.dropna(thresh=300000)
    sample_df = sample_df.dropna(axis=1, thresh=5)
    sample_df = sample_df.astype(float)

    shape = sample_df.shape
    print('Final number of cases in methylation dataset :', shape[0])

    return(sample_df)


def metaDataFormat(df, all_files):
    """Returns meta data as df. Index column is file_id, metadata = days_to_death/days_in_study, vital_status"""
    # makes list of days to death dataframes matching file_ids with beta arrays passing missing data conditions
    index = list(df.index)
    meta_df_list = []
    for entry in all_files:
        if entry['file_id'] in index:
            # Alive cases handling
            if entry['cases'][0]['demographic']['vital_status'] == "Alive":
                try:
                    meta_df_list.append(
                        {                                
                            'file_id': entry['file_id'],
                            'days_to_death': entry['cases'][0]['diagnoses'][0]['days_to_last_follow_up'],
                            'vital_status': entry['cases'][0]['demographic']['vital_status']
                        }
                    )
                except KeyError:
                    meta_df_list.append(
                        {
                            'file_id': entry['file_id'],
                            'days_to_death': None,
                            'vital_status': entry['cases'][0]['demographic']['vital_status']
                        }
                    )
            # Dead cases handling 
            elif entry['cases'][0]['demographic']['vital_status'] == "Dead":
                try:
                    meta_df_list.append(
                        {
                            'file_id': entry['file_id'],
                            'days_to_death': entry['cases'][0]['demographic']['days_to_death'],
                            'vital_status': entry['cases'][0]['demographic']['vital_status']
                        }
                    )
                except KeyError:
                    meta_df_list.append(
                        {
                            'file_id': entry['file_id'],
                            'days_to_death': None,
                            'vital_status': entry['cases'][0]['demographic']['vital_status']
                        }
                    )

    # merging meta datasets and formatting for use in glm
    meta_df = pd.DataFrame(meta_df_list)
    meta_df = meta_df.set_index(list(meta_df)[0])
    meta_df.index.name = None
    meta_df = meta_df.squeeze()

    # prints final number of sets that pass missing data conditions
    shape = meta_df.shape
    print('Final number of cases in metadata dataset :', shape[0])

    return(meta_df)