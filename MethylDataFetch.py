import requests
import json
import re
import os
import pandas as pd
from tqdm.auto import tqdm

"""Contains functions related to pulling metadata and downloading methylation beta arrays, running file runs download function"""


def getMethylMetaData(primary_site, file_type="Methylation Beta Value"):
    """Returns list of dictonaries of specfic meta data fields corresponding to input priamry site"""

    # filters to return on files with methylation data from chosen primary site
    filters = {
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "cases.primary_site",
                "value": [str(primary_site)]
                }
            },
            # methylation data filter
            {
            "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["dna methylation"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_type",

                ### "Masked Intensities" selects for .IDAT raw methylation files
                ### "Methylation Beta Values" selects for .txt beta arrays files
                "value": [str(file_type)]
                }
            }
        ]
    }

    # meta data fields requested
    ### could add more data catagories here (bmi, sex... ) https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#file-fields
    fields = [
        "file_id",
        "file_name",
        "file_size",
        "platform",

        "cases.case_id",

        "cases.diagnoses.created_datetime",
        "cases.diagnoses.days_to_last_follow_up",

        "cases.demographic.vital_status",
        "cases.demographic.days_to_birth",
        "cases.demographic.days_to_death",

        "cases.samples.sample_id",
        "cases.samples.created_datetime"
        ]
    fields = ",".join(fields)

    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "JSON",
        "size": "1000"
        }

    response = requests.get("https://api.gdc.cancer.gov/files", params = params)

    data_array = json.loads(response.content.decode("utf-8"))["data"]["hits"]
    return data_array


def checkDownloadSize(primary_site, file_type="Methylation Beta Value"):
    """Prints size of all methylation files in filtered search"""

    size = 0
    counter = 0 
    for file in getMethylMetaData(str(primary_site), str(file_type)):
        size += file["file_size"]
        counter += 1

    print("Number of matching files = " + str(counter))
    print("File download size is " + str((size/1000000)) + " megabytes")
    print(str((size/1000000000)) + " gigabytes")


def downloadMethylFiles(primary_site, file_type="Methylation Beta Value"):
    """Downloads methylation files of requested type and primary site into folder beta_arrays"""

    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        os.mkdir(script_dir + "/beta_arrays")
    except OSError as error:
        print(error)


    checkDownloadSize(primary_site, file_type)
    proceed = input("Would you like to download? (Y/N): ")

    if proceed != "Y":
        return
    else:

        uuid_list = []
        for file in getMethylMetaData(str(primary_site), str(file_type)):
            uuid_list.append(file['id'])

        params = {"ids": uuid_list}

        response = requests.post("https://api.gdc.cancer.gov/data",
                                data = json.dumps(params),
                                headers={
                                    "Content-Type": "application/json"
                                        })

        response_head_cd = response.headers["Content-Disposition"]
        file_name = re.findall("filename=(.+)", response_head_cd)[0]

        save_name = script_dir + "/beta_arrays/" + file_name
        with open(save_name, "wb") as output_file:
            output_file.write(response.content)


def getMethylBetaArrays(primary_site):
    """Returns dataframe with cpg index as index column and each file_id as a column with beta values, 
        method to avoid downloading files onto PC, still takes as long"""

    # Creates directory for temporary file storage
    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        os.mkdir(script_dir + "/data")
    except OSError as error:
        print(error)

    # Get all releavent files and compile into uuid_list
    uuid_list = []
    for file in getMethylMetaData(str(primary_site)):
        uuid_list.append(file['id'])

    df = pd.DataFrame()

    # For each file UUID
    for file_uuid in tqdm(uuid_list):

        # API call
        response = requests.post("https://api.gdc.cancer.gov/data",
                        data = json.dumps({"ids": str(file_uuid)}),
                        headers={
                            "Content-Type": "application/json"
                            })
        
        # File meta data retrieval
        response_head_cd = response.headers["Content-Disposition"]
        file_name = re.findall("filename=(.+)", response_head_cd)[0]
        save_name = script_dir + "/data/" + file_name

        # downloading file
        with open(save_name, "wb") as output_file:
            output_file.write(response.content)

        # reading in file
        read_file = pd.read_table(save_name, header = None)
        read_file.rename(columns = {0:'cpg',1:str(file_uuid)}, inplace = True)
        
        # merging read file into main df
        if df.empty != True:
            df = df.merge(read_file, on='cpg')
        else:
            df = read_file

        # deleting downloaded file
        os.remove(save_name)

    df = df.set_index('cpg')
    
    return(df)


if __name__ == "__main__":
    downloadMethylFiles(input("Enter primary site name (eg. pancreas): "))
