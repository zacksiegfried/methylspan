import pandas as pd
import os
from MethylDataFetch import getMethylBetaArrays
from MethylMapping import methylMap450
from MethylMapping import averageGene


def compareGeneMethylContent(list_of_primary_sites):
    '''Returns dataframe'''

    script_dir = os.path.dirname(os.path.abspath(__file__))

    try:
        os.mkdir(script_dir + "/analysis")
    except OSError as error:
        print(error)

    # populating list of dfs
    list_of_frames = []
    for site in list_of_primary_sites:
        list_of_frames.append(getMethylBetaArrays(str(site)))

    # building complete df
    df = pd.DataFrame()
    for i in list_of_frames:
        
        ### Uses array 450 for all files (could improve in MethylMapping.py)
        mapped_frame = methylMap450(i['data'].transpose())

        # gene averaing method
        av_mapped_frame = averageGene(mapped_frame)

        # setup formatting
        av_mapped_frame = av_mapped_frame.transpose()
        av_mapped_frame.index.name = None
        av_mapped_frame.columns.names = ['gene']

        # taking mean of all patients for each gene + dropping all other columns
        av_mapped_frame[str(i['primary_site'])] = av_mapped_frame.mean(axis=1, skipna=True, numeric_only=True)
        av_mapped_frame = av_mapped_frame.reset_index()
        frame_done = av_mapped_frame[['index', str(i['primary_site'])]]

        # merging read file into main df
        if df.empty != True:
            df = df.merge(frame_done, on='index')
        else:
            df = frame_done

    df = df.set_index('index')
    df.index.name = None

    df.to_csv(str(script_dir + '/analysis/comp.csv')) # overwrites itself (could fix)

    return(df)

 

