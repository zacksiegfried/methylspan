import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from lifelines import CoxPHFitter
import numpy as np


def survivalPCAFormat(sample_data, meta_data, no_PCA):
    """Takes methylation beta array and performs PCA combines meta data from MSM package + formats"""

    # drops all CpG sites with at least 1 missing value
    sample_data.index.name = 'id'
    sample_data.dropna(axis=1, inplace=True)

    # standardizes dataframe
    scaled_df = sample_data.copy()
    scaled_df = pd.DataFrame(StandardScaler().fit_transform(scaled_df), index=scaled_df.index, columns=scaled_df.columns)

    print("Number of CpG sites with complete data from all cases: " + str(scaled_df.shape[1]))

    # converts CpG sites to PCA components
    pca = PCA(n_components=no_PCA)      # setting number of PCAs to 15 (decided using PCA notebook's scree plot)
    pca.fit(scaled_df)

    output = pca.transform(scaled_df)
    x = pd.DataFrame(output, index=scaled_df.index)

    meta_data.index.name = 'id'

    # merges dfs and formats 
    full_data = pd.merge(meta_data, x, on = 'id')
    full_data.replace(['Dead', 'Alive'],[1, 0], inplace=True)
    full_data = full_data[full_data['days_to_death'].notna()]

    print("Final number of cases used: " + str(full_data.shape[0]))

    return(full_data)


def k_fold_cross_survival(k, data):
    """Performs cross validation for the mean difference between predicted and observed time to event"""
    
    scores = []
    for i in range(k):

        # splitting data set
        X_train = data.sample(frac = 0.7)
        X_test = data.drop(X_train.index)
        X_test_a = X_test.copy(deep=True)

        # fitting Cox model using X_train
        cph = CoxPHFitter()
        cph.fit(X_train, duration_col='days_to_death', event_col='vital_status')

        # predicting values using X_test, formatting
        median_table = cph.predict_median(X_test)
        median_table = median_table.to_frame()
        median_table.columns = ['pred']
        median_table.index.name = 'id'

        # merges predicted and actual values
        pred_data = pd.merge(median_table, X_test_a, on='id')
        # calculates difference between predicted and actual 
        pred_data['difference'] = pred_data['pred'] - pred_data['days_to_death']
        pred_data['difference'].replace([np.inf, -np.inf], 0, inplace=True)

        pred_data2 = pred_data.copy(deep=True)

        ### graphing could be improved (red = death, blue = alive)
        pred_data2.reset_index().plot(kind='scatter', x='id', y='difference', c=np.where(pred_data2['vital_status'], '#fc4f30', '#008fd5')).tick_params(bottom=False, labelbottom=False)

        scores.append(abs(pred_data['difference']).median())
    
    return(scores)