import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from scipy.stats import rankdata


#variables
antibody_list = ['4E-BP1 (Ab-65)', 'AKT1 (Phospho-Ser473)', 'E2F1 (Phospho-Thr433)', 
                 'IGF1R (Phospho-Tyr1165/1166)', 'IkB-alpha (Ab-32/36)', 
                 'IKK-beta (Ab-188)', 'IRS-1 (Phospho-Ser636)', 
                 'MKK6/MAP2K6 (Phospho-Ser207)', 'mTOR (Phospho-Ser2481)', 
                 'p38 MAPK (Phospho-Tyr322)', 'p53 (Phospho-Ser15)', 
                 'P70S6K (Phospho-Ser411)', 'Rb (Phospho-Thr821)', 
                 'Tuberin/TSC2 (Phospho-Ser939)']
#timepoints_list = [np.array([0, .5, 1, 3, 6, 24]) / 24, 
timepoints_list = [np.array([0, 1, 6, 24]) / 24, 
                  np.array([0, 1, 6, 24, 96]) / 96, 
                  np.array([0, 1, 6, 24, 96]) / 96]
condition_list = ['Proliferation', 'Quiescence', 'Senescence']
z_score = True

#import data
rppa_solution = pd.read_csv('rppa_auc_solution.csv', index_col=0)
rppa_datas = []
for i in range(6):
    _rppa_data = pd.read_csv('rppa_fc_data_{}.csv'.format(i), index_col=0)
    rppa_datas.append(_rppa_data)

rppa_auc_df_list = []
for _, rppa_data in enumerate(rppa_datas):
    rppa_data.index = [_ab.strip() for _ab in rppa_data.index]
    original = False
    if _ == 0 or _ == 1:
        original = True
    rppa_auc_df = pd.DataFrame(index=antibody_list, columns=condition_list, dtype=float)
    for antibody in antibody_list:
        antibody_data = rppa_data.loc[antibody]
        for timepoints, condition in zip(timepoints_list, condition_list):
            datapoints = np.concatenate(([1], antibody_data.loc[[condition in _ for _ in antibody_data.index]].values))
            if original and condition == 'Proliferation':
                datapoints = np.delete(datapoints, 3)
                datapoints = np.delete(datapoints, 1)
            rppa_auc_df[condition][antibody] = auc(timepoints, datapoints)
    rppa_auc_df = rppa_auc_df[['Quiescence', 'Senescence', 'Proliferation']]
    if z_score:
        for ab, rowdata in rppa_auc_df.iterrows():
            rppa_auc_df.loc[ab] = (rowdata - rowdata.mean())/rowdata.std()
    rppa_auc_df_list.append(rppa_auc_df)

rppa_auc_df_mean = pd.concat(rppa_auc_df_list).groupby(level=0).mean()
rppa_auc_df_mean = rppa_auc_df_mean.loc[antibody_list]

rppa_auc_df_median = pd.concat(rppa_auc_df_list).groupby(level=0).median()
rppa_auc_df_median = rppa_auc_df_median.loc[antibody_list]



#rank solution
rppa_solution_ranked = pd.DataFrame(np.zeros_like(rppa_solution), index=rppa_solution.index, columns=rppa_solution.columns, dtype=int)
for _idx, _rowdata in rppa_solution.iterrows():
    rppa_solution_ranked.loc[_idx] = rankdata(_rowdata, method='max')

rppa_auc_names = ['data_0', 'data_1', 'data_2','data_3', 'data_4', 'data_5', 'data_mean', 'data_median']
for _rppa_auc_name, _rppa_auc_data in zip(rppa_auc_names, rppa_auc_df_list + [rppa_auc_df_mean, rppa_auc_df_median]):
    #rank data
    _ranked = pd.DataFrame(np.zeros_like(_rppa_auc_data), index=_rppa_auc_data.index, columns=_rppa_auc_data.columns, dtype=float)
    for _idx, _rowdata in _rppa_auc_data.iterrows():
        _ranked.loc[_idx] = rankdata(_rowdata, method='max')
    
    #calculate congruency
    _congruency = np.sum(2**(_ranked - rppa_solution_ranked), axis=1) == 2.5
    _rppa_auc_data['Congruency'] = _congruency
    _rppa_auc_data.to_csv('results/{}.csv'.format(_rppa_auc_name))































