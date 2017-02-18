import pandas as pd 
from ipdb import set_trace
from os.path import dirname,join
import os 
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.style.use('ggplot')

def run_postproc(datafile, processed):
    df = pd.read_csv(datafile, sep='\t')

    output = [i for i in df.columns if i[0:5] == 'Prob[']
    output_2 = [i.replace('Prob[','').replace(']','').split('--') for i in output]

    output_3 = [] 
    for i in output_2: 
        output_3 += i 

    output_3 = set(output_3)

    # node names 
    cols = {}
    for c in output_3:
        cols[c] = [i for i in output if i.find(c) >= 0]

    unnamed_cols = [i for i in df.columns if i[0:7] == 'Unnamed']
    for unnamed in unnamed_cols: 
        df.drop(unnamed, axis=1, inplace=True)
        
    for c in cols: 
        df[c] = df[cols[c]].sum(axis=1)

    df.to_csv(processed, index=False)


if __name__ == '__main__':

    # MaBoSS의 결과는 state vector 기준으로 결과가 나온다. 그러므로 이것을 
    # 개별 state value로 전환해 줄 필요가 있다. 이것을 실행하기 전에 아래 명령을 
    # 먼저 실행해 주어야 한다.
    # > MBSS_FormatTable.pl Loic2016-model.bnd Loic2016-model.cfg

    import sys
    datadir = sys.argv[1]

    datafile = '%s/%s_probtraj_table.csv' % (datadir, datadir)
    processed = '%s/%s_probtraj_table_processed.csv' % (datadir, datadir)

    run_postproc(datafile=datafile, processed=processed)

    df = pd.read_csv(processed) 
    df2 = df[['Time', 'p21', 'mTORC1_S6K1', 'AKT', 'Insulin', 'Senescence', '<nil>']]
    df2.set_index('Time', inplace=True)
    
    plt.figure(); 
    df2.plot()
    # plt.show()
    plt.savefig(datadir+'_fig.jpg')
    # set_trace()

