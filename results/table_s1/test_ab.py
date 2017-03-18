from ipdb import set_trace
import pandas as pd
from sbie_aging.results import table_s1
from os.path import dirname

a_netdata = dirname(table_s1.__file__)+'/a/network-updated.csv' 
a_processed = dirname(table_s1.__file__)+'/a/network-processed.csv' 
b_logic = dirname(table_s1.__file__)+'/b/logical-network-generated.txt'

def test_a():
    df = pd.read_csv(a_netdata, names=['source','interaction','target'])
    nodes = df['source'].values.tolist() + df['target'].values.tolist()
    nodes_set = set(nodes)
    dictdata = dict([(node, node) for node in nodes_set])
    dictdata = {'4EBP1': 'FOUR_EBP1',
     'AKT': 'AKT',
     'AMP/ATP, ADP/ATP': 'AXP',
     'AMPK': 'AMPK',
     'ATG13': 'ATG13',
     'ATG5': 'ATG5',
     'ATP': 'ATP',
     'BAX': 'BAX',
     'BIM': 'BIM',
     'BMI1': 'BMI1',
     'Bcl-2': 'Bcl2',
     'CDK4,6': 'CDK',
     'DNA damage': 'DNAdamage',
     'E2F': 'E2F',
     'ERK': 'ERK',
     'FOXO': 'FOXO',
     'G6pase': 'G6pase',
     'GLUT4': 'GLUT4',
     'Glucose': 'Glucose',
     'Glycolysis': 'Glycolysis',
     'IGF1': 'IGF1',
     'IGF1R': 'IGF1R',
     'IRS': 'IRS',
     'LC3-2': 'LC',
     'Low nutrient,exercise': 'LowNuEx',
     'MDM2': 'MDM2',
     'MKK3,MKK6': 'MKK',
     'MnSOD': 'MnSOD',
     'NAD+': 'NAD',
     'PDK1,2': 'PDK',
     'PGC1a': 'PGC1a',
     'PI3K': 'PI3K',
     'PP2A': 'PP2A',
     'PTEN': 'PTEN',
     'RAS': 'RAS',
     'ROS': 'ROS',
     'Rheb': 'Rheb',
     'S6K': 'S6K',
     'SGK': 'SGK',
     'SIRT1': 'SIRT1',
     'TCA cycle': 'TCAcyc',
     'TSC1,TSC2': 'TSC',
     'ULK1': 'ULK1',
     'caspase-3': 'caspase3',
     'mTOR': 'mTOR',
     'p16': 'p16',
     'p21': 'p21',
     'p27': 'p27',
     'p38MAPK': 'p38MAPK',
     'p53': 'p53',
     'pRB1': 'pRB1'}
     
    df['source'] = df['source'].apply(lambda x: dictdata[x])
    df['target'] = df['target'].apply(lambda x: dictdata[x])

    # save network with short nodes names 
    df.to_csv(a_processed, index=False)

def test_b():
    df = pd.read_csv(a_processed)
    nodes = set(df['source'].values.tolist() + df['target'].values.tolist())
    # set_trace()
    grp = df.groupby('target')
    with open(b_logic,'w') as fout: 
        for node in nodes: 
            fout.write('%s = Random\n' % node)

        for target, df0 in grp: 
            plus = [] 
            minus = []
            for i in df0.index: 
                src = df0.loc[i, 'source']
                itype = df0.loc[i, 'interaction']
                if itype > 0: 
                    plus.append(src)
                else: 
                    minus.append(src)

            s = target + ' *= '
            if plus != []:
                s += "("+" or ".join(plus)+")"
            if minus != []:
                if plus != []: 
                    s += " and not " + "(" + " or ".join(minus) + ")"
                else: 
                    s += " not " + "(" + " or ".join(minus) + ")"

            fout.write(s+'\n')
            
