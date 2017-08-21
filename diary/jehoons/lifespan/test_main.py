import pandas as pd 

def test():
    df_fromto = pd.read_csv(fromto,sep=r'\s*')
    
    df_lifespan = pd.read_csv(lifespan)
    
    df_network = pd.read_csv(network,sep=r'\t')

    df_network = df_network[ df_network['combined_score'] > 0.9] 

    df2 = df_lifespan[['lifespan', 'uniprot_id']].merge(df_fromto,
        how='left', left_on='uniprot_id', right_on='From').drop_duplicates()

    df2 = df2.dropna()

    df2[['To','lifespan']].to_csv(propfile, index=False)

    df_network2 = df_network.merge(df2[['To', 'lifespan']], how='left', left_on='#node1', right_on='To')

    df_network3 = df_network2.merge(df2[['To', 'lifespan']], how='left', left_on='node2', right_on='To')

    df_network3.drop('To_x', 1, inplace=True)
    
    df_network3.drop('To_y', 1, inplace=True)

    df_network3['edgetype'] = df_network3['lifespan_x'] + ';' +df_network3['lifespan_y']

    df_network3.to_csv(network_updated)


# database 
lifespan = 'lifespanofmice_v3.csv'

# id mapping 
fromto = 'fromto.txt' 

# string-db 
network = 'network.txt'

# lifespan +, - ?
propfile = 'node-prop.csv'

network_updated = 'network-updated.txt'
