import pandas as pd
df = pd.read_csv('agingrppa.csv', index_col=0)
col_headers = list(df.columns)
row_headers = list(df.index)
# print(len(col_headers),len(row_headers))
# q=df[1:4]
# print(col_headers)
q = df[['Quiescent 1h/control', 'Quiescent 6h/control', 'Quiescent 24h/control', 'Quiescent 96h/control']]
rapa = df[['Rapamycin 1h/control', 'Rapamycin 6h/control', 'Rapamycin 24h/control', 'Rapamycin 96h/control']]
igf = df[['IGF 0.5h/control', 'IGF 1h/control', 'IGF 3h/control', 'IGF 6h/control', 'IGF 24h/control']]
apop = df[['Apoptosis 0.5h/control', 'Apoptosis 1h/control', 'Apoptosis 3h/control', 'Apoptosis 6h/control']]
s = df[['Senescent 1h/control', 'Senescent 6h/control', 'Senescent 24h/control', 'Senescent 96h/control']]
# print(row_headers)
list_q =[]
list_rapa =[]
list_igf =[]
list_apop =[]
list_s=[]
for i in row_headers:
    zscores = list(q.ix[i])
    # print(zscores)
    if len([x for x in zscores if abs(x)>=2]) >=2:
        list_q.append(i)
    zscores = list(rapa.ix[i])
    if len([x for x in zscores if abs(x) >= 2]) >= 2:
        list_rapa.append(i)
    zscores = list(igf.ix[i])
    if len([x for x in zscores if abs(x) >= 2]) >= 2:
        list_igf.append(i)
    zscores = list(apop.ix[i])
    if len([x for x in zscores if abs(x) >= 2]) >= 2:
        list_apop.append(i)
    zscores = list(s.ix[i])
    if len([x for x in zscores if abs(x) >= 2]) >= 2:
        list_s.append(i)
dic_node = {}
dic_node['q'] = list_q
dic_node['rapa'] = list_rapa
dic_node['igf'] = list_igf
dic_node['apop'] = list_apop
dic_node['s'] = list_s
f = open('essential_node_rppa.txt','w')
f.write('condition\tnode\n')
for key in dic_node:
    f.write('%s\t%s\n'%(key,'\t'.join([x for x in dic_node[key]])))
# print(len(set(list_igf).intersection(set(list_s))))
# print(len(list_s),len(list_apop),len(list_q))
#
# print(len(set(list_apop).intersection(set(list_s).intersection(set(list_q)))))
# print(len(set(list_apop).intersection(set(list_s))))
# # print(len(set(list_rapa).intersection(set(list_s))))
# print(len(set(list_q).intersection(set(list_s))))
# print(len(set(list_q).intersection(set(list_apop))))
# print(len(set(list_q).intersection(set(list_s))))
set_total = set(list_apop + list_q+list_rapa+list_s+list_igf)
node_list = []
for i in set_total:
    new = i[:i.find(' (')]
    node_list.append(new)
# for i in node_list:
#     print(i)

list_jh =['fcgr2a', 'map2k1', 'tank', 'ptk2b', 'jund', 'met','cxcr4', 'hdac1', 'eif4ebp1', 'gsk3a', 'ppp2r4', 'slc6a13'
    ,'casp9', 'elk1', 'ankrd23', 'hnf4a', 'myc']
list_jh_cap = []
for i in list_jh:
    list_jh_cap.append(str(i.upper()))


list_sg = ['FKBP4','BLNK','DDX5','CCND1','MAPK14','FLT3','PLCG1','ELK1','RPS6KA5','ITGB3','ICK','STAT3','BCL2A1','TP53'
    ,'AKT1','FCGR2A','CALM2','ASRGL1','IL7R','ATRIP','A2M','PTK2B','LIMK1','EIF4EBP1','CDH5','NF2','PLCG2','PEA15'
    ,'PCMT1','RGS16','IKBKG','GLMN','CHUK','ABL1','HDAC1','SNRPE','LYN','RPS6KA1','APAF1']
print(len(set(list_jh_cap)), len(set(list_sg)),set(list_jh_cap).intersection(set(list_sg)))
print(len(set(list_jh_cap)), len(set(list_s)),len(set(list_jh_cap).intersection(set(list_s))))