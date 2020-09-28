import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx


'''
Effect of IGF-1's absence on senescence phenotype and positive feedback's contribution
'''

input_nodes = ['DNAdamage', 'IGF1', 'lowNutrition']#100 vs 110; 4-6
positive_feedback_nodes = ['PDK1', 'AKT1', 'IKBKB', 'PTEN']
output_nodes = ['E2F1', 'EIF4EBP1', 'IL1B', 'IL6', 'S6K1', 'ULK1']

attractor_avg_df_list = []
for pdk1_inhibit, pdk1_controlled_network_only in zip([False, True, False, True], [False, False, True, True]):
    for pdk1_force in [(False, False), (False, True), (True, False), (True, True)]:        
        if pdk1_inhibit:
            pdk1_force_str = 'inhibit'
        else:
            pdk1_force_str = 'activate'
        if pdk1_controlled_network_only:
            netparams_str = '_pdk1_network_only'
        else:
            netparams_str = ''        
        #save organized and total result
        attractor_avg_df = pd.read_csv('results/average_attractor_pdk1_{}_bf_{}_af_{}{}.csv'.format(pdk1_force_str, int(pdk1_force[0]), int(pdk1_force[1]), netparams_str), index_col=0)
        attractor_avg_df_list.append(attractor_avg_df)

rep_case = attractor_avg_df_list[0]
igf1_removed = rep_case['4'] - rep_case['6']
act_delta_vec = abs(igf1_removed).sort_values()

#find shorest path to igf-1
network = []
with open('PKN_24.sif', 'r') as f:
    for link in f:
        network.append(link.strip().split('\t'))
node_list = []
for link in network:
    if not link[0] in node_list:
        node_list.append(link[0])
    if not link[2] in node_list:
        node_list.append(link[2])
node_list = sorted(node_list)

G = nx.DiGraph()
G.add_nodes_from(node_list)
for link in network:
    if link[1] == 'activate':
        weight = 1
    else:
        weight = -1
    G.add_edge(link[0], link[2], weight=weight)

shortest_path_igf1 = {}
for node in node_list:
    try:
        shortest_path_igf1[node] = nx.shortest_path_length(G, source='IGF1', target=node)
    except nx.NetworkXNoPath:
        shortest_path_igf1[node] = -1
shortest_path_igf1 = pd.Series(shortest_path_igf1)


colors = ['r' if node in ['PDK1', 'AKT', 'IKBKB', 'PTEN'] else 'b' for node in node_list]
sizes = [100 if node in ['PDK1', 'AKT', 'IKBKB', 'PTEN'] else 10 for node in node_list]

#show results
fig, ax = plt.subplots()
ax.scatter(shortest_path_igf1[node_list], act_delta_vec[node_list], label=node_list, c=colors, s=sizes)
for node in node_list:
    ax.annotate(node, (shortest_path_igf1[node], act_delta_vec[node]))



















