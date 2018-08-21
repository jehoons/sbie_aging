import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#parameters for simulation
mode = 1#0: mean of fitness, 1: mode of fitness under threshold
thres = .5#fitness cutoff; only needed if the mode == 1
boolrecinit = 500#iteration number at which to start recording state history
booltotit = 1000#total iteration of boolean simulation
dosagebinsize = .01#bin size of drug dosage
#nodes/control targets to inhibit/activate; upto 2
#[['nodename1', 0(inhibition) or 1(activation), default_dosage=0],
# ['nodename2', 0(inhibition) or 1(activation), default_dosage=0]]
inhibparams = [['PDK1', 0, 0]]

#seed the random function for reproducability
np.random.seed(1)

#import network information
netfilename = "PKN_22.sif"
netdata = []
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

#import pd.DataFrame with input and output condition information
datafilename = 'obj_attractor_sq.csv'
conddata = pd.read_csv(datafilename, index_col = 0)

#record positions and values of input condition
iptidx = []#indicies of the input conditions
iptval = []#input value for code revision
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)
        if state == 1:
            iptval.append(1)
        else:
            iptval.append(0)

#record values of objective conditions
optval = []
for idx, state in enumerate(list(conddata.iloc[1])):
    if state != 9:
        if state == 1:
            optval.append(1)
        else:
            optval.append(0)
optval = np.array(optval)

#basic variables
node_list = list(conddata.columns)
num_nodes = len(node_list)
num_edges = len(netdata)

#import fitted model parameters calculated from previous simulations
netparamslist = []
paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

#the data is in a specific format which the code under deconstructs
#redistribute parameters to weights and basal activities
totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.float64)
bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.float64)
for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)

#indexing input basal activities as it does not serve a purpose; unnecessary
tempbs = np.transpose(bs).copy()
bs = np.zeros((len(node_list), totdatano), dtype=np.float64)
i = 0
for idx in range(len(node_list)):
    if idx in iptidx:
        bs[idx] = np.array([iptval[iptidx.index(idx)]] * totdatano)
    else:
        bs[idx] = tempbs[i]
        i += 1
bs = np.transpose(bs)

#repackage weights and basal activities parameters by models
for w, b in zip(ws, bs):
    netparamslist.append((w, b))


###the Boolean network model
def _func(w, b, inistates, inhibparams):
    #generate weight matrix from the global variable; not the best code
    weight_matrix = np.zeros((num_nodes, num_nodes))
    weight_matrix[ind_rows, ind_cols] = signs*np.round(w)
    
    #record the state transition history
    state_history = [] # Record of state transition
    state_history = np.zeros((booltotit + 1, len(node_list)))
    
    x_t1 = np.array(inistates)
    state_history[0] = x_t1#record the first state
    for it in range(1, booltotit + 1):#state transition from t to t + 1
        x_t2 = np.dot(x_t1, weight_matrix) + b #transitioned state
        
        #forcing the initial states
        for idx in iptidx:             
            if inistates[idx] == 0:
                x_t2[idx] = 0
            elif inistates[idx] == 1:
                x_t2[idx] = 1

        #forcing control target states with appropriate probability
        for inhibparam in inhibparams:
            if np.random.rand() < inhibparam[2]:
                x_t2[node_list.index(inhibparam[0])] = int(inhibparam[1] > 0)

        #activation function
        x_t2[x_t2 > 0] = 1
        x_t2[x_t2 <= 0] = 0
        
        #move on
        state_history[it] = x_t2.copy()
        x_t1 = x_t2
    return state_history[boolrecinit + 1:]
    
#generate initial states and objective states from the imported data
inistates = np.array(conddata.loc['INPUT'].tolist())#initial states
inistates[inistates == 9] = 0#the unspecified initial states are set to 0
objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
objboolmask = objstates != 9#objective states mask for only phenotypic nodes
objphen = objstates[objboolmask]
if len(inhibparams) <= 1:
    drugeff1d = []
    for dose in np.arange(0, 1 + dosagebinsize, dosagebinsize):
        inhibparams[0][2] = dose#setting dose for the drug
        drugeffbymodel = []#vector with mean node activation by models
        for idx, (w, b) in enumerate(netparamslist):#simulate each model
            signs = np.ones((num_edges,), dtype=np.int32)
            ind_rows = np.zeros((num_edges,), dtype=np.int32)
            ind_cols = np.zeros((num_edges,), dtype=np.int32)
            
            #setting up rows and columns vector for weight matrix
            for idx, edge in enumerate(netdata):
                if edge[1] == 'inhibit':
                    signs[idx] = -1
                ind_rows[idx] = node_list.index(edge[0])
                ind_cols[idx] = node_list.index(edge[2])
            
            #calculate state history of phenotype nodes for this model
            state_history = _func(w, b, inistates, inhibparams)
            phen_hist = state_history[:, objboolmask]
            
            #calculate effectiveness of control by measuring fitness
            #fitness is the euclidian distance b/w model and desired state
            #the drug is considered effective if the fitness is close to 0
            meanact = np.mean(phen_hist, axis = 0)
            fitness = np.sum((optval - meanact)**2)
            if mode > 0 and fitness > thres:
                continue
            drugeffbymodel.append(fitness)
        if mode <= 0:
            drugeff = np.mean(np.array(drugeffbymodel))
        else:
            drugeff = len(drugeffbymodel)
        drugeff1d.append(drugeff)
    drugeffarr = np.array(drugeff1d)
    
    #plot histogram of drug effectiveness
    fig = plt.figure()
    plt.bar(np.arange(len(drugeffarr)) * dosagebinsize, drugeffarr,
            dosagebinsize)
    if inhibparams[0][1] > 0:
        xlabel = inhibparams[0][0] + ' activation'
    else:
        xlabel = inhibparams[0][0] + ' inhibition'
    plt.xlabel(xlabel, fontsize=18)
    title = 'Senescence to Quiescence Control'
    #title = 'Senescence Escape'
    if mode <= 0:
        title = title + ': Mean Fitness'
        plt.ylabel('Mean fitness', fontsize=18)
    else:
        title = title + ': cutoff ' + str(thres)
        plt.ylabel('Mode of fitness less than ' + str(thres), fontsize=18)
        if thres < 1:
            plt.ylabel('Successfully controlled models', fontsize=18)
    if datafilename[-5] == 'q':
        title = title + ' (Quiescence phenotype)'
    elif datafilename[-5] == 'p':
        title = title + ' (Proliferation phenotype)'
    else:#for senescence phenotype
        title = title + ' (Senescence phenotype)'
    plt.title(title, fontsize=20)
    plt.tick_params(labelsize=14)
    plt.ylim(0, len(netparamslist))
    plt.show()
else:
    drugeff2d = []
    for dose1 in np.arange(0, 1 + dosagebinsize, dosagebinsize):
        inhibparams[0][2] = dose1#setting dose for drug1
        drugeff1d = []
        for dose2 in np.arange(0, 1 + dosagebinsize, dosagebinsize):
            inhibparams[1][2] = dose2#setting dose for drug2
            drugeffbymodel = []#vector with mean node activation by models
            for idx, (w, b) in enumerate(netparamslist):#simulate each model
                signs = np.ones((num_edges,), dtype=np.int32)
                ind_rows = np.zeros((num_edges,), dtype=np.int32)
                ind_cols = np.zeros((num_edges,), dtype=np.int32)
                
                #setting up rows and columns vector for weight matrix
                for idx, edge in enumerate(netdata):
                    if edge[1] == 'inhibit':
                        signs[idx] = -1
                    ind_rows[idx] = node_list.index(edge[0])
                    ind_cols[idx] = node_list.index(edge[2])
                
                #calculate state history of phenotype nodes for this model
                state_history = _func(w, b, inistates, inhibparams)
                phen_hist = state_history[:, objboolmask]
                
                #calculate effectiveness of control by measuring fitness
                #fitness is the euclidian distance b/w model and desired state
                #the drug is considered effective if the fitness is close to 0
                meanact = np.mean(phen_hist, axis = 0)
                fitness = np.sum((optval - meanact)**2)
                if mode > 0 and fitness > thres:
                    continue
                drugeffbymodel.append(fitness)
            if mode <= 0:
                drugeff = np.mean(np.array(drugeffbymodel))
            else:
                drugeff = len(drugeffbymodel)
            drugeff1d.append(drugeff)
        drugeff2d.append(drugeff1d)
    drugeffarr = np.array(drugeff2d)
    
    #plot 3D histogram of drug effectiveness by dosage combinations
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x_data, y_data = np.meshgrid(np.arange(drugeffarr.shape[1]),
                                 np.arange(drugeffarr.shape[0]))
    x_data = x_data.flatten() * dosagebinsize
    y_data = y_data.flatten() * dosagebinsize
    z_data = drugeffarr.flatten()
    ax.bar3d(x_data, y_data, np.zeros(len(z_data)),
             dosagebinsize, dosagebinsize, z_data)
    if inhibparams[0][1] > 0:
        xlabel = inhibparams[0][0] + ' activation'
    else:
        xlabel = inhibparams[0][0] + ' inhibition'
    ax.set_xlabel(xlabel)
    if inhibparams[1][1] > 0:
        ylabel = inhibparams[1][0] + ' activation'
    else:
        ylabel = inhibparams[1][0] + ' inhibition'
    ax.set_ylabel(ylabel)
    title = 'Senescence to Quiescence Control'
    if mode <= 0:
        title = title + ': Mean Fitness'
        ax.set_zlabel('Mean fitness')
    else:
        title = title + ': cutoff ' + str(thres)
        ax.set_zlabel('Mode of fitness less than ' + str(thres))
    if datafilename[-5] == 'q':
        title = title + ' (Quiescence phenotype)'
    elif datafilename[-5] == 'p':
        title = title + ' (Proliferation phenotype)'
    else:#for senescence phenotype
        title = title + ' (Senescence phenotype)'
    plt.title(title)
    plt.show()














