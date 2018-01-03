import numpy as np
import pandas as pd
import time
from PyGMO import algorithm, population, island, util
from optweightandbasal import OptWeightandBasal

#variables relevant to genetic algorithm performance
totnetno = 8192#number of networks; this should be consistent with actual number of PKN_20 variants
n = 20#number of GA simulations
num_gens = 10000#number of generations
popsize = 500#population size

fitpkn20s = []#list of pkn20 networks with perfect fitness values
for netno in range(totnetno):
    if np.random.uniform() > .01:
        continue#only perform part of the simulation
    
    netno += 1
    print('%d of %d networks' % (netno, totnetno))
    netfilename = 'PKN_20_' + str(netno) + '.sif'
    
    #import network and objective network conditions with inputs
    conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
    netdata = []#list with network information
    netfile = open(netfilename, 'r')
    for link in netfile:
        netdata.append(link.strip().split('\t'))
    dim = len(netdata) + len(conddata.columns)#genetic algorithm parameter dimension; number of edges(weights) + number of nodes(basal activities)
    
    champfitlist = []#list of best fit parameters' fitness
    champparamslist = []#list of best fit parameters
    for it in range(n):#iterate GA simulation
        #implementation of genetic algorithm
        print('%d of %d iterations' % (it + 1, n))
        prob = OptWeightandBasal(dim = dim, i_dim = dim, netdata = netdata, conddata = conddata)#dim = edgeno + nodeno; weightno + basalno
        algo = algorithm.jde(gen=num_gens)#self adaptive differential evolution algorithm
        pop = population(prob, popsize)
        isl = island(algo, pop)#island class ready to evolve
        isl.evolve(1)
        isl.join()#evolve until end of generations
        champfitlist.append(isl.population.champion.f[0])#fitness of best fit parameters
        champparamslist.append(isl.population.champion.x)#best fit parameters
        
        #experimental codes for trying out the analysis procedure in the package
        #inspector = util.analysis(pop, output_to_file=False)
        #inspector.f_distribution()
        #inspector.f_linearilty_convexity()
        #inspector.f_regression(degree=[1,1,2,2,3,3],interaction=[False,True,False,True,False,True])
    
    algo_name = algo.get_name()
    now = time.localtime()
    #fstr_out = "output_%s_%d%02d%02d_%02dh%02dm%02ds.txt" % ('fittedparams', now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
    fstr_out = 'output_fittedparams_' + netfilename[:-4] + '.txt'
    with open(fstr_out, "wt") as f_out:#documenting results
        f_out.write("Algorithm: %s\nNumber of Generations: %d\nPopulation Size: %d\nNumber of Iterations: %d\n\n" % (algo_name, num_gens, popsize, n))
        wslist = []
        bslist = []
        for itno, (champparams, champfit) in enumerate(zip(champparamslist, champfitlist)):
            itno += 1
            f_out.write("%d. Best fitness: %f\n" % (itno, champfit))
            ws = np.array(champparams[:prob.num_edges])#list of weights for each links
            ws = ((ws + 1) / 2).astype(np.int32)
            bs = np.array(champparams[prob.num_edges:])#list of basal activities for each nodes
            bmask = [1 if s % 2 == 1 else -1 for s in bs]
            bs = ((bs + 1) / 2).astype(np.int32)
            bs = np.multiply(bs, bmask)
            f_out.write('Weights:\n')
            for w in ws:
                f_out.write('%d\t' % int(w))
            f_out.write('\nBasal Activities:\n')
            for b in bs:
                f_out.write('%d\t' % int(b))
            f_out.write('\n\n')
            wslist.append(ws)
            bslist.append(bs)
        bestchampfit = min(champfitlist)
        champfitfreq = champfitlist.count(bestchampfit) / n
        f_out.write('Best Best Fitness: %f\n' % bestchampfit)
        f_out.write('Best Fitness Frequency: %f%%\n' % champfitfreq)
        mws = np.mean(wslist, axis = 0)
        mbs = np.mean(bslist, axis = 0)
        f_out.write('Mean Weights:\n')
        for mw in mws:
            f_out.write('%f\t' % mw)
        f_out.write('\nMean Basal Activities:\n')
        for mb in mbs:
            f_out.write('%f\t' % mb)
        f_out.write('\n')
    
    if bestchampfit == 0.0:
        fitpkn20s.append(netfilename)

with open('PKN_20_with_perfect_fitness.txt', "wt") as f_out:#export list of pkn20s with perfect fitness values
    for fitpkn20 in fitpkn20s:
        f_out.write(fitpkn20 + '\n')
