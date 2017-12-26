import numpy as np
import pandas as pd
import time
from PyGMO import algorithm, population, island
from optweightandbasal import OptWeightandBasal

#variables relevant to genetic algorithm performance
n = 50#number of GA simulations
num_gens = 10000#number of generations
popsize = 500#population size

#import network and objective network conditions with inputs
conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
netdata = []#list with network information
netfile = open('PKN_19.sif', 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))
dim = len(netdata) + len(conddata.columns)#genetic algorithm parameter dimension; number of edges(weights) + number of nodes(basal activities)

champfitlist = []#list of best fit parameters' fitness
champparamslist = []#lsit of best fit parameters
for it in range(n):#iterate GA simulation
    #implementation of genetic algorithm
    prob = OptWeightandBasal(dim = dim, i_dim = dim, netdata = netdata, conddata = conddata)#dim = edgeno + nodeno; weightno + basalno
    algo = algorithm.jde(gen=num_gens)#self adaptive differential evolution algorithm
    pop = population(prob, popsize)
    isl = island(algo, pop)#island class ready to evolve
    isl.evolve(1)
    isl.join()#evolve until end of generations
    champfitlist.append(isl.population.champion.f[0])#fitness of best fit parameters
    champparamslist.append(isl.population.champion.x)#best fit parameters

algo_name = algo.get_name()
now = time.localtime()
fstr_out = "output_%s_%d%02d%02d_%02dh%02dm%02ds.txt" % ('fittedparams', now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
with open(fstr_out, "wt") as f_out:#documenting results
    f_out.write("Algorithm: %s\nNumber of Generations: %d\nPopulation Size: %d\nNumber of Iterations: %d\n\n" % (algo_name, num_gens, popsize, n))
    for itno, (champparams, champfit) in enumerate(zip(champparamslist, champfitlist)):
        itno += 1
        f_out.write("%d. Best fitness: %f\n" % (itno, champfit))
        for param in champparams:
            f_out.write("%f " % param)
        f_out.write("\n\n")
