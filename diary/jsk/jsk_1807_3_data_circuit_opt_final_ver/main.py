import datetime
import os
import numpy as np
from PyGMO import algorithm, population, island#, util
from optweightandbasal import OptWeightandBasal

#variables relevant to genetic algorithm performance
num_gens = 3500 #number of generations
popsize = 300 #population size
iteration = 300

champfitlist = []#list of best fit parameters' fitness
champparamslist = []#list of best fit parameters
attractorlist = []#list of attractors for best fit parameters
incornodeslist = []#list of incorrectly modeled nodes

for itno in range(iteration):
    prob = OptWeightandBasal()#dim = edgeno + nodeno; weightno + basalno
    algo = algorithm.jde(gen=1)#self adaptive differential evolution algorithm
    #algo = algorithm.vega(gen=1)#multivariable simple genetic algorithm
    pop = population(prob, popsize)
    isl = island(algo, pop)#island class ready to evolve
    
    for i in range(num_gens):
        isl.evolve(1)
        isl.join()#evolve until end of generations
        champfit = isl.population.champion.f[0]
        #champfit = isl.population.champion.f # for multivariable genetic algorithm
        champparams = isl.population.champion.x
        print("[Generation #%d]The current best fitness: %f" % (i+1, champfit))
        if champfit == 0:
            break
        #sum_fit = sum(champfit)
        #print("[Generation #%d]The current best fitness: %s (%d)"%(i+1, champfit, sum_fit))
        #print("- params:", champparams[:3])
        #if sum_fit == 0:
        #    break
    # end of for
    
    #gathering attractors for each conditions on the best fitness parameters
    ws = np.array(champparams[:prob.num_edges - len(prob.forcew)])#weight parameters
    bs = np.array(champparams[prob.num_edges - len(prob.forcew):])#basal activity parameter
    
    for fw in prob.forcew:#fw[0] is the index and fw[1] is the value to append
        ws = np.insert(ws, fw[0], fw[1])
    for fb in prob.forceb:#fb[0] is the index and fb[1] is the value to append
        bs = np.insert(bs, fb[0], fb[1])
    
    attractors = []
    incornodes = []
    for condno in range(prob.num_conds):
        condno += 1
        inistates = prob.conddata.loc['INPUT: ' + str(condno)].tolist()
        for idx in prob.iptidx:#fixing input nodes using extreme basal activity settings
            if inistates[idx] == 0:
                bs[idx] = -prob.wr
            elif inistates[idx] == 1:
                bs[idx] = prob.wr
        inistates = np.array(inistates)
        inistates[inistates == 9] = 0
        attractor = prob._func(ws, bs, inistates)#calculated states form the Boolean model
        attractors.append(attractor)
        objstates = prob.conddata.loc['OUTPUT: ' + str(condno)].tolist()#objective states
        objidx = [idx for idx, state in enumerate(objstates) if state != 9]
        incornode = []
        for attstate in attractor:
            for idx in objidx:
                if objstates[idx] != attstate[idx]:
                    incornode.append(prob.node_list[idx])
            incornode.append('\t')
        incornodes.append(incornode)
    champfitlist.append(champfit)#fitness of best fit parameters
    champparamslist.append(champparams)#best fit parameters
    attractorlist.append(attractors)#attractors for best fit parameters
    incornodeslist.append(incornodes)#incorrect nodes for best fit parameters; none if perfect
    
    #record the fitted weight and basal activity only for the networks with fitness 0
    with open('netparamsdump.txt', 'a') as f:
        if os.stat('netparamsdump.txt').st_size == 0:#if file is empty
            f.write('Format: tab separated\nWeights\nBasal Activities\n\n')
        if champfit == 0:
            ws = np.array(champparams[:prob.num_edges - len(prob.forcew)], dtype = np.float)#weight parameters
            bs = np.array(champparams[prob.num_edges - len(prob.forcew):], dtype = np.float)#basal activity parameter
            for w in ws:
                f.write('%f\t' % w)
            f.write('\n')
            for b in bs:
                f.write('%f\t' % b)
            f.write('\n\n')

algo_name = algo.get_name()
dt = datetime.datetime.now()
str_tp = dt.strftime("%y%m%d_%Hh%Mm%Ss")
fstr_out = 'output_fittedparams_' + prob.netfilename[:-4] + '_' + str_tp + '.txt'
with open(fstr_out, "wt") as f_out:#documenting results
    f_out.write("Algorithm: %s\nNumber of Generations: %d\nPopulation Size: %d\n" % (algo_name, num_gens, popsize))
    wslist = []
    bslist = []
    for itno, (champparams, champfit, attractors, incornodes) in enumerate(zip(champparamslist, champfitlist, attractorlist, incornodeslist)):
        itno += 1
        f_out.write("%d. Best fitness: %f\n" % (itno, champfit))
        ws = np.array(champparams[:prob.num_edges - len(prob.forcew)])#weight parameters
        bs = np.array(champparams[prob.num_edges - len(prob.forcew):])#basal activity parameter
        f_out.write('Weights:\n')
        for w in ws:
            f_out.write('%f\t' % w)
        f_out.write('\nBasal Activities:\n')
        for b in bs:
            f_out.write('%f\t' % b)
        f_out.write('\nAttractors:\n')
        for attractor in attractors:
            f_out.write(str(attractor) + '\n')
        f_out.write('\nIncorrect nodes:\n')
        for incornode in incornodes:
            f_out.write(str(incornode) + '\n')
        f_out.write('\n\n')
        wslist.append(ws)
        bslist.append(bs)
    bestchampfit = min(champfitlist)
    champfitfreq = 100 * champfitlist.count(bestchampfit) / iteration
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

#record the best fitness for this particular pkn
#with open('pknsbestfitness.txt', 'a') as f:
#    f.write(prob.netfilename)
#    f.write(': %f\n' % min(champfitlist))

#if the simulation contained a perfect fitness, add to the list of successful pkns
#if 0 in champfitlist:
#    with open('workinigpkn.txt', 'a') as f:
#        f.write(prob.netfilename + '\n')



