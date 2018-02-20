import numpy as np
from PyGMO import algorithm, population, island#, util
from optweightandbasal import OptWeightandBasal
from randalgo import RandAlgo

#variables relevant to genetic algorithm performance
n = 1#number of GA simulations
num_gens = 300#number of generations
popsize = 1000#population size

rec = open('runrec.txt', 'w')

champfitlist = []#list of best fit parameters' fitness
champparamslist = []#list of best fit parameters
attractorlist = []#list of attractors for best fit parameters
incornodeslist = []
for it in range(n):#iterate GA simulation
    #implementation of genetic algorithm
    rec.write('%d of %d iterations' % (it + 1, n))
    rec.write('import and delcare OptWeightandBasal class')#
    rec.write('\n')#
    prob = OptWeightandBasal()#dim = edgeno + nodeno; weightno + basalno
    #algo = algorithm.jde(gen=num_gens)#self adaptive differential evolution algorithm
    rec.write('import and declare ended')#
    rec.write('\n')#
    rec.write('import and declare RandAlgo class')#
    rec.write('\n')#
    algo = RandAlgo(num_gens)
    rec.write('import and declare ended')#
    rec.write('\n')#
    rec.write('GA evolution commence')#
    rec.write('\n')#
    pop = population(prob, popsize)
    isl = island(algo, pop)#island class ready to evolve
    isl.evolve(1)
    isl.join()#evolve until end of generations
    rec.write('GA evolution ended')#
    rec.write('\n')#
    champfit = isl.population.champion.f[0]
    rec.write('champfit:')#
    rec.write(str(champfit))#
    rec.write('\n')#
    champparams = isl.population.champion.x
    rec.write('champparams')#
    rec.write(str(champparams))#
    rec.write('\n')#

    #gathering attractors for each conditions on the best fitness parameters
    ws = np.array(champparams[:prob.num_edges - len(prob.forcew)])#weight parameters
    ws = ((ws + 1) / 2).astype(np.int32)#adjusting the boundaries to suit the weight parameters
    rec.write('weights:')#
    rec.write(str(ws))#
    rec.write('\n')#
    bs = np.array(champparams[prob.num_edges - len(prob.forcew):])#basal activity parameter
    bs = np.array([a - prob.wr / 2 if a > prob.wr / 2 else a - prob.wr / 2 - 1 for a in bs]).astype(np.int32)#adjusting the boundaries to cover negative basal activity levels
    rec.write('basal activities:')#
    rec.write(str(bs))#
    rec.write('\n')#
    
    for fw in prob.forcew:#fw[0] is the index and fw[1] is the value to append
        ws = np.insert(ws, fw[0], fw[1])
    rec.write('weights with forced:')#
    rec.write(str(ws))#
    rec.write('\n')#
    for fb in prob.forceb:#fb[0] is the index and fb[1] is the value to append
        bs = np.insert(bs, fb[0], fb[1])
    rec.write(str(prob.node_list))#
    rec.write('basal activities with forced:')#
    rec.write(str(bs))#
    rec.write('\n')#
    
    attractors = []
    incornodes = []
    for condno in range(prob.num_conds):
        condno += 1
        rec.write('condition number:')#
        rec.write(str(condno))#
        rec.write('\n')#
        inistates = prob.conddata.loc['INPUT: ' + str(condno)].tolist()
        for idx in prob.iptidx:#fixing input nodes using extreme basal activity settings
            if inistates[idx] == 0:
                bs[idx] = -prob.maxindeg
            elif inistates[idx] == 1:
                bs[idx] = prob.maxindeg
        rec.write('initial state:')#
        rec.write(str(inistates))#
        rec.write('\n')#
        inistates = list(map(lambda x: int(str(x).replace('9', '0')), inistates))
        rec.write('initial state with nine replaced:')#
        rec.write(str(inistates))#
        rec.write('\n')#
        attractor = prob._func(ws, bs, inistates)#calculated states form the Boolean model
        rec.write('attractor:')#
        rec.write(str(attractor))#
        rec.write('\n')#
        attractors.append(attractor)
        objstates = prob.conddata.loc['OUTPUT: ' + str(condno)].tolist()#objective states
        rec.write('objective states:')#
        rec.write(str(objstates))#
        rec.write('\n')#
        objidx = [idx for idx, state in enumerate(objstates) if state != 9]
        rec.write('objective state indices:')#
        rec.write(str(objidx))#
        rec.write('\n')#
        incornode = []
        for idx in objidx:
            if objstates[idx] != attractor[0][idx]:
                incornode.append(prob.node_list[idx])
        rec.write('incorrect nodes:')#
        rec.write(str(incornode))#
        rec.write('\n')#
        incornodes.append(incornode)
    champfitlist.append(champfit)#fitness of best fit parameters
    champparamslist.append(champparams)#best fit parameters
    attractorlist.append(attractors)#attractors for best fit parameters
    incornodeslist.append(incornodes)

    #experimental codes for trying out the analysis procedure in the package
    #inspector = util.analysis(pop, output_to_file=False)
    #inspector.f_distribution()
    #inspector.f_linearilty_convexity()
    #inspector.f_regression(degree=[1,1,2,2,3,3],interaction=[False,True,False,True,False,True])

algo_name = algo.get_name()
fstr_out = 'output_fittedparams_' + prob.netfilename[:-4] + '.txt'
with open(fstr_out, "wt") as f_out:#documenting results
    f_out.write("Algorithm: %s\nNumber of Generations: %d\nPopulation Size: %d\nNumber of Iterations: %d\n\n" % (algo_name, num_gens, popsize, n))
    wslist = []
    bslist = []
    for itno, (champparams, champfit, attractors, incornodes) in enumerate(zip(champparamslist, champfitlist, attractorlist, incornodeslist)):
        itno += 1
        f_out.write("%d. Best fitness: %f\n" % (itno, champfit))
        ws = np.array(champparams[:prob.num_edges - len(prob.forcew)])#list of weights for each links
        ws = ((ws + prob.maxindeg) / (prob.maxindeg + 1)).astype(np.int32)
        bs = np.array(champparams[prob.num_edges - len(prob.forcew):])#list of basal activities for each nodes
        bmask = [1 if s % 2 == 1 else -1 for s in bs]
        bs = ((bs + 1) / 2).astype(np.int32)
        bs = np.multiply(bs, bmask)
        f_out.write('Weights:\n')
        for w in ws:
            f_out.write('%d\t' % int(w))
        f_out.write('\nBasal Activities:\n')
        for b in bs:
            f_out.write('%d\t' % int(b))
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
    champfitfreq = 100 * champfitlist.count(bestchampfit) / n#in %
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
with open('pknsbestfitness.txt', 'a') as f:
    f.write(prob.netfilename)
    f.write(': %f\n' % min(champfitlist))

#if the simulation contained a perfect fitness, add to the list of successful pkns
if 0 in champfitlist:
    with open('workinigpkn.txt', 'a') as f:
        f.write(prob.netfilename + '\n')

rec.close()