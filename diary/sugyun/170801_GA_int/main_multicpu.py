from PyGMO import *
import numpy as np
from optweight import OptWeight
# archi 8 island (ind pop) 20 gen 100 evolve 10
prob = OptWeight(dim = 9, i_dim = 9)
algo = algorithm.sga(gen = 100)
archi = archipelago(algo,prob,8,100) #constructs an archipelago having 8 islands with populations of 100 individuals.
# print min([isl.population.champion.f for isl in archi])
archi.evolve(5) #starts the asynchronous generalized island model.
                 #each of the 8 islands will call algo 10 times and try to migrate in between calls
archi.join()
print([isl.population.champion.f for isl in archi])
fchamp = isl.population.champion.f
# xchamp = ', '.join([str(elem) for elem in isl.population.champion.x])
print("- Champion fitness: %f" % (fchamp[0]))
print(np.array(isl.population.champion.x))
# import time
#
# now = time.localtime()
#
# fstr_out = "output_%s_%d%02d%02d_%02dh%02dm%02ds.txt" % ('singleobj',
#                                                          now.tm_year,
#                                                          now.tm_mon,
#                                                          now.tm_mday,
#                                                          now.tm_hour,
#                                                          now.tm_min,
#                                                          now.tm_sec)
#
# with open(fstr_out, "wt") as f_out:
#     f_out.write("100	%f	" % (isl.population.champion.f[0]))
#
#     xchamp = isl.population.champion.x
#     for i, param in enumerate(xchamp):
#         f_out.write("%f\t" % np.round(param))
        # end of for
# end of with