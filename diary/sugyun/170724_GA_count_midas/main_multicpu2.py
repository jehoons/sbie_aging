from PyGMO import *
import numpy as np
from optweight import OptWeight
# archi 8 island (ind pop) 20 gen 100 evolve 10
prob = OptWeight(dim = 9)
algo = algorithm.sga(gen = 500)
archi = archipelago(algo,prob,20,50)
print min([isl.population.champion.f for isl in archi])
archi.evolve(10)
isl.join()
print min([isl.population.champion.f for isl in archi])
fchamp = isl.population.champion.f
# xchamp = ', '.join([str(elem) for elem in isl.population.champion.x])
print("- Champion fitness: %f" % (fchamp[0]))
print(np.array(isl.population.champion.x))
import time

now = time.localtime()

fstr_out = "output_%s_%d%02d%02d_%02dh%02dm%02ds.txt" % ('singleobj',
                                                         now.tm_year,
                                                         now.tm_mon,
                                                         now.tm_mday,
                                                         now.tm_hour,
                                                         now.tm_min,
                                                         now.tm_sec)

with open(fstr_out, "wt") as f_out:
    f_out.write("0.1	%f	" % (isl.population.champion.f[0]))

    xchamp = isl.population.champion.x
    for i, param in enumerate(xchamp):
        f_out.write("%f\t" % np.round(param))
        # end of for
# end of with