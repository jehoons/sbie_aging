from PyGMO import algorithm, population, island
import numpy as np

from optweight import OptWeight

num_gens = 100
prob = OptWeight(dim=9, i_dim=9)
algo = algorithm.sga(gen=num_gens)  # de, sga
pop = population(prob, 50)  # the last argument is population size.
isl = island(algo, pop) # island # != pop #
algo_name = algo.get_name()

isl.evolve(10) # evolve island 10 times
isl.join()
fchamp = isl.population.champion.f
# xchamp = ', '.join([str(elem) for elem in isl.population.champion.x])
print("[%s #%d]" % (algo_name, num_gens))
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
    f_out.write("Algorithm: %s\n" % algo_name)
    f_out.write("Best fitness: %f\n" % (isl.population.champion.f[0]))

    xchamp = isl.population.champion.x
    for i, param in enumerate(xchamp):
        f_out.write("%f,\n" % np.round(param))
        # end of for
# end of with
