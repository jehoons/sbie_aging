from PyGMO import *
import numpy as np
from optweight import OptWeight
# archi 8 island (ind pop) 20 gen 100 evolve 10
prob = OptWeight(dim = 9,i_dim = 9)
algo = algorithm.sea(gen = 100)
archi = archipelago(algo,prob,8,20)
print min([isl.population.champion.f for isl in archi])
archi.evolve(10)
isl.join()
print min([isl.population.champion.f for isl in archi])
fchamp = isl.population.champion.f
# xchamp = ', '.join([str(elem) for elem in isl.population.champion.x])
print("- Champion fitness: %f" % (fchamp[0]))
print(np.array(isl.population.champion.x))