import itertools
import numpy as np
def c(List):
    com = []
    for L in [1,2,3]:
        sam = []
        for x in itertools.combinations(List, L):
            sam.append(list(x))
        if len(sam) > 10000:
            choice_index = (list(np.random.choice(len(sam),10000)))
            com.extend([sam[x] for x in choice_index])
        else:
            com.extend(sam)
    return(com)