import itertools
def c(List):
    com = []
    for L in range(len(List)+1):
        for x in itertools.combinations(List, L):
            com.append(list(x))
    return(com)