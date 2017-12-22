import numpy as np

sampleno = 2

data=[]
f = open('obj_attractor.txt', 'r')
for line in f:
    data.append(line.strip().split('\t'))

gendatalist = []

for cond in data:
    for i in range(sampleno):
        tmpipt = list(cond[0])
        tmpopt = cond[1]
        for idx, nd in enumerate(tmpipt):
            if nd == '-':
                tmpipt[idx] = str(np.random.choice(2))
            else:
                pass
        tmpipt = ''.join(str(s) for s in tmpipt)
        gendatalist.append([tmpipt, tmpopt])

with open('obj_attractor_sample.txt', "wt") as f_out:
    for attset in gendatalist:
        f_out.write("%s\t%s\n" % (attset[0], attset[1]))


