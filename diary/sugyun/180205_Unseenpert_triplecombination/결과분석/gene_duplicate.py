import operator

gene_duplicate = {}
with open('score_result.txt','r') as f:
    for i,line in enumerate(f):
        if i == 500:
            break
        else:
            line = line.split('\t')
            gene = line[0][:-1].split('/')
            for i in gene:
                if i in gene_duplicate:
                    gene_duplicate[i] += 1
                else:
                    gene_duplicate[i] = 1

print(gene_duplicate)
with open('gene_dup.txt', 'w') as g:
    sorted_x = sorted(gene_duplicate.items(), key=operator.itemgetter(1))
#    print(sorted_x)
    for key, value in sorted_x:
        g.write('%s\t%s\n'%(key, value))                    
                