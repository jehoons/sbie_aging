import operator

score_dic = {}
with open('ab_senR_result .txt','r') as f:
    for i, line in enumerate(f):
        if line.startswith('total'):
            break
        mod = (i+1)%5
#        print(mod)
        if mod == 1:
            continue
        elif mod == 2:
            line = line.strip().split(':')
            genelist = line[1].strip()
            for i in genelist:
                genecom = ''
                new_genelist = genelist.split('\t')
                for k in new_genelist:
                    genecom = genecom+k+'/'
#            print(genecom)
        elif mod == 3:
            continue
        elif mod == 4:
            line = line.strip().split(':')
            score_dic[genecom] = float(line[1].strip())
with open('score_result.txt', 'w') as g:
    sorted_x = sorted(score_dic.items(), key=operator.itemgetter(1))
#    print(sorted_x)
    for key, value in sorted_x:
        g.write('%s\t%s\n'%(key, value))