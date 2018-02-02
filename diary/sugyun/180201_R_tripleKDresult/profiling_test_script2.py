score_dic = {}
with open('ab_senR_result .txt','r') as f:
    for i, line in enumerate(f):
        mod = (i+1)%5
#        print(mod)
        if mod == 1:
            continue
        elif mod == 2:
            line = line.strip().split(':')
            genelist = line[1].strip().split('\t')
        elif mod == 3:
            continue
        elif mod == 4:
            line = line.strip().split(':')
            score_dic[genelist] = line[1].strip()
with open('score_result.txt','w') as g:
    for i in score_dic:
        g.write('%s\t%s\n'%(i,score_dic[i]))