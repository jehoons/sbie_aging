import json
dic = {}
dic['simulation_condition'] = {'sample_size':1000,
                               'repeat':500,
                               'lam':0.1,
                               'step_limit':100}
# dic['input_condition'] = [{'node':'x1', 'state':0}]
dic['input_condition'] = []
# dic['PKN'] = [{'source': 'x4', 'interaction': '((inhibit))', 'target': 'x1'},
#               {'source': 'x1', 'interaction': '((inhibit))', 'target': 'x2'},
#               {'source': 'x2', 'interaction': '((activate))', 'target': 'x3'},
#               {'source': 'x3', 'interaction': '((inhibit))', 'target': 'x4'},
#               {'source': 'x1', 'interaction': '((activate))', 'target': 'x3'},
#               {'source': 'x4', 'interaction': '((inhibit))', 'target': 'x3'}]

dic['PKN'] = [{'source': 'x4', 'interaction': '((inhibit))', 'target': 'x1'},
              {'source': 'x1', 'interaction': '((inhibit))', 'target': 'x2'},
              {'source': 'x2', 'interaction': '((activate))', 'target': 'x3'},
              {'source': 'x3', 'interaction': '((inhibit))', 'target': 'x4'},
              {'source': 'x1', 'interaction': '((activate))', 'target': 'x3'}]

dic["Node"] = [{'name':'x1', 'num':0},
               {'name':'x2', 'num':1},
               {'name':'x3', 'num':2},
               {'name':'x4', 'num':3}]
dic["train_graph"] = [{'before':'0010','pert': 0,'state':0,'after':'0110'},
                      {'before':'1010','pert': 0,'state':0,'after':'0110'},
                      {'before':'0010','pert': 1,'state':0,'after':'0010'},
                      {'before':'1100','pert': 1,'state':0,'after':'1000'},
                      {'before':'1011','pert': 1,'state':0,'after':'1011'},
                      {'before':'1011','pert': 1,'state':0,'after':'0010'},
                      {'before':'1011','pert': 1,'state':0,'after':'1000'},
                      {'before':'1010','pert': 1,'state':0,'after':'1010'},
                      {'before':'0010','pert': 2,'state':0,'after':'0101'},
                      {'before':'1100','pert': 2,'state':0,'after':'0101'},
                      {'before':'1011','pert': 2,'state':0,'after':'0101'},
                      {'before':'1010','pert': 2,'state':0,'after':'0101'},
                      {'before':'0010','pert': 3,'state':0,'after':'1010'},
                      {'before':'1010','pert': 3,'state':0,'after':'1010'}]
with open('data.txt', 'w') as outfile:
    # json.dump(data, outfile)
    json.dump(dic, outfile)