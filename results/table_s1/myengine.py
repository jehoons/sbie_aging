import sys
import json 
from ipdb import set_trace
from boolean3_addon import attr_cy

def run_engine(infile, outfile): 
    with open(infile, 'r') as fin: 
        indata = json.load(fin)
    res_list = []
    for i,ind in enumerate(indata): 
        if i == 0:             
            attr_cy.build("\n".join(ind['model']))
            import pyximport;
            pyximport.install(reload_support=True)
        res = attr_cy.run(
            samples=ind['samples'], steps=ind['steps'], debug=ind['debug'],
            on_states=ind['on_states'], off_states=ind['off_states'], 
            progress=True)

        # set_trace()
        res['input_condition'] = ind['input_condition']
        res_list.append(res)

    json.dump(res_list, open(outfile, 'w'), indent=4)

if __name__ == '__main__': 
    infile = sys.argv[1]
    outfile = sys.argv[2]
    run_engine(infile, outfile)
