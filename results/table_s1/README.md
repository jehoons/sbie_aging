### Table 1. Network for aging study 

#### (**A**) Network 

A regulation network is a network that contains information about activation and inhibition, but it does not contain information about logic functions. `A, +1, B` means that A activates B, and` B, -1, C` means that B inhibits C.

*A1. Regulatory network, initial version*

*A2. Regulatory network, processed*

*A3. Logical model, generated by simple rule*

In the regulation network (A1), the regulatory relation between two nodes is described, but when the target node is controlled by two or more source nodes, its logical function has not yet been defined. Therefore, I created a logical network through simple assumptions. The rules are as follows: When two or more source nodes activate, they are assumed to be joined with OR operation. When two or more source nodes inhibit, they are assumed to be joined with OR operation. Activation group and inhibition group are joined with AND operation .

*A4. Logical model*

* updated version 1 - initial version 

* updated version 2 - Delete PP2A node, DNA damage as input 

* updated version 3 - Delete PTEN node, DNA damage as input 

* updated version 4 - Delete S6K -IRS link, DNA damage as input

#### (**B**) Analysis the network

With python module `boolean3` I simulated the generated logical network (A3,A4) and then identified attractors and its basin sizes from the network. 

*B1. Input combination generation*

Nodes `IGF1` and `NowNutEx` are inputs. Then, I generated input combinations as (False,False), (False,True), (True,False), (True,True). 

*B2. Results for the input combinations shown in B1*

Simulation is performed using the input combinations shwon in B1. The result is stored with JSON style.

*B3. Summarized results of B2*

*B4. Similarity matrix of state vectors*

*B5. Summarized results with detailed information about attractors of B2*

#### (**C**) Analysis the network - Phenotype 

Original comment from 안수균:

*C1. Phenotype analysis from the attractors(B5)*
> Open the "sim_result_analysis-v2.py" and run. This python file imports a module which is "phenotype_detector.py". In addition, "phenotype_detector.py" imports a module which is "attractor_phenotype.py". You have to put the file, gene symbol list, and the logic to get the phenotype of the attractor. The result will be provided as txt file. by A.S.K.

The code `test_c.py` uses previous simulation result(B5) to identify phenotypes of the calculated attractors.

input:

`c/b5-simul-log_v2.txt` 

output:

`c/b5-updated1-boolsim-results-summary_phenotype_v2.txt`


