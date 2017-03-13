### Table 1. Network for aging study 

#### (**A**) Regulation network
A regulation network is a network that contains information about activation and inhibition, but it does not contain information about logic functions. `A, +1, B` means that A activates B, and` B, -1, C` means that B inhibits C.

#### (**B**) Logical network generation 
In the regulation network (A), the regulatory relation between two nodes is described, but when the target node is controlled by two or more source nodes, its logical function has not yet been defined. Therefore, I created a logical network through simple assumptions. The rules are as follows: When two or more source nodes activate, they are assumed to be joined with OR operation. When two or more source nodes inhibit, they are assumed to be joined with OR operation. Activation group and inhibition group are joined with AND operation .

#### (**C**) Hello attractor analysis 
With python module `boolean3` I simulated the generated logical network (B) and then identified attractors and its basin sizes from the network.