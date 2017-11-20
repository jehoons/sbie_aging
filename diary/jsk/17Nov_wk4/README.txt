PKN_v15 has a slight modification including exclusion of INSR and RB-ABL link and the addition of DNAdamagex which represents DNA damage strong enough to cause apoptosis
data2midas.py will take the raw data and the list of manually selected antibody list in order to make total of 40960 MIDAS data.
pkn15_cnordt.r takes the MIDAS generated and conducts a CNORdt procedure on all of them; it saves the fitness score and output data for further analysis

There are several problems with this method:
-the computation time required to run this is estimated to be a month
-the fitting done on these data are not significant enough for comparison

In order to work around the computational time limitation, I selected a sample of 500 to run this simulation
The r code can be run after specifying the location suited to your device.