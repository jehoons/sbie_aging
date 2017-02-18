### Table 2. Test aging model (Loic et al., 2016)
This result is related to [#11](https://github.com/jehoons/sbie_aging/issues/11). 
#### (**A**) Test toy model

This model is a simple toy model provided by MaBoSS. The `.bnd` file defines the simulation model, and `.cfg` defines the configuration how the simulation model runs. How to run simulation:
```
MaBoSS toymodel.bnd -c toymodel.cfg -o toymodel.out
# or 
MBSS_FormatTable.pl toymodel.bnd toymodel.cfg
MBSS_TrajectoryFig.py toymodel
```

If you run `MBSS_FormatTable.pl`, you can change the format of the output file into a more understandable table format.

#### (**B**) Test Loic2016 model 

**Step 1**. Mutate original network to set constant node such as input or mutation.
```
MBSS_MutBndCfg.pl Loic2016-model.bnd Loic2016-model.cfg 'Insulin'
```
**Step 2**. `vim Loic2016-model.cfg`, and set the value of the low_insulin or high_insulin variable to a value between 0 and 1.

**Step 3**. Run simulation.
```
MBSS_FormatTable.pl Loic2016-model_mut.bnd Loic2016-model_mut.cfg
```

**Step 4**. Postprocess the output 
```
python postproc.py Loic2016-model_mut
```


