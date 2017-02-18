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
**Step 1**. Mutate original network
In order to create effects such as external stimulation or mutation, it is necessary to fix certain values from the original network. This is called mutation in MaBoSS. By doing this, you can, for example, add the effect of adding insulin into the model network.

```
MBSS_MutBndCfg.pl Loic2016-model.bnd Loic2016-model.cfg 'Insulin'
```

This will create `_mut.bnd` and` _mut.cfg` files. The `MBSS_MutBndCfg.pl` can help you with simple tasks such as setting the value of a specific variable, but adding complex logic should be done by manual work. In the case of T2D patients, the `IRS_PIK3CA` depends on both the `Insulin` and the` mTORC1_S6K1`, so the logic equation is modified as follows.

```
Node IRS_PIK3CA {
  logic = $T2D_PATIENT ? Insulin & !mTORC1_S6K1 : Insulin;
  rate_up = @logic ? $u_IRS_PIK3CA : 0;
  rate_down = @logic ? 0 : $d_IRS_PIK3CA;
}
```

For convenience, names `_mut.bnd` and` _mut.cfg` are `_mut_normal.bnd` and` _mut_normal.cfg` respectively. Then copy it and save it as `_mut_t2d.bnd`,` _mut_t2d.cfg`. These are the configuration files for normal and diabetic patients, respectively. Set `$ T2D_PATIENT` variable to 0 in` _mut_normal.cfg` file and 1 in `_mut_t2d.cfg` file with` vim`, respectively.

**Step 3**. Run simulation.

```
MBSS_FormatTable.pl Loic2016-model_mut_normal.bnd Loic2016-model_mut_normal.cfg
MBSS_FormatTable.pl Loic2016-model_mut_t2d.bnd Loic2016-model_mut_t2d.cfg
```

**Step 4**. Postprocess the output.
```
python postproc.py Loic2016-model_normal
python postproc.py Loic2016-model_t2d
```
