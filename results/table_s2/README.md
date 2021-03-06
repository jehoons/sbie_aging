### Table 2. Test geroconversion model

This table is related to [#11](https://github.com/jehoons/sbie_aging/issues/11). [(Verlingue et al., 2016)][(Verlingue et al., 2016)] developed a boolean network model to describe type 2 diabetes. See also [supporting information][Verlingue et al., 2016-SI009].

#### (**A**) Test toy model

This model is a simple toy model provided by MaBoSS paper, and its structure is shown in (**S2FIG1**). This model consists of [toymodel.bnd](https://github.com/jehoons/sbie_aging/blob/master/result/table_s2/a/toymodel.bnd), [toymodel_slow.cfg](https://github.com/jehoons/sbie_aging/blob/master/result/table_s2/a/toymodel_slow.cfg), and [toymodel_fast.cfg](https://github.com/jehoons/sbie_aging/blob/master/result/table_s2/a/toymodel_fast.cfg). The `.bnd` file defines the simulation model, and `.cfg` defines the configuration how the simulation model runs. One or more `.cfg` files may exist for a ` .bnd` file.

<p align=center>
<img src="../../assets/img/maboss-toymodel.png" width="300" /><br>
S2FIG1. Toy model with three nodes
</p>

The code that describes the toy model is as follows:
```
Node C
{
    rate_up=0.0;
    rate_down=((not A) and (not B)) ? $escape : 0.0;
}

Node A
{
    rate_up=(C and (not B)) ? $Au : 0.0;
    rate_down=B ? $Ad : 0.0;
}

Node B
{
    rate_up=A ? $Au : 0.0;
    rate_down=A ? 0.0 : $Ad;
}
```

How to run simulation:

```
MaBoSS toymodel.bnd -c toymodel.cfg -o toymodel.out
# or
MBSS_FormatTable.pl toymodel.bnd toymodel.cfg
# to make a result figure run:
MBSS_TrajectoryFig.py toymodel
```

If you run `MBSS_FormatTable.pl` instead of `MaBoSS`, you can change the format of the output file into a more understandable table format.

The results for the two simulation configurations can be found in [toymodel_slow](a/toymodel_slow) and [toymodel_fast](a/toymodel_fast), respectively.

#### (**B**) Reproduce Figure 2 of the paper

Model network is shown as below.

<p align=center>
<img src="../../assets/img/verlingue2016-2.png" width="600" /><br>
S2FIG2. Model diagram from (Verlingue et al., 2016)
</p>

**S2FIG2** shows the response of the normal and diabetic networks to insulin stimulation. Follow the steps below to perform the simulation.

**Step 1. Mutate original network**

In order to create effects such as external stimulation or mutation, it is necessary to fix certain values from the original network. This is called mutation in MaBoSS. By doing this, you can, for example, add the effect of insulin into the model network.

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

**Step 2. Run simulation**

This step simulates a normal model and a diabetes model.

The normal model consists of the `Loic2016-model_mut_normal.bnd` file that defines the model structure and the `Loic2016-model_mut_normal.cfg` file that defines the simulation method. Run the following command to simulate the normal model:
```
MBSS_FormatTable.pl Loic2016-model_mut_normal.bnd Loic2016-model_mut_normal.cfg
```

The diabetes model consists of the `Loic2016-model_mut_t2d.bnd` file that defines the model structure and the `Loic2016-model_mut_t2d.cfg` file that defines the simulation method. Run the following command to simulate the diabetes model:
```
MBSS_FormatTable.pl Loic2016-model_mut_t2d.bnd Loic2016-model_mut_t2d.cfg
```

**Step 3. Postprocess the output**

By default, the output is given as the probability that the attractor will occur. The `postproc.py` module converts this to the probability that each node value will be activated. The conversion result is stored in the `_probtraj_table_processed.csv` file. In order to post-process the simulation result, the path where the result is stored should be input to the program as follows.
```
python postproc.py Loic2016-model_normal
python postproc.py Loic2016-model_t2d
```


[(Verlingue et al., 2016)]: https://github.com/jehoons/sbie_aging/files/792007/Verlingue.et.al.-.2016.-.A.comprehensive.approach.pdf
[Verlingue et al., 2016-SI009]: https://github.com/jehoons/sbie_aging/files/792011/Verlingue.et.al.-.2016.-.A.comprehensive.approach-supInfo009.docx
