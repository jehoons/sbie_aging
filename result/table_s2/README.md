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

이렇게 해줍으로서 인슐린자극을 시뮬레이션 할 수 있습니다.

또한, 원본의 네트워크를 다음과 같이 T2D환자의 경우일때 `Insulin`과 `mTORC1_S6K1`에 모두 의존하도록 로직을 다음과 같이 수정합니다.

```
Node IRS_PIK3CA {
  logic = $T2D_PATIENT ? Insulin & !mTORC1_S6K1 : Insulin;
  rate_up = @logic ? $u_IRS_PIK3CA : 0;
  rate_down = @logic ? 0 : $d_IRS_PIK3CA;
}
```

수정된 네트워크는 `_mut.bnd` 파일에 저장합니다. 다음으로 우리는 정상과 당뇨환자의 경우를 각각 조사하고자 합니다. `$T2D_PATIENT` 변수는 `.cfg` 파일에서 변수 값을 0과 1로 설정해 주어야 하고 각각의 경우에 대해서 `_mut_normal.cfg`와 `_mut_t2d.cfg` 파일로 저장합니다.

**Step 2**. `vim Loic2016-model.cfg`, and set the value of the low_insulin or high_insulin variable to a value between 0 and 1.

**Step 3**. Run simulation.

```
MBSS_FormatTable.pl Loic2016-model_mut.bnd Loic2016-model_mut.cfg
```

**Step 4**. Postprocess the output

```
python postproc.py Loic2016-model_mut
```
