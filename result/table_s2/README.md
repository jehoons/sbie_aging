### Table 2. Test aging model (Loic et al., 2016)
#### (**A**) Test toy model 
이 모델은 MaBoSS에서 제공하는 간단한 토이모델입니다. `*.bnd` 파일은 시뮬레이션 모델을 정의하고, `*.cfg`는 시뮬레이션 모델을 실행할 조건을 정의합니다. 구성파일은 아래와 같습니다. 

```
Four_cycle.bnd
Four_cycle_FEscape.cfg
Four_cycle_SEscape.cfg
Four_cycle_SEscape_long.cfg
```

MaBoSS 시뮬레이션을 실행하는 방법은: 

```
MaBoSS Four_cycle.bnd -c Four_cycle_FEscape.cfg -o Four_cycle.out
```

#### (**B**) Test Loic2016 model

```
Loic2016-model.bnd
Loic2016-model.cfg
```

MaBoSS 시뮬레이션을 실행하는 방법은: 

```
MaBoSS Loic2016-model.bnd -c Loic2016-model.cfg -o Loic2016-model.out
```

