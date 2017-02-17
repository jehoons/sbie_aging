### Table 2. Test aging model (Loic et al., 2016)
#### (**A**) Test toy model 
이 모델은 MaBoSS에서 제공하는 간단한 토이모델입니다. `*.bnd` 파일은 시뮬레이션 모델을 정의하고, `*.cfg`는 시뮬레이션 모델을 실행할 조건을 정의합니다. 구성파일은 아래와 같습니다. 

시뮬레이션 실행방법: 
```
MaBoSS toymodel.bnd -c toymodel.cfg -o toymodel.out
# or 
MBSS_FormatTable.pl toymodel.bnd toymodel.cfg
MBSS_TrajectoryFig.py toymodel
```

`MBSS_FormatTable.pl`을 실행하는 경우에는 출력파일의 형식을 보다 이해하기 쉬운 테이블형식으로 바꾸어 잔다는 장점이 있다. 

#### (**B**) Test Loic2016 model 
시뮬레이션 실행방법: 
```
MBSS_FormatTable.pl Loic2016-model.bnd Loic2016-model.cfg
MBSS_TrajectoryFig.py Loic2016-model
```

