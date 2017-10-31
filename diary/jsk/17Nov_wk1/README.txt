1. MI_calc은 raw data의 input condition을 여러가지 조합으로 (2,3,4개) 나눈 후에 mutual information을 계산한다. 시간이 매우 오래걸린다.
2. positive_pruning_v3 and negative_pruning_v3은 기존의 pruning 알고리즘을 모듈화 시켰고 여러가지 cutoff를 시도해보는 대신, 총 샘플 개수의 mutual information median에 가까운 mutual information 값을 cutoff으로 잡은 후 pruning을 하였다.
3. pruning_analysis은 MI_calc에서 나온 데이터셋을 positive pruning과 negative pruning을 돌려주는 mother function이다. 이것도 조금 오래걸린다.
4. MI_results에는 MI_calc에 대한 결과들이 들어 있고 network_pruning_result에는 pruning_analysis에 대한 결과들이 들어있다.
5. agingnetworkp는 network_pruning_results에 있는 모든 네트워크를 반영하는 최종 네트워크이다. 이를 만드는 코드가 net_network_gen이다