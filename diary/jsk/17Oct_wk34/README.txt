1. network_pruning은 첫번째 시도; PKN을 배제한 시도여서 패스해도됨.(tempnetwork가 이 시도의 result)
2. positive_pruning은 PKN(Aging_network)에 존재하는 link들 중에서 MI_raw의 데이터를 이용하여 주어진 cutoff 이상의 mutual information을 가지는 link들만 선별하는 network pruning임. positive_pruning_v2 코드를 더 깔끔하게 정리한 버전; 하는 일은 같음.
3. negative_pruning은 PKN(Aging_network)에 존재하는 link들 중에서 MI_raw의 데이터를 이용하여 주어진 cutoff 이하의 mutual information을 가지는 link들을 지워낸 후 남는 PKN을 반환하는 network pruning임. negative_pruning_v2 코드를 더 깔끔하게 정리한 버전; 하는 일은 같음.
positive_pruning과 negative_pruning 둘다 cutofflist에 들어있는 cutoff 숫자들만 조절해주면 원하는 cutoff에서의 결과를 도출할 수 있음.