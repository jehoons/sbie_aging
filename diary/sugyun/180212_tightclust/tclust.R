# install.packages('tightClust')
library(tightClust)
library(Rserve)
Rserve(args='--vanilla')
data(tclust.test.data)
## find 10 tight clusters
ptm<-proc.time()
## k.min=25, tighter clusters will be found
## target=1 is used to save time, target=10 is recommended
tclust1<-tight.clust(MI_diff_complete, target=10, k.min=25, random.seed=12345)
proc.time()-ptm
## plot the heat map of cluster result
plot(tclust1)
## write the cluster result
write.tight.clust(tclust1)
write.table(tclust1$cluster,'tcluster1.txt')
ptm<-proc.time()
## k.min=10, looser clusters will be found
## target=1 is used to save time, target=5 is recommended
tclust2<-tight.clust(MI_diff_complete, target=5, k.min=10, random.seed=12345)
proc.time()-ptm
## plot the heat map of cluster result
plot(tclust2)
## write the cluster result
write.tight.clust(tclust2)
write.table(tclust2$cluster,'tcluster2.txt')