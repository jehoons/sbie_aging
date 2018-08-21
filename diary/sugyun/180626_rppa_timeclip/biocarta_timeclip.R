# source("https://bioconductor.org/biocLite.R")
# biocLite("clipper")
# biocLite("BiocStyle")
# biocLite("car")
# biocLite("RSQLite")
# biocLite("RUnit")
# install.packages("stringi")
# install.packages("rrcov")
# install.packages("XML")
# install.packages("RCurl")
# install.packages("XMLRPC", repos="http://R-Forge.R-project.org")
# biocLite("RCytoscape")
# # biocLite("RCy3")
# install.packages("graph")
# 
# ## http://stat.columbia.edu/~gelman/bugsR/alternate_install.html  <- zip 파일의 package 설치하기
# install.packages("stringi")
# library(stringi)
# library(timeClip)
# library(RCytoscape)
# source("https://bioconductor.org/biocLite.R")
# biocLite("graphite")
# biocLite("org.Hs.eg.db")
# biocLite("RCy3")

###################################################
### read biocarta pathway 
###################################################
library(graphite)
library(RCy3)
pathwayDatabases()
biocarta <- pathways("hsapiens", "biocarta")
biocarta <- convertIdentifiers(biocarta, "symbol")
bi <- as.list(biocarta)

ke <- pathways("hsapiens", "kegg")
ke <- convertIdentifiers(ke, "symbol")
ke <- as.list(ke)

k <- bi
###################################################
### load data
###################################################
s <- read.csv(file="s_0.1_var_cutoff.csv", header=FALSE, sep=",",stringsAsFactors=FALSE)
s2 <- s[-1,]
colnames(s2) <- s[1,]
s =s2
s <- s[!duplicated(s[,1]), ]
s2 <- s[,-1]
rownames(s2) <- s[,1]
s=s2
mydata = data.matrix(s)


###################################################
### code chunk number 5: pathwayTimeCourseUsage
###################################################
library(timeClip)
times = as.numeric(colnames(s))
result = matrix(nrow=length((k)), ncol=length(times)+2)
for (i in c(1:length((k)))){
  g = pathwayGraph(k[[i]])
#  f <- nodes(graph)
#  f <- intersect(f, rownames(mydata))
#  if (length(f) == 0){
#    print('no intersect')
#    next}
  try({
    pathwayAnalysis <- pathwayTimeCourse(mydata, times, g,npc = 1,eqids=c(1,2,3,4))
    a = k[[i]]
    pathwayAnalysis$title = a@title
    result[i,1] <- a@title
    result[i,2] <- pathwayAnalysis$alpha
    result[i,3:6] <- pathwayAnalysis$pc
  },  silent = T)
}

write.csv(result, file = "s_0.1_biocarta.csv")


#timeClipped <- pathwayTimeCourse(hifExpression, times, hifGraph, npc=1, eqids=c(2,3,4,6,8,10,12,14,16,18,19,20,21,22))

# ###################################################
# ### code chunk number 6: timeClipUsage
# ###################################################
# #times <- as.numeric(colnames(hifExpression))
# #timeClipped <- timeClipSpaced(hifExpression, times, hifGraph, npc=1, eqids=c(2,3,4,6,8,10,12,14,16,18,19,20,21,22))
# g = pathwayGraph(k$`nf-kb signaling pathway`)
# timeClipped <- timeClipSpaced(mydata, times, g, npc=1, eqids=c(1,2,3,4))
# 
# timeClipped$clipped[,1:5]
# 
# 
# ###################################################
# ### code chunk number 7: sintetizeYourResult
# ###################################################
# clipped <- prunePaths(timeClipped$clipped, thr=0.2)
# clipped[,1:5]
# 
# plotTimeInCytoscape(mydata,timeClipped, 1)

