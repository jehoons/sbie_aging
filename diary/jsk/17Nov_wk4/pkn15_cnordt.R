library(CellNOptR)
library(CNORdt)

pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk4\\PKN_v15.sif")

optlist <- c()
scorelist <- c()
#datano = 0
#while (datano < 40960) {
sp = sort(sample(0:40959, 500, replace = F))
for (spno in sp) {
  ptm <- proc.time()#mark the starting time of each simulation
  cnodata = readMIDAS(paste("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk4\\MIDAS_v15_", toString(spno), ".csv", sep = ""))
  #datano = datano + 1
  cnolist = makeCNOlist(cnodata, subfield = FALSE)
  cnolistn = normaliseCNOlist(CNOlist = cnolist,
                              mode = "time",
                              verbose = FALSE)
  
  model = preprocessing(cnolistn, pknmodel)
  initBstring <- rep(1, length(model$reacID))
  
  opt1 <- gaBinaryDT(CNOlist = cnolistn,
                     model = model,
                     initBstring = initBstring,
                     verbose = FALSE,
                     boolUpdates = 5,
                     lowerB = .99,
                     upperB = 1.01,
                     popSize = 1000, 
                     maxGens = 10000,
                     stallGenMax = 500)
  
  # cutAndPlotResultsDT(model = model,
  #                     CNOlist = cnolistn,
  #                     bString = opt1$bString,
  #                     plotPDF = FALSE,
  #                     boolUpdates = 5,
  #                     lowerB = .99,
  #                     upperB = 1.01)
  
  score = computeScoreDT(cnolistn, 
                         model, 
                         bString = opt1$bString, 
                         boolUpdates = 5, 
                         lowerB = .99, 
                         upperB = 1.01)
  scorelist = c(scorelist, score)
  #optlist = c(optlist, opt1)# the opts can be obtained after data pruning using the score list
  print(proc.time() - ptm)#checking the time taken to run this simulation
}
# writeScaffold(
#   modelComprExpanded=model,
#   optimResT1=opt1,
#   optimResT2=NA,
#   modelOriginal=pknmodel,
#   CNOlist=cnolistn
# )
# 
# writeNetwork(
#   modelOriginal=pknmodel,
#   modelComprExpanded=model,
#   optimResT1=opt1,
#   CNOlist=cnolistn
# )




