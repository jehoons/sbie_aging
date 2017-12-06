library(CellNOptR)

pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Dec_wk1\\PKN_v16_1.sif")

optlist <- c()
scorelist <- c()
#datano = 0
#while (datano < 204800) {
sp = sort(sample(0:204799, 500, replace = F))
for (spno in sp) {
  ptm <- proc.time()#mark the starting time of each simulation
  cnodata = readMIDAS(paste("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Dec_wk1\\MIDAS_v16_", toString(spno), ".csv", sep = ""))
  #datano = datano + 1
  cnolist = makeCNOlist(cnodata, subfield = FALSE)
  cnolistn = normaliseCNOlist(CNOlist = cnolist,
                              mode = "time",
                              verbose = FALSE)
  
  model = preprocessing(cnolistn, pknmodel)
  initBstring <- rep(1, length(model$reacID))
  
  opt1 <- gaBinaryT1(CNOlist = cnolistn,
                     model = model,
                     initBstring = initBstring,
                     verbose = FALSE,
                     popSize = 1000, 
                     maxGens = 10000,
                     stallGenMax = 500)
  
  # cutAndPlotResultsT1(model = model,
  #                     CNOlist = cnolistn,
  #                     bString = opt1$bString,
  #                     plotPDF = FALSE)
  
  score = computeScoreT1(CNOlist = cnolistn, 
                         model = model, 
                         bString = opt1$bString)
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




