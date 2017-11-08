library(CellNOptR)
library(CNORdt)

#the cnordt FOR PKNv14
pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\17Nov_wk2\\string_PKN.sif")
cnodata = readMIDAS("C:\\Users\\JSK\\Desktop\\17Nov_wk2\\string_PKN.csv")
cnolist = makeCNOlist(cnodata, subfield = FALSE)

model = preprocessing(cnolist, pknmodel)
initBstring <- rep(1, length(model$reacID))

opt1 <- gaBinaryDT(CNOlist = cnolist,
                   model = model,
                   initBstring = initBstring,
                   verbose = FALSE,
                   boolUpdates = 4,
                   maxTime = 30,
                   lowerB = .8,
                   upperB = 10,
                   popSize = 1000000)

cutAndPlotResultsDT(model = model,
                    CNOlist = cnolist,
                    bString = opt1$bString,
                    plotPDF = FALSE,
                    boolUpdates = 4,
                    lowerB = .8,
                    upperB = 10)

writeScaffold(
  modelComprExpanded=model,
  optimResT1=opt1,
  optimResT2=NA,
  modelOriginal=pknmodel,
  CNOlist=cnolist
)

writeNetwork(
  modelOriginal=pknmodel,
  modelComprExpanded=model,
  optimResT1=opt1,
  CNOlist=cnolist
)




