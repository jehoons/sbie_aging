library(CellNOptR)
library(CNORdt)

#the cnordt FOR PKNv14
pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk3\\PKN_v14.sif")
cnodata = readMIDAS("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk3\\MIDAS_v14.csv")
cnolist = makeCNOlist(cnodata, subfield = FALSE)

model = preprocessing(cnolist, pknmodel)
initBstring <- rep(1, length(model$reacID))

opt1 <- gaBinaryDT(CNOlist = cnolist,
                   model = model,
                   initBstring = initBstring,
                   verbose = FALSE,
                   boolUpdates = 5,
                   maxTime = 30,
                   lowerB = 0,
                   upperB = 3,
                   popSize = 10000)

cutAndPlotResultsDT(model = model,
                    CNOlist = cnolist,
                    bString = opt1$bString,
                    plotPDF = FALSE,
                    boolUpdates = 5,
                    lowerB = 0,
                    upperB = 3)

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




