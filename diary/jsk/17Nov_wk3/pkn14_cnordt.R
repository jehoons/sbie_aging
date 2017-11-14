library(CellNOptR)
library(CNORdt)

#the cnordt FOR PKNv14
pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk3\\PKN_v14.sif")
cnodata = readMIDAS("C:\\Users\\JSK\\Desktop\\Git\\sbie_aging\\diary\\jsk\\17Nov_wk3\\MIDAS_v14_test.csv")
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
                   popSize = 100000)

cutAndPlotResultsDT(model = model,
                    CNOlist = cnolistn,
                    bString = opt1$bString,
                    plotPDF = FALSE,
                    boolUpdates = 5,
                    lowerB = .99,
                    upperB = 1.01)

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




