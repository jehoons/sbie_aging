library(CellNOptR)

pknmodel = readSIF("C:\\Users\\JSK\\Desktop\\17Dec_wk1\\PKN_v16_1.sif")
ptm <- proc.time()#mark the starting time of each simulation
cnodata = readMIDAS("C:\\Users\\JSK\\Desktop\\17Dec_wk1\\MIDAS_v16_test.csv")
cnolist = makeCNOlist(cnodata, subfield = FALSE)

# cnolistn = normaliseCNOlist(CNOlist = cnolist,
#                             mode = "time",
#                             verbose = FALSE)

model = preprocessing(cnolist, pknmodel)
initBstring <- rep(1, length(model$reacID))

opt1 <- gaBinaryT1(CNOlist = cnolist,
                   model = model,
                   initBstring = initBstring,
                   verbose = FALSE,
                   popSize = 50000, 
                   maxGens = 100000,
                   stallGenMax = 1000)

cutAndPlotResultsT1(model = model,
                    CNOlist = cnolist,
                    bString = opt1$bString,
                    plotPDF = FALSE)

print(proc.time() - ptm)#checking the time taken to run this simulation





