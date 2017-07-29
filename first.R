#comparing t-test of data 
#manifest data contains information on each cpg
#pheno contains data on tissue type for BETRNET or MEMO >?
#mset.swan contains beta values for the samples given
#what do I need to do construct a table with cpgs as rows 
#were trying to determine 

#prolem 1 beternet and memo arent different

oncol <- getM(mset.swan[,"9979553161_R06C01"])
tocol <- getM(mset.swan[,"9985131138_R02C01"])
t.test(oncol, tocol)
 c(94,260,261,268,273)
 
 #look at cpgs that are present in both data sets
 cpgs = intersect(rownames(BETRNETmset),rownames(MEMOset))
 BETRNETmset = BETRNETmset[cpgs,]
 ids = c(1:367)[pheno$Tissue_Type == "SQ"]
 #ids contains the column indices with SQ tissue
 ids = setdiff(ids, c(94,260,261,268,273))
 #dum con
 dum = as.character(pheno$X[ids])
 SQMval = getM(BETRNETmset[,dum])
 all(colnames(SQMval) == dum)
 
 
 #ids contains pure SQuamus samples
 #get 100 cpgs markers 
 #get 
 
 cardiaids = c(1:367)[pheno$Tissue_Type == "cardia"]
 
 cardiaLAB = as.character(pheno$X[cardiaids])
 cardVal = getM(BETRNETmset[,cardiaLAB])
 dif = cardVal - SQMval
 meanrow = rowMeans(cardVal)
 meanrowSQM = rowMeans(SQMval)
 
 difmeans = meanrow - meanrowSQM
 
 d <- density(difmeans)
 plot(d)
 a = list(cardia = cardVal[2,], SQ = SQMval[2,])
 boxplot(a)
 
 
 
 sum(abs(meanrow - meanrowSQM) > 5)
 #dmp finder
 inputclass = c(rep("SQ",52),rep("CARD",9))
 total = cbind(SQMval, cardVal)
 find = dmpFinder(dat = total ,  pheno = inputclass,  type = "categorical",shrinkVar = TRUE)
 sum(abs(meanrow - meanrowSQM) > 5)
 sigRows = which(abs(meanrow - meanrowSQM) > 5)
 
 
 sigMeans = difmeans[sigRows]
 d <- density(sigMeans)
 plot(d)
 
 #mean m value for difference larger then five
 #squamus is fairly large
 largesig = which(meanrow - meanrowSQM > 5)
 val = meanrowSQM[largesig]
 val2 = meanrow[largesig]
 
 
 #BEids 
 BEids = c(1:367)[pheno$Tissue_Type == "BE"& pheno$Patient_Dx == "BE" & pheno$Project_ID == "BETRNet"]
 BEval = getM(BETRNETmset[,BEids])
 #rows of barret values at significant locations
 
 BEcpgs = BEval[sigRows, ]
 BEsd = apply(BEcpgs, 1, sd)
 CARDsd = apply(cardVal[sigRows, ],1,sd)
 SQsd = apply(SQMval[sigRows,],1,sd)
 
 
 #order makers by methylation value

 
 
 
 #MEMO interpretation
 
 
 
#volcano plot
 deltaM = difmeans[sigRows]
 
 qval = -log(find$qval[sigRows])
 
 plot(deltaM, qval, pch = 19, cex = 0.5)

 