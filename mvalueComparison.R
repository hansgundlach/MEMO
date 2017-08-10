

#retrieving msets and pheno data
BETRNETmset = load("C:/Users/hansg17/Documents/methData2017/mset.swan.BETRNet.rda")
BETRNETmset = mset.swan
MEMOset =  load("C:/Users/hansg17/Documents/methData2017/msetMEMO.rda")
MEMOset = mset.swan
BLOODset = load("C:/Users/hansg17/Documents/methData2017/bloodM.rda",verbose = T)
BLOODset = bloodM
MEMOpheno <- read.csv("~/methData2017/new_pdat_072516.csv") 

pheno <- save("C:/Users/hansg17/Documents/methData2017/pheno.BETRNet." , verbose = T)
phenoset = load("C:/Users/hansg17/Documents/methData2017/pheno.BETRNet.rda" , verbose = T)
pheno = pheno

msetHans22 <- read.csv("~/methData2017/msetHans22.csv")
recMset = load("C:/Users/hansg17/Documents/methData2017/recM.rda",verbose = T)
recMset = recM
#contains objects DNAsoln and age.BE30
matched30BE = load("C:/Users/hansg17/Documents/methData2017/30matchedBE.rda",verbose = T)




 #look at cpgs that are present in both BETRNET and MEMOset and augment manifest as well
 intercpgs = intersect(rownames(BETRNETmset),rownames(MEMOset))
 BETRNETmset = getM(BETRNETmset[intercpgs,])
 MEMOset = getM(MEMOset[intercpgs,])
 cpic = intersect(rownames(manifestData),rownames(BETRNETmset))
 manifestData = manifestData[cpic,]
 
 
 #ids contains indices of SQamnus samples
 ids = c(1:367)[pheno$Tissue_Type == "SQ"]
 #setdiff removes sqamus examples contaminated with barrets
 ids = setdiff(ids, c(94,260,261,268,273))
 #dum contains the patient sample names for given indices
 sampName = as.character(pheno$X[ids])
 SQMval = getM(BETRNETmset[,sampName])
 
 
 cardiaids = c(1:367)[pheno$Tissue_Type == "cardia"]
 #sample names for cardia samples
 cardiaLAB = as.character(pheno$X[cardiaids])
 cardVal = getM(BETRNETmset[,cardiaLAB])
 
 #retrieve fundus values from BETRNET
 fundusids = c(1:367)[pheno$Tissue_Type == "fundus"]
 fundusVal = getM(BETRNETmset[,fundusids])
 fundusSig = fundusVal[sigRows,]
 
 
 #look at cpgs that are common to both blood and normal tissue
 #blood methylation data
 universal = intersect(cpic,rownames(BLOODset))
 bloodVal = BLOODset[universal,]
 
 
 
 
 #========================================================================================================
 
 #dmp finder
 inputclass = c(rep("SQ",52),rep("CARD",9))
 total = cbind(SQMval, cardVal)
 # finds contains the cpgs that are most differentialy methylate dbetween SQamus a Cardia
 find = dmpFinder(dat = total ,  pheno = inputclass,  type = "categorical",shrinkVar = TRUE)
 sum(abs(meanrow - meanrowSQM) > 5)
 #look at the rows where mean cardia and squamus mean differ by 5 in M-values 
 sigRows =  rownames(find[1:200,])
   
 #BEids 
 BEids = c(1:367)[pheno$Tissue_Type == "BE"& pheno$Patient_Dx == "BE" & pheno$Project_ID == "BETRNet"]
 BEval = getM(BETRNETmset[,BEids])
 #rows of barret values at significant locations
 
 
 #SQsd conatins the row variation for a paricular cpg
 BEsig = BEval[sigRows, ]
 BEsd = apply(BEsig, 1, sd)
 CARDsd = apply(cardVal[sigRows, ],1,sd)
 SQsd = apply(SQMval[sigRows,],1,sd)
 
 BEmean = rowMeans(BEsig)
 
 #boxplot
 box = list(BE = BEsig[1, ], CARD = cardVal[sigRows[1],], SQ =  SQMval[sigRows[1],], FUND = fundusSig[1,] )
 boxplot(box)
 t.test(BEsig[2, ], cardVal[sigRows[2],])
 

#testing for drift cpgs in BE
corBE = numeric(265)
corSQ = numeric(265)
corCARD = numeric(265)

#for loops test if cpgs correlate with age ie are clock cpgs 
for(i in 1:nrow(sigBEval)) {
  corBE[i] = cor.test(sigBEval[i, ], pheno$Age_at_biopsy[BEids])$p.value
  
}
for(i in 1:nrow(sigSQMval)){
corSQ[i] = cor.test(sigSQMval[i, ], pheno$Age_at_biopsy[ids])$p.value
}
#correlation between age for SQ
for(i in 1:nrow(sigSQMval)){
  corSQ[i] = cor.test(sigSQMval[i, ], pheno$Age_at_biopsy[ids])$p.value
}

#correlation between age for cardia
for(i in 1:nrow(sigCARD)){
  corCARD[i] = cor.test(sigCARD[i, ], pheno$Age_at_biopsy[cardiaids])$p.value
}


#filter out q-value below 0.5

#don't we want to look at most significantly different indicators instead of ones that differ by 5


#do side by side comparison for matched samples
sigdifBeta = which( abs(meanrow(cardValB)- meanrow(SQMvalB))>.5)
 

test2 = cor.test(rowMeans(fundusSig), rowMeans(BEval[sigRows,]))

#different tissue cluster withhemselves
#body exon1 tss 5primeUTR

#look at cpgs associated with TP63
#checking for specific genes in markers
som = manifestData[sigRows, "USC_RefGene_Name"]
duck = strsplit(som, ";")
unlistduck = unlist(duck)
geneL = unique(unlistduck)
sort(geneL)


#PCA Analysis for BETRNET methylation levels at sigRow cpgs
#main = rbind(t(SQMval[sigRows,]),t(BEsig), t(cardVal[sigRows,]), t(funduVal[sigRows,]),t(BEnine),t(cardsamp))

main = rbind(t(SQMval[finalcpgs,]),t(BEval[finalcpgs,]), t(cardVal[finalcpgs,]), t(fundusVal[finalcpgs,]))
maincor = cor(main)
eigenmain = eigen(maincor)
proj = main %*% eigenmain$vector
color = c(rep(2,52),rep(5,64), rep(4,9),rep(6,12))
plot(proj[,1],proj[,2],col = color,pch = 19,xlab = "PCA1",ylab = "PCA2")
#this is going to seperate 
legend("topleft", legend = c("SQ","BE","CARD","Fundus"), col = c(2,5,4,6,7),pch = 19 ,bty = "n")
title("Tissue Methylation PCA400 Analysis")
#mds analysis (very similar)
mdsPlot(t(main), sampGroups =  color,pch = 19)



#look at relations in the significant cpgs identified
#do all the cpgs go up in unison or is their variety

islands = manifestData[sigRows, ]$Islands_Name
#quite a large amount of OpenSea in islandRelat 
islandRelat = manifestData[names(sigRows),]$Relation_to_Island
#relation to island for cpgs in sigRows 
#seems like an abnormal number are open sea 
opengenenames = manifestData[names(sigRows), "UCSC_RefGene_Name"][manifestData[names(sigRows),]$Relation_to_Island == "OpenSea"]
stuff = manifestData[opencpgnames,"Relation_to_Island" ]   
#gene names cleaned and sorted
cleaned = sort(unique(unlist(strsplit(opengenenames,";"))))


#look at correlation up esophagus in MEMO data
datsix = MEMOpheno[MEMOpheno$Patient_ID == 619,]

MEMOsmall = MEMOset[sigRows,]
#set of samples
#memcol = MEMOMval[,sixsamp]
#cardsamp contains cardia data 
cardsamp = as.vector(MEMOsmall[,datsix$X[13]])

sixsamp = as.character(datsix$X[-13])

BEnine  = MEMOsmall[,sixsamp[1:13]]
#correlation with all samples
meanBETRNET = rowMeans(sigCARD)

distlist = c(35,34,33,32,31,30,29,28,27,26,25,24,20)

#correlation with carida for BE samples in patient 619
#look at varience around mean 
meanCARD = rowMeans(cardVal[sigRows,])
meanBE = rowMeans(BEval[sigRows,])
corarraynot = c()
#adjust for z-score
BEsds = apply(BEval[sigRows,],1,sd)
CARDsds = apply(cardVal[sigRows,],1,sd)
#patient 619 has 13 samples 
#if I change cardia samples correlation goes away but correlation pattern does not 

corarray = c()
  #7 random 6 higly correlated 5 random 4 uhh 3 good 
for(i in 1:13)
{
  print(i)
  #cardVal[sigRows]
  #index = sample(1:265,265,replace = FALSE)
  corarray[i] = cor((cardVal[sigRows, 2] - meanCARD)/CARDsds, 
                    as.vector((MEMOsmall[,sixsamp[i] ]- meanCARD)/CARDsds))
  #corarray[i] = cor(avgMEMCARD, as.vector(MEMOsmall[,sixsamp[i] ]))
  #corarraynot[i] = cor(cardsamp, as.vector(MEMOsmall[,sixsamp[i] ]))
}

plot(distlist[1:13], corarray,pch = 19)
title("Mean Correlation Patient 619 Cardia ")








#488 analysis
#######################################################################33333
pat488 = MEMOpheno[MEMOpheno$Patient_ID == 488,]
pat489 = MEMOpheno[MEMOpheno$Patient_ID == 489,]

samples = as.character(pat488$X)

dist488 = c(33,31, 0,37,35,0, 0 )
corarray488 = numeric(length(dist488))
for(i in 1:7)
{
  #index = sample(1:265,265,replace = FALSE)
  #corarray[i] = cor((cardsamp -meanCARD)/CARDsds, as.vector((MEMOsmall[,sixsamp[i] ]- meanBE)/BEsds))
  corarray488[i] = cor(memeCARD, as.vector(MEMOsmall[,samples[i] ]  - meanBE))
 # corarraynot488[i] = cor(cardsamp, as.vector(MEMOsmall[,sixsamp[i] ]))
}

plot(dist488, corarray488,pch = 19)



#squamus carida 489 488 look at BETRNET 8 and 22
#look at EAC in patient 8 
pat8 = pheno[pheno$Patient_ID == 8,]
pat8Mval = getM(BETRNETmset[sigRows,as.character(pat8$X)])
cor.test(pat8Mval[,1]-meanCARD,pat8Mval[,2])
cor.test(pat8Mval[,2],pat8Mval[,3])
cor.test(pat8Mval[,1],pat8Mval[,4])
cor.test(pat8Mval[,1]-meanCARD,pat8Mval[,4]-meanCARD)

#8cardia EAC    HGD    fundus
#22cardia HGD    EAC    BE     SQ     fundus BE     BE    
pat22 = pheno[pheno$Patient_ID == 22,]
pat22Mval = getM(BETRNETmset[sigRows,as.character(pat22$X)])
cor.test(pat8Mval[,1],pat8Mval[,3])
cor.test()

colnames(pat8Mval)  <- c("cardia","EAC","HGD","fundus")
res <- cor(pat8Mval)
round(res, 2)

colnames(pat22Mval)  <- c("cardia","HGD","EAC","BE","SQ","fundus","BE","BE")
res22 <- cor(pat22Mval)
round(res22,2)
plot(pat22Mval[,1],pat22Mval[,6],abline)


install.packages(corrplot)
library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)




#rectal analysis to compare 
recSig = recMset[sigRows,]
chosen = sigRows[125]
box = list(BE = BEval[chosen,], CARD = cardVal[chosen,], SQ =  SQMval[chosen,], FUND = fundusVal[chosen,], REC = recMset[chosen,],BL = bloodVal[chosen,])
boxplot(box)
#BE and rectal have decent correlation 
cor.test(rowMeans(BEval[sigRows,]),rowMeans(recMset[sigRows,]))
#Blood does not have high correlation 
cor.test(rowMeans(BEval[sigRows,]),rowMeans(bloodVal[sigRows,]))



#code below deals with finding diffential methylation markers between fundus and cardia
cardMEMOval = getM(MEMOset[,c("200325530003_R01C01","200394970056_R04C01","200325530006_R02C02")])
#cardStrom contains the fundus values 
cardStom = getM(MEMOset[,c("9979553161_R04C01", "9985131138_R06C01","9992571073_R01C01")])
#cardFun = cbind(cardVal,cardMEMOval,fundusVal,cardStom) when memo was combined
cardFun = cbind(cardVal,fundusVal)
inputFC = c(rep("CARD",9),rep("FUND",12))
#fundus cardia markers using combined MEMO 
fundfind = dmpFinder(dat = cardFun ,  pheno = inputFC,  type = "categorical",shrinkVar = TRUE)
#rownames(fundfind)[fundfind$qval < 5*1e-4]
meanFUND = rowMeans(fundusVal)
meanCARD = rowMeans(cardVal)
meanBE = rowMeans(BEval)

difFC = meanFUND[sigFC] - meanCARD[sigFC]



sigFC = rownames(fundfind[1:200,])

sigFC[grepl("LGR5", manifestData[sigFC, "UCSC_RefGene_Name"])]
cor.test(meanCARD[sigFC], meanBE[sigFC])
manifestData[sigFC,"UCSC_RefGene_Name"]
cor.test(meanCARD[realsigFC],meanFUND[realsigFC])
plot(meanCARD[sigFC],meanFUND[sigFC])
abline(0,1)
plot(c(rep(1,9),rep(2,12)), c(cardVal["cg15909056",], fundusVal["cg15909056",]),pch =19)

#experiment 
#BE and cardia actually differe significantly 
inputBE = c(rep("CARD",9),rep("BE",64))
beCARD = cbind(cardVal[,],BEval[,1:64])
cardBE = dmpFinder(dat = beCARD,  pheno = inputBE,  type = "categorical",shrinkVar = TRUE)
cardBEsig = rownames(cardBE[1:200,])

#none of these BE cARDIA differentialy methylated cpgs seem like drifters 
corBECARD = numeric(200)
#look if those are clock cpgs 
for(i in 1:length(cardBEsig)) {
  corBECARD[i] = cor(sigBEval[i, ], pheno$Age_at_biopsy[BEids])
  
}



#hierarchical clustering of patient 8 samples 
colnames(pat8Mval)
clusters <-hclust(dist(t(pat8Mval)))
plot(clusters)

#final list of cpgs combination of fundusSQ cpgs and caridaSQ cpgs 400 in length
finalcpgs = c(sigRows, sigFC)

sigpat22 = pat22Mval[finalcpgs,]
colnames(sigpat22) <- c("cardia","HGD","EAC","BE","SQ","fundus","BE","BE")
clust3 = hclust(dist(t(sigpat22)))   
plot(clust3)

mod = getM(MEMOset)
dug  <- c(as.character(MEMOset$Tissue_location),row.names =  NULL)
colnames(mod) <- dug
clust4 = hclust(dist(t(mod[finalcpgs,])))
plot(clust4)

#look at dwell time vs correlation to see if it
#induces some form of phylogenetic radiation 

#load mset with dwell times for each sample
#looking at all samples from MEMO that have a dwell time in msetHans
#msetHans22 <- read.csv("~/msetHans22.csv")
patDWELL = msetHans22$tissue_age[ msetHans22$BE == "yes" & !is.na(msetHans22$tissue_age) & msetHans22$HGD.LGD == "yes"]
namesDWELL = msetHans22$X[ msetHans22$BE == "yes" & !is.na(msetHans22$tissue_age)& msetHans22$HGD.LGD == "yes"]
#be careful using which and indices 
#serious = which(!is.na(patDWELL))
#serious2 = which(msetHans22$BE == "yes")
#serious3 = intersect(serious,serious2)
#DWELL = patDWELL
#namesDWELL = as.character(as.vector(patDWELL$X))
allDWELL = MEMOsmall[,namesDWELL]
allDWELL = MEMOsmall[,namesDWELL]
#msetHans22[serious,]$tissue_age
#corDWELL will contain the correlation with mean squamus for all BE MEMO samples 
corDWELL = numeric(length(namesDWELL))
meanSQM = meanrowSQM[sigRows]

#memeCARD is the mean cardia sample from MEMO
for(i in 1:length(namesDWELL)) {
  
    corDWELL[i] = cor(c(), c(MEMOsmall[,namesDWELL[i]]))
}


plot(patDWELL, corDWELL, pch=19)
plot()
title("Corr with Carida vs dwell time MEMO")

#plot points above and below correlation level
ding = corDWELL[corDWELL > 0.6]
stuff = DWELL[corDWELL > 0.6]
plot(stuff,ding,pch=19)


points(stuff,ding, col = 2, pch = 19)
cor.test(stuff,ding)
ding = which(corDWELL < .3)


#looking at commonality between markers 
#fundus cardia 
a = rownames(fundfind[1:200,])
#cardia BE diff makers
b = rownames(cardBE[1:200,])
#carida SQ
c = rownames(find[1:200,])
#seems to be no intersects in cpg sets
ab = intersect(a,b)
ac = intersect(a,c)
bc = intersect(b,c)

#list of MEMO cardia samples
MEMOcard= c("200325530006_R02C02", "200325530006_R06C01", "200394970056_R04C01", "200325530003_R01C01", "200325530006_R05C01")
avgMEMCARD = rowMeans(MEMOsmall[,MEMOcard])
plot(avgMEMCARD, meanCARD , pch =19)


BETRNET200 = getM(BETRNETmset[sigRows, ])
#analysis for BETRNET matched
matBETnames = pheno[pheno$DNA.soln %in% DNAsoln,]
#allBETmach has all the BETRNET samples with dwell times 
allBETmach  = BETRNET200[,names(matBETnames)]


#index = which(pheno$DNA.soln %in% DNAsoln)
#indexSQM = index + 1

#3allSQMmatch = BETRNET200[,indexSQM]

corBETdwell = numeric(ncol(allBETmach))

#compare squamus with cardia and barrett

for(i in 1:length(corBETdwell)) {
  
  corBETdwell[i] = cor(c(meanSQM), c(allBETmach[,matBETnames[i]]))
  
}
#why are squamus markers so well correlated 
plot(age.BE30 - age,corBETdwell,pch = 19)
cor.test(age.BE30 - age,corBETdwell)


age = c(8.679, 32.44, 12.77, 1.967, 4.664, 9.01 , 47.54, 26.76,  46.82 ,  49.29 , 28.66 , 17.96 ,19.82 ,59.04 ,31.75 , 34.67 ,42.03 ,49.24 ,26.58,
50.07, 27.53, 38.95, 25.03, 50.41, 40.71, 56.01, 29.19, 58.68, 40.77, 46.78)


#analysis based on seans packages
install.packages("pvclust")
library(pvclust)
try = pvclust(sigpat22)
#ape part of analysis 
install.packages("ape")
library(ape)

apeTRY = fastme.bal(dist(sigpat22))


#try doing the  data but for different markers all of it 
#look at fine structure correlation across samples
















#extenous old code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ask georg about taking out a few cpgs for regularization between samples


#color by age and see if they cluster

#blood meth value correlation blood takes up a lot of damages stuff

#centroid distance vs GEJ distance correlation 

#dwell time vs correlation with cardia

#test matched cardia squamus

#redo correlation graph with only BE

#familial BE



dif = cardVal - SQMval
meanrowCARD = rowMeans(cardVal)
meanrowSQM = rowMeans(SQMval)

difmeans = meanrowCARD - meanrowSQM
d <- density(difmeans)
plot(d)
a = list(cardia = cardVal[2,], SQ = SQMval[2,])
boxplot(a)

# old sigrows which(abs(meanrow - meanrowSQM) > 5)

sigMeans = difmeans[sigRows]
d <- density(sigMeans)
plot(d)

sum(find$qval[sigRows] > 0.01)
sort(find$qval[sigRows])
#volcano plot
deltaM = difmeans[sigRows]

qval = -log(find$qval[sigRows])

plot(deltaM, qval, pch = 19, cex = 0.5)

#mean m value for difference larger then five
#squamus is fairly large
largesig = which(meanrow - meanrowSQM > 5)
val = meanrowSQM[largesig]
val2 = meanrow[largesig]

#plot
sigCARD = cardVal[sigRows,]
sigSQMval = SQMval[sigRows,]
sigBEval = BEval[sigRows,]

plot(SQMval[sigRows,1], pch = 19, cex = 0.6, , ylim = c(-6,6) )
cor.test(meanrowSQM[sigRows], BEmean)

