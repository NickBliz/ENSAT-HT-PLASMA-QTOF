sparrow$`sample name 3`=paste(sparrow$`sample name`, sparrow$group, sparrow$RUN)
sparrow$`sample name 3`[which(sparrow$group=="QC plasmamonster" | sparrow$group=="Qc Nick")]=paste(sparrow$`sample name`[which(sparrow$group=="QC plasmamonster" | sparrow$group=="Qc Nick")], sparrow$group[which(sparrow$group=="QC plasmamonster" | sparrow$group=="Qc Nick")])

#BRUNIUS
library(batchCorr)
sparrow1=sparrow[which(sparrow$RUN=="run1"),]
cm=checkmate.m[which(sparrow$RUN=="run1"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor1=QC1_DegrCorr$TestFeatsCorr
uncor1=cm
sparrow1=sparrow[which(sparrow$RUN=="run2"),]
cm=checkmate.m[which(sparrow$RUN=="run2"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor2=QC1_DegrCorr$TestFeatsCorr
uncor2=cm
sparrow1=sparrow[which(sparrow$RUN=="run3"),]
cm=checkmate.m[which(sparrow$RUN=="run3"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor3=QC1_DegrCorr$TestFeatsCorr
uncor3=cm
sparrow1=sparrow[which(sparrow$RUN=="run4"),]
cm=checkmate.m[which(sparrow$RUN=="run4"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor4=QC1_DegrCorr$TestFeatsCorr
uncor4=cm
sparrow1=sparrow[which(sparrow$RUN=="run5"),]
cm=checkmate.m[which(sparrow$RUN=="run5"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor5=QC1_DegrCorr$TestFeatsCorr
uncor5=cm
sparrow1=sparrow[which(sparrow$RUN=="run6"),]
cm=checkmate.m[which(sparrow$RUN=="run6"),]
cm=cm[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
sparrow1=sparrow1[order(as.numeric(as.character(sparrow1$`QTOF No`))),]
QC1_seq=as.numeric(as.character(sparrow1$`QTOF No`))
QC1_group=as.character(sparrow1$`sample name 3`)
set.seed(1)
QC1_DegrCorr <- correctDrift(peakTable = cm, injections = QC1_seq, sampleGroups = QC1_group, QCID = 'QC plasmamonster QC plasmamonster', RefID="Qc Nick Qc Nick", G = 1:10)
Brunius_cor6=QC1_DegrCorr$TestFeatsCorr
uncor6=cm
Bru=rbind(Brunius_cor1, Brunius_cor3, Brunius_cor4, Brunius_cor2, Brunius_cor5, Brunius_cor6)
sparrow_s=sparrow[order(as.numeric(as.character(sparrow$BATCH)),as.numeric(as.character(sparrow$`QTOF No`))),]
Bru=Bru[order(as.numeric(as.character(sparrow_s$TrueNo))),]
set.seed(1)
normData <- normalizeBatches(Bru, checkmate.m, batches = sparrow$BATCH, sampleGroup = as.character(sparrow$group), refGroup = "QC plasmamonster", population = 'all', CVlimit = 0, FCLimit = 0, medianZero = 'min')
Bru1 <- normData$peakTable
set.seed(1)
normData <- normalizeBatches(Bru, checkmate.m, batches = sparrow$BATCH, sampleGroup = as.character(sparrow$group), refGroup = "QC plasmamonster", population = 'all', medianZero = 'min')
Bru2 <- normData$peakTable


#RUVs
library(RUVSeq)
sparrow1=sparrow
QC1_group=as.character(sparrow1$`sample name 3`)
seq <- newSeqExpressionSet(as.matrix(t(checkmate.m)))
controls <- rownames(seq)
forty=rep(-1,38)
each_name=function(x){
  c(which(QC1_group==x), forty)
}
all_name=lapply(QC1_group[-which(QC1_group=="QC plasmamonster QC plasmamonster" | QC1_group=="Qc Nick Qc Nick")], each_name)
all_name=all_name[-which(duplicated(all_name))]
differences=do.call(rbind,all_name)
differences=rbind(which(QC1_group=="QC plasmamonster QC plasmamonster"), which(QC1_group=="Qc Nick Qc Nick"), differences)
differences[2,39]=-1
differences[2,40]=-1
set.seed(1)
seqRUVs <- RUVs(seq, controls, k=1, differences)
RUVs_cor=t(normCounts(seqRUVs))

#Wehrens
library(BatchCorrMetabolomics)
each_doBC=function(x){
  M2c <- doBC(checkmate.m[,x], ref.idx = which(sparrow$group=="QC plasmamonster"),
              batch.idx = sparrow$BATCH, method = "lm",
              seq.idx = as.numeric(as.character(sparrow$`QTOF No`)), imputeVal = 0, correctionFormula = formula("X ~ S * B"))
}
set.seed(1)
all_doBC=lapply(1:as.numeric(ncol(checkmate.m)), each_doBC)
all_BC=do.call(cbind, all_doBC)
all_BC=cbind.data.frame(sparrow$group, all_BC)

each_doBC=function(x){
  M2c <- doBC(checkmate.m[,x], ref.idx = which(sparrow$group=="QC plasmamonster"),
              batch.idx = sparrow$BATCH, method = "rlm",
              seq.idx = as.numeric(as.character(sparrow$`QTOF No`)), imputeVal = 0, correctionFormula = formula("X ~ S * B"))
}
set.seed(1)
all_doBC=lapply(1:as.numeric(ncol(checkmate.m)), each_doBC)
all_BC_rlm=do.call(cbind, all_doBC)
all_BC_rlm=cbind.data.frame(sparrow$group, all_BC_rlm)

#QC_RSC
#in R 3.4.1
library(pmp)
#classes <- MTBLS79$Class  #"QC", "C", OR "S" character
#batch <- MTBLS79$Batch    #"1" character
#order <- c(1:ncol(MTBLS79))  #QTOF NO AS NUMERIC integer

#out <- QCRSC(df = MTBLS79[1:10, ], order = order, batch = MTBLS79$Batch, classes = MTBLS79$Class, spar = 0, minQC = 4)
#mdat <- matrix(rep(0:10, 1720), nrow = 172, ncol = 10, byrow = TRUE)
#out <- QCRSC(df = mdat, order = order, batch = MTBLS79$Batch, classes = MTBLS79$Class, spar = 0, log=FALSE, minQC = 4, qc_label = "QC")

group <- as.character(sparrow$group)  #"QC", "C", OR "S" character
BATCH <- as.character(sparrow$BATCH)    #"1" character
QTOFNo <- as.numeric(as.character(sparrow$`QTOF No`))  #QTOF NO AS NUMERIC integer
#set.seed(42), identical with setseed 1 and 42
set.seed(1)
QCRSC_checkmate <- QCRSC(df = checkmate.m, order = QTOFNo, batch = BATCH, classes = group, spar = 0, log=FALSE, minQC = 4, qc_label = "QC plasmamonster")
tq=t(QCRSC_checkmate)
tq[is.na(tq)]=0

#AVERAGING: LOAD PREPROCESSED DATA. AV=cbind.data.frame(sparrow$group, checkmate.m), AV_sparrow=sparrow, AV0=checkmate
AV_i=AV[which(row.names(AV) %in% row.names(RUVs_cor)),]
AV_sparrow_i=AV_sparrow[which(row.names(AV) %in% row.names(RUVs_cor)),]
AV0_i=AV0[which(row.names(AV) %in% row.names(RUVs_cor)),]
RUVs_cor_i=RUVs_cor[which(row.names(RUVs_cor) %in% row.names(AV)),]
Brunius_cor_i1=Bru1[which(row.names(Bru1) %in% row.names(AV)),]
Brunius_cor_i2=Bru2[which(row.names(Bru2) %in% row.names(AV)),]
uncor_i=checkmate.m
uncor_i=uncor_i[which(row.names(uncor_i) %in% row.names(AV)),]
all_BC_i=all_BC[which(row.names(all_BC) %in% row.names(AV)),]
all_BCr_i=all_BC_rlm[which(row.names(all_BC_rlm) %in% row.names(AV)),]
sparrow_i=sparrow[which(row.names(sparrow) %in% row.names(AV)),]
checkmate_i=checkmate[which(row.names(sparrow) %in% row.names(AV)),]
all_BC_i$`sparrow$group`=all_BCr_i$`sparrow$group`=NULL
tq_i=tq[which(row.names(tq) %in% row.names(AV)),]
all_BC_i=as.matrix(all_BC_i)
all_BCr_i=as.matrix(all_BCr_i)

#BETWEEN
bn=c(rep(1,57), rep(2,62), rep(3,61), rep(4,35), rep(5,61), rep(6,58))
bn=as.factor(bn)

#checkmate.m=as.matrix(checkmate)
#checkmate.m[which(checkmate.m==0)]=min(checkmate.m[which(checkmate.m>0)])

#checkmate.m=Brunius_cor_i  #starting dataset. One of: AV, uncor_i, all_BC_i, all_BCr_i, Brunius_cor_i, RUVs_cor_i
#checkmate.m$`sparrow$group`=NULL

#normalize by batch average
c_cor1=sweep(AV_i[which(bn==1),], 2, colMeans(AV_i[which(bn==1),]) ,'/')
c_cor2=sweep(AV_i[which(bn==2),], 2, colMeans(AV_i[which(bn==2),]) ,'/')
c_cor3=sweep(AV_i[which(bn==3),], 2, colMeans(AV_i[which(bn==3),]) ,'/')
c_cor4=sweep(AV_i[which(bn==4),], 2, colMeans(AV_i[which(bn==4),]) ,'/')
c_cor5=sweep(AV_i[which(bn==5),], 2, colMeans(AV_i[which(bn==5),]) ,'/')
c_cor6=sweep(AV_i[which(bn==6),], 2, colMeans(AV_i[which(bn==6),]) ,'/')
#checkmate.m=rbind(c_cor1, c_cor2, c_cor3, c_cor4, c_cor5, c_cor6)

batcor=rbind(c_cor1, c_cor2, c_cor3, c_cor4, c_cor5, c_cor6)
#checkmate=checkmate[which(row.names(checkmate) %in% row.names(AV)),]   #NON-AV CHECKMATE
#CM=EITHER CHECKMATE OR AV0 to find original NAs.
#cm=checkmate

c_cor11=sweep(AV[which(bn==1),], 2, colMeans(AV[which(bn==1 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
c_cor21=sweep(AV[which(bn==2),], 2, colMeans(AV[which(bn==2 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
c_cor31=sweep(AV[which(bn==3),], 2, colMeans(AV[which(bn==3 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
c_cor41=sweep(AV[which(bn==4),], 2, colMeans(AV[which(bn==4 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
c_cor51=sweep(AV[which(bn==5),], 2, colMeans(AV[which(bn==5 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
c_cor61=sweep(AV[which(bn==6),], 2, colMeans(AV[which(bn==6 & sparrow_i$group=="QC plasmamonster"),]) ,'/')
#checkmate.m=rbind(c_cor1, c_cor2, c_cor3, c_cor4, c_cor5, c_cor6)

batcor1=rbind(c_cor11, c_cor21, c_cor31, c_cor41, c_cor51, c_cor61)

#Combat-seq
gr=sparrow_i$group
levels(gr)=c("1", "2", "3", "4", "5", "6", "7", "8")
batch=as.numeric(as.character(bn))
count_matrix=t(AV_i)
Combat_seq_AV_G <- ComBat_seq(t(AV_i), batch=as.numeric(as.character(bn)), group=as.numeric(as.character(gr)), covar_mod = NULL, full_mod=TRUE, shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL)
Combat_seq_AV <- ComBat_seq(t(AV_i), batch=as.numeric(as.character(bn)), group=NULL, full_mod=FALSE)
#Error in dispCoxReid(y, design = design, offset = offset, subset = subset,  : no data rows with required number of counts
#CONCLUSION: COMBAT-SEQ NOT SUITABLE AS VALUES ARE TOO SMALL AT THIS STAGE.


RUVs_cor_i[checkmate_i==0]=Brunius_cor_i1[checkmate_i==0]=Brunius_cor_i2[checkmate_i==0]=all_BC_i[checkmate_i==0]=all_BCr_i[checkmate_i==0]=tq_i[checkmate_i==0]=0
batcor[AV0_i==0]=0
batcor1[AV0_i==0]=0
AV_mc=batcor
AV_mc1=batcor1

# BE CAREFUL: FRPA2-0215 SHOULD BE EXCLUDED! GYDR-0145 ALSO. PLWW-0927 ARE BOTH PRE PPGL. CHOOSE EARLIEST.



library(mixOmics)
X=AV_mc1
#X$`sparrow$group`=NULL
pca.ENP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
Y=sparrow_i$BATCH
plotIndiv(pca.ENP, group = Y, title = "QTOF", ind.names = FALSE,legend = F,cex=7, comp = c(1,2))


pI=plotIndiv(pca.ENP, group = Y, title = "QTOF", ind.names = FALSE,legend = F,cex=7, comp = c(1,2), centroid = T)
#allD4=c(round(sqrt((pI$df$x0[which(pI$df$pch==1)][1]-pI$df$x0[which(pI$df$pch==2)][1])^2 + (pI$df$y0[which(pI$df$pch==1)][1]-pI$df$y0[which(pI$df$pch==2)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==1)][1]-pI$df$x0[which(pI$df$pch==3)][1])^2 + (pI$df$y0[which(pI$df$pch==1)][1]-pI$df$y0[which(pI$df$pch==3)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==1)][1]-pI$df$x0[which(pI$df$pch==4)][1])^2 + (pI$df$y0[which(pI$df$pch==1)][1]-pI$df$y0[which(pI$df$pch==4)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==1)][1]-pI$df$x0[which(pI$df$pch==5)][1])^2 + (pI$df$y0[which(pI$df$pch==1)][1]-pI$df$y0[which(pI$df$pch==5)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==1)][1]-pI$df$x0[which(pI$df$pch==6)][1])^2 + (pI$df$y0[which(pI$df$pch==1)][1]-pI$df$y0[which(pI$df$pch==6)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==2)][1]-pI$df$x0[which(pI$df$pch==3)][1])^2 + (pI$df$y0[which(pI$df$pch==2)][1]-pI$df$y0[which(pI$df$pch==3)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==2)][1]-pI$df$x0[which(pI$df$pch==4)][1])^2 + (pI$df$y0[which(pI$df$pch==2)][1]-pI$df$y0[which(pI$df$pch==4)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==2)][1]-pI$df$x0[which(pI$df$pch==5)][1])^2 + (pI$df$y0[which(pI$df$pch==2)][1]-pI$df$y0[which(pI$df$pch==5)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==2)][1]-pI$df$x0[which(pI$df$pch==6)][1])^2 + (pI$df$y0[which(pI$df$pch==2)][1]-pI$df$y0[which(pI$df$pch==6)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==3)][1]-pI$df$x0[which(pI$df$pch==4)][1])^2 + (pI$df$y0[which(pI$df$pch==3)][1]-pI$df$y0[which(pI$df$pch==4)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==3)][1]-pI$df$x0[which(pI$df$pch==5)][1])^2 + (pI$df$y0[which(pI$df$pch==3)][1]-pI$df$y0[which(pI$df$pch==5)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==3)][1]-pI$df$x0[which(pI$df$pch==6)][1])^2 + (pI$df$y0[which(pI$df$pch==3)][1]-pI$df$y0[which(pI$df$pch==6)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==4)][1]-pI$df$x0[which(pI$df$pch==5)][1])^2 + (pI$df$y0[which(pI$df$pch==4)][1]-pI$df$y0[which(pI$df$pch==5)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==4)][1]-pI$df$x0[which(pI$df$pch==6)][1])^2 + (pI$df$y0[which(pI$df$pch==4)][1]-pI$df$y0[which(pI$df$pch==6)][1])^2), 7)
#        ,round(sqrt((pI$df$x0[which(pI$df$pch==5)][1]-pI$df$x0[which(pI$df$pch==6)][1])^2 + (pI$df$y0[which(pI$df$pch==5)][1]-pI$df$y0[which(pI$df$pch==6)][1])^2), 7))

#dist_AV_combat_NP=allD4
#wilcox.test(dist_uncor, dist_AV_combat_NP, paired = T, alternative = "greater")
kruskal.test(pI$df$x~sparrow_i$BATCH)
kruskal.test(pI$df$y~sparrow_i$BATCH)

X=AV_mc1[which(sparrow_i$BATCH=="1"),]
pca.ENP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#Y=sparrow_i$group[which(sparrow_i$BATCH=="1")]
pI=plotIndiv(pca.ENP, title = "QTOF", ind.names = FALSE,legend = F,cex=7, comp = c(1,2), centroid = T)

cor.test(pI$df$x, as.numeric(as.character(sparrow_i$`QTOF No`[which(sparrow_i$BATCH=="1")])), method = "kendall")
cor.test(pI$df$y, as.numeric(as.character(sparrow_i$`QTOF No`[which(sparrow_i$BATCH=="1")])), method = "kendall")


#IS
plot(X[which(sparrow_i$group=="QC plasmamonster"),as.numeric(ncol(X))],col=Y[which(sparrow_i$group=="QC plasmamonster")], main= "QC plasmamonster 126.04999",xlim=c(0,25))
plot(X[which(sparrow_i$group=="QC plasmamonster"),as.numeric(ncol(X))-1],col=Y[which(sparrow_i$group=="QC plasmamonster")], main= "QC plasmamonster 169.10330",xlim=c(0,25))
plot(X[which(sparrow_i$group=="QC plasmamonster"),as.numeric(ncol(X))-2],col=Y[which(sparrow_i$group=="QC plasmamonster")], main= "QC plasmamonster 183.08242",xlim=c(0,25))



#RSD
hay=function(x) {
  sd(X[which(sparrow_i$group=="QC plasmamonster"),x], na.rm = T)
}
sdQC_plasmamonster=sapply((1:(as.numeric(ncol(X)))), hay)
bay=function(x) {
  mean(X[which(sparrow_i$group=="QC plasmamonster"),x], na.rm = T)
}
meanQC_plasmamonster=sapply((1:(as.numeric(ncol(X)))), bay)
cvQC_plasmamonster=sdQC_plasmamonster/meanQC_plasmamonster
#cvQC_plasmamonster[]
length(cvQC_plasmamonster)
length(which(cvQC_plasmamonster<0.3))
median(cvQC_plasmamonster, na.rm = T)
mad(cvQC_plasmamonster, na.rm = T)
median(cvQC_plasmamonster[cvQC_plasmamonster<0.3], na.rm = T)
mad(cvQC_plasmamonster[cvQC_plasmamonster<0.3], na.rm = T)
#RSD
hay=function(x) {
  sd(X[which(sparrow_i$group=="Qc Nick"),x], na.rm = T)
}
sdQC_plasmamonster=sapply((1:(as.numeric(ncol(X)))), hay)
bay=function(x) {
  mean(X[which(sparrow_i$group=="Qc Nick"),x], na.rm = T)
}
meanQC_plasmamonster=sapply((1:(as.numeric(ncol(X)))), bay)
cvQC_plasmamonster=sdQC_plasmamonster/meanQC_plasmamonster
length(cvQC_plasmamonster)
length(which(cvQC_plasmamonster<0.3))
median(cvQC_plasmamonster, na.rm = T)
mad(cvQC_plasmamonster, na.rm = T)
median(cvQC_plasmamonster[cvQC_plasmamonster<0.3], na.rm = T)
mad(cvQC_plasmamonster[cvQC_plasmamonster<0.3], na.rm = T)
