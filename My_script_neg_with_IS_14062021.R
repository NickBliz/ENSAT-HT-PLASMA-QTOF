setwd("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/ALL_NEG/BLANKS")
EQP=read.table("EQP.tsv", sep="\t", header = T)
EQP_names <- read_csv("EQP_names.csv")
EQP_names=t(EQP_names)
EQP_names=EQP_names[-(1:2),]
colnames(EQP_names)=colnames(EQP)
EQP_b1d1=EQP[,which(EQP_names[10,]=="run1" & (EQP_names[9,]==1))]
EQP_b1d2=EQP[,which(EQP_names[10,]=="run1" & (EQP_names[9,]==2))]
EQP_b2d1=EQP[,which(EQP_names[10,]=="run3" & (EQP_names[9,]==1))]
EQP_b2d2=EQP[,which(EQP_names[10,]=="run3" & (EQP_names[9,]==2))]
EQP_b3d1=EQP[,which(EQP_names[10,]=="run4" & (EQP_names[9,]==1))]
EQP_b3d2=EQP[,which(EQP_names[10,]=="run4" & (EQP_names[9,]==2))]
EQP_b4d1=EQP[,which(EQP_names[10,]=="run2" & (EQP_names[9,]==1))]
EQP_b4d2=EQP[,which(EQP_names[10,]=="run2" & (EQP_names[9,]==2))]
EQP_b5d1=EQP[,which(EQP_names[10,]=="run5" & (EQP_names[9,]==1))]
EQP_b5d2=EQP[,which(EQP_names[10,]=="run5" & (EQP_names[9,]==2))]
EQP_b6d1=EQP[,which(EQP_names[10,]=="run6" & (EQP_names[9,]==1))]
EQP_b6d2=EQP[,which(EQP_names[10,]=="run6" & (EQP_names[9,]==2))]

EQP_b1d1=cbind(EQP[,(1:9)], EQP_b1d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b1d2=cbind(EQP[,(1:9)], EQP_b1d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b2d1=cbind(EQP[,(1:9)], EQP_b2d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b2d2=cbind(EQP[,(1:9)], EQP_b2d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b3d1=cbind(EQP[,(1:9)], EQP_b3d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b3d2=cbind(EQP[,(1:9)], EQP_b3d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b4d1=cbind(EQP[,(1:9)], EQP_b4d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b4d1=EQP_b4d1[,-29]
EQP_b4d2=cbind(EQP[,(1:9)], EQP_b4d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
#REMOVE ALL SAMPLES NOT RERUN, INCLUDING QCs AT BEGINNING!
EQP_b4d2=EQP_b4d2[,-(32:35)]
EQP_b5d1=cbind(EQP[,(1:9)], EQP_b5d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b5d2=cbind(EQP[,(1:9)], EQP_b5d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b6d1=cbind(EQP[,(1:9)], EQP_b6d1, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])
EQP_b6d2=cbind(EQP[,(1:9)], EQP_b6d2, EQP[,c(ncol(EQP), ncol(EQP)-1, ncol(EQP)-2)])

bl=read.table("EQP_blanks.tsv", sep="\t", header = T)
bl_b1d1=cbind(bl[,(1:9)], bl[,c(10,13)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b1d2=cbind(bl[,(1:9)], bl[,c(11,12)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b2d1=cbind(bl[,(1:9)], bl[,c(14,17)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b2d2=cbind(bl[,(1:9)], bl[,c(15,16)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b3d1=cbind(bl[,(1:9)], bl[,c(18,21)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b3d2=cbind(bl[,(1:9)], bl[,c(19,20)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d1=cbind(bl[,(1:9)], bl[,c(22,23)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d1=bl_b4d1[,-11]
bl_b4d2=cbind(bl[,(1:9)], bl[,c(22,23)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b4d2=bl_b4d2[,-12]
bl_b5d1=cbind(bl[,(1:9)], bl[,c(24,27)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b5d2=cbind(bl[,(1:9)], bl[,c(25,26)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b6d1=cbind(bl[,(1:9)], bl[,c(28,31)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])
bl_b6d2=cbind(bl[,(1:9)], bl[,c(29,30)], bl[,c(ncol(bl), ncol(bl)-1, ncol(bl)-2)])


write.table(EQP_b1d1, file='EQP_b1d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b1d2, file='EQP_b1d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b2d1, file='EQP_b2d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b2d2, file='EQP_b2d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b3d1, file='EQP_b3d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b3d2, file='EQP_b3d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b4d1, file='EQP_b4d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b4d2, file='EQP_b4d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b5d1, file='EQP_b5d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b5d2, file='EQP_b5d2.tsv', quote=FALSE, sep='\t')
write.table(EQP_b6d1, file='EQP_b6d1.tsv', quote=FALSE, sep='\t')
write.table(EQP_b6d2, file='EQP_b6d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b1d1, file='bl_b1d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b1d2, file='bl_b1d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b2d1, file='bl_b2d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b2d2, file='bl_b2d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b3d1, file='bl_b3d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b3d2, file='bl_b3d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b4d1, file='bl_b4d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b4d2, file='bl_b4d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b5d1, file='bl_b5d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b5d2, file='bl_b5d2.tsv', quote=FALSE, sep='\t')
write.table(bl_b6d1, file='bl_b6d1.tsv', quote=FALSE, sep='\t')
write.table(bl_b6d2, file='bl_b6d2.tsv', quote=FALSE, sep='\t')

add_blank_matches('EQP_b1d1.tsv', 'bl_b1d1.tsv')
add_blank_matches('EQP_b1d2.tsv', 'bl_b1d2.tsv')
add_blank_matches('EQP_b2d1.tsv', 'bl_b2d1.tsv')
add_blank_matches('EQP_b2d2.tsv', 'bl_b2d2.tsv')
add_blank_matches('EQP_b3d1.tsv', 'bl_b3d1.tsv')
add_blank_matches('EQP_b3d2.tsv', 'bl_b3d2.tsv')
add_blank_matches('EQP_b4d1.tsv', 'bl_b4d1.tsv')
add_blank_matches('EQP_b4d2.tsv', 'bl_b4d2.tsv')
add_blank_matches('EQP_b5d1.tsv', 'bl_b5d1.tsv')
add_blank_matches('EQP_b5d2.tsv', 'bl_b5d2.tsv')
add_blank_matches('EQP_b6d1.tsv', 'bl_b6d1.tsv')
add_blank_matches('EQP_b6d2.tsv', 'bl_b6d2.tsv')

EQP_b1d1_pbmatched=read.table("EQP_b1d1_pbmatched.tsv", sep="\t", header = T)
EQP_b1d2_pbmatched=read.table("EQP_b1d2_pbmatched.tsv", sep="\t", header = T)
EQP_b2d1_pbmatched=read.table("EQP_b2d1_pbmatched.tsv", sep="\t", header = T)
EQP_b2d2_pbmatched=read.table("EQP_b2d2_pbmatched.tsv", sep="\t", header = T)
EQP_b3d1_pbmatched=read.table("EQP_b3d1_pbmatched.tsv", sep="\t", header = T)
EQP_b3d2_pbmatched=read.table("EQP_b3d2_pbmatched.tsv", sep="\t", header = T)
EQP_b4d1_pbmatched=read.table("EQP_b4d1_pbmatched.tsv", sep="\t", header = T)
EQP_b4d2_pbmatched=read.table("EQP_b4d2_pbmatched.tsv", sep="\t", header = T)
EQP_b5d1_pbmatched=read.table("EQP_b5d1_pbmatched.tsv", sep="\t", header = T)
EQP_b5d2_pbmatched=read.table("EQP_b5d2_pbmatched.tsv", sep="\t", header = T)
EQP_b6d1_pbmatched=read.table("EQP_b6d1_pbmatched.tsv", sep="\t", header = T)
EQP_b6d2_pbmatched=read.table("EQP_b6d2_pbmatched.tsv", sep="\t", header = T)

EQP_matches=cbind(EQP_b1d1_pbmatched, EQP_b1d2_pbmatched, EQP_b2d1_pbmatched, EQP_b2d2_pbmatched, EQP_b3d1_pbmatched, EQP_b3d2_pbmatched, EQP_b4d1_pbmatched, EQP_b4d2_pbmatched, EQP_b5d1_pbmatched, EQP_b5d2_pbmatched, EQP_b6d1_pbmatched, EQP_b6d2_pbmatched)
colnames(EQP_matches)=make.unique(colnames(EQP_matches))
EQP_matches=cbind(EQP_matches$sample_min0, EQP_matches$blank_matches, EQP_matches$bf, EQP_matches$bf2, 
                  EQP_matches$sample_min0.1, EQP_matches$blank_matches.1, EQP_matches$bf.1, EQP_matches$bf2.1,
                  EQP_matches$sample_min0.2, EQP_matches$blank_matches.2, EQP_matches$bf.2, EQP_matches$bf2.2,
                  EQP_matches$sample_min0.3, EQP_matches$blank_matches.3, EQP_matches$bf.3, EQP_matches$bf2.3,
                  EQP_matches$sample_min0.3, EQP_matches$blank_matches.3, EQP_matches$bf.3, EQP_matches$bf2.3,
                  EQP_matches$sample_min0.4, EQP_matches$blank_matches.4, EQP_matches$bf.4, EQP_matches$bf2.4,
                  EQP_matches$sample_min0.5, EQP_matches$blank_matches.5, EQP_matches$bf.5, EQP_matches$bf2.5,
                  EQP_matches$sample_min0.6, EQP_matches$blank_matches.6, EQP_matches$bf.6, EQP_matches$bf2.6,
                  EQP_matches$sample_min0.7, EQP_matches$blank_matches.7, EQP_matches$bf.7, EQP_matches$bf2.7,
                  EQP_matches$sample_min0.8, EQP_matches$blank_matches.8, EQP_matches$bf.8, EQP_matches$bf2.8,
                  EQP_matches$sample_min0.9, EQP_matches$blank_matches.9, EQP_matches$bf.9, EQP_matches$bf2.9,
                  EQP_matches$sample_min0.10, EQP_matches$blank_matches.10, EQP_matches$bf.10, EQP_matches$bf2.10,
                  EQP_matches$sample_min0.11, EQP_matches$blank_matches.11, EQP_matches$bf.11, EQP_matches$bf2.11)
jolly=colnames(EQP_matches)=c("sample_min0", "blank_matches", "bf", "bf2", 
                              "sample_min0.1", "blank_matches.1", "bf.1", "bf2.1",
                              "sample_min0.2", "blank_matches.2", "bf.2", "bf2.2",
                              "sample_min0.3", "blank_matches.3", "bf.3", "bf2.3",
                              "sample_min0.3", "blank_matches.3", "bf.3", "bf2.3",
                              "sample_min0.4", "blank_matches.4", "bf.4", "bf2.4",
                              "sample_min0.5", "blank_matches.5", "bf.5", "bf2.5",
                              "sample_min0.6", "blank_matches.6", "bf.6", "bf2.6",
                              "sample_min0.7", "blank_matches.7", "bf.7", "bf2.7",
                              "sample_min0.8", "blank_matches.8", "bf.8", "bf2.8",
                              "sample_min0.9", "blank_matches.9", "bf.9", "bf2.9",
                              "sample_min0.10", "blank_matches.10", "bf.10", "bf2.10",
                              "sample_min0.11", "blank_matches.11", "bf.11", "bf2.11")
EQP_matches=data.frame(EQP_matches)
EQP_new <- filter_relevant(read.table('EQP.tsv',sep = "\t",header=TRUE))
jelly=EQP_new[-which((EQP_matches$sample_min0==0 & EQP_matches$blank_matches>=1) | (EQP_matches$blank_matches==1 & EQP_matches$bf<0) | (EQP_matches$bf<0 & EQP_matches$bf2<0) |
                       (EQP_matches$sample_min0.1==0 & EQP_matches$blank_matches.1>=1) | (EQP_matches$blank_matches.1==1 & EQP_matches$bf.1<0) | (EQP_matches$bf.1<0 & EQP_matches$bf2.1<0) |
                       (EQP_matches$sample_min0.2==0 & EQP_matches$blank_matches.2>=1) | (EQP_matches$blank_matches.2==1 & EQP_matches$bf.2<0) | (EQP_matches$bf.2<0 & EQP_matches$bf2.2<0) |
                       (EQP_matches$sample_min0.3==0 & EQP_matches$blank_matches.3>=1) | (EQP_matches$blank_matches.3==1 & EQP_matches$bf.3<0) | (EQP_matches$bf.3<0 & EQP_matches$bf2.3<0) |
                       (EQP_matches$sample_min0.4==0 & EQP_matches$blank_matches.4>=1) | (EQP_matches$blank_matches.4==1 & EQP_matches$bf.4<0) | (EQP_matches$bf.4<0 & EQP_matches$bf2.4<0) |
                       (EQP_matches$sample_min0.5==0 & EQP_matches$blank_matches.5>=1) | (EQP_matches$blank_matches.5==1 & EQP_matches$bf.5<0) | (EQP_matches$bf.5<0 & EQP_matches$bf2.5<0) |
                       (EQP_matches$sample_min0.6==0 & EQP_matches$blank_matches.6>=1) | (EQP_matches$blank_matches.6==1 & EQP_matches$bf.6<0) | (EQP_matches$bf.6<0 & EQP_matches$bf2.6<0) |
                       (EQP_matches$sample_min0.7==0 & EQP_matches$blank_matches.7>=1) | (EQP_matches$blank_matches.7==1 & EQP_matches$bf.7<0) | (EQP_matches$bf.7<0 & EQP_matches$bf2.7<0) |
                       (EQP_matches$sample_min0.8==0 & EQP_matches$blank_matches.8>=1) | (EQP_matches$blank_matches.8==1 & EQP_matches$bf.8<0) | (EQP_matches$bf.8<0 & EQP_matches$bf2.8<0) |
                       (EQP_matches$sample_min0.9==0 & EQP_matches$blank_matches.9>=1) | (EQP_matches$blank_matches.9==1 & EQP_matches$bf.9<0) | (EQP_matches$bf.9<0 & EQP_matches$bf2.9<0) |
                       (EQP_matches$sample_min0.10==0 & EQP_matches$blank_matches.10>=1) | (EQP_matches$blank_matches.10==1 & EQP_matches$bf.10<0) | (EQP_matches$bf.10<0 & EQP_matches$bf2.10<0) |
                       (EQP_matches$sample_min0.11==0 & EQP_matches$blank_matches.11>=1) | (EQP_matches$blank_matches.11==1 & EQP_matches$bf.11<0) | (EQP_matches$bf.11<0 & EQP_matches$bf2.11<0)),]

jelly=rbind(jelly,EQP_new[545,], EQP_new[1531,], EQP_new[1881,])
jel=t(jelly)
colnames(jel)=jel[1,]
tzel=jel[-c(1:9, 707:709),]

library(readr)
EQP_cl=read.delim(file="EQP_pb_cleaned_withMetaData.txt",sep = "\t",header=TRUE, dec = ".")
EQP_cle=t(EQP_cl)
colnames(EQP_cle)=EQP_cle[1,]
EQP_cle=EQP_cle[-1,]
EQP_cle=as.data.frame(EQP_cle)
EQP_cle$TrueNo=as.numeric(as.character(EQP_cle$TrueNo))
EQPc=EQP_cle[order(EQP_cle$TrueNo),]

tzelis=tzel[which(row.names(tzel) %in% row.names(EQPc)),]

num=function(x){
  which(row.names(tzelis)==row.names(EQPc)[x])
}
num_all=lapply(1:as.numeric(nrow(EQPc)), num)
numall=unlist(num_all)
tzelis=tzelis[numall,]

fab=function(x){
  batchX=EQPc[which(EQPc$BATCH==x),]
  #b1=batchX[,21:ncol(batchX)]
  b1=tzelis[row.names(tzelis) %in% row.names(batchX),]
  
  for (i in 1:ncol(b1)) {
    mi=b1[,i]
    mi=as.numeric(as.character(mi))
    for (j in 1:length(which(batchX$`F/B`=="forward"))) {
      if (mi[j]==0 | mi[length(mi)-j+1]==0){
        mi[j]=mi[length(mi)-j+1]=0
      } else {
        mi[j]=mean(c(mi[j], mi[length(mi)-j+1]))
      }}
    b1[,i]=mi
  }
  
  return(b1)
}

EQP_av=lapply(1:6, fab)
EQP_av_m=rbind(EQP_av[[1]],EQP_av[[2]],EQP_av[[3]],EQP_av[[4]],EQP_av[[5]],EQP_av[[6]])

make.num=function(x) {
  tzelis[,x]=as.numeric(as.character(tzelis[,x]))
}
tzelis=sapply(1:ncol(tzelis), make.num)
#replace wrongly corrected end of batch QCs with their true values
EQP_av_m[c(56,57,173,174,294,295,406,407,518,519,635,636),]=tzelis[c(56,57,173,174,294,295,406,407,518,519,635,636),1:ncol(tzelis)]
EQP_av_m=cbind(EQPc[,(1:20)], EQP_av_m)
EQP_av_m_m=EQP_av_m[-which(EQP_av_m$`F/B`=="backwards"),]

test=jel[which(row.names(EQP_av_m_m) %in% row.names(jel)),]
EQP_bp=t(EQP_av_m_m) #REMEMBER! <jel> has your variable names!

jack=t(EQP_bp)
jack=as.data.frame(jack)
pearl=jack[,(21:ncol(jack))]
indx <- sapply(pearl, is.factor)
pearl[indx] <- lapply(pearl[indx], function(x) as.numeric(as.character(x)))
pearl[pearl==0]=NA
#remove peaks not present in at least 80% of samples in each group
require("plyr")
count(jack, c("RUN"))
count(jack, c("group"))
#captain=pearl[,-which(colSums(is.na(pearl[which(jack$group=="QC plasmamonster"),]))>=3)]
captain2=pearl[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=11 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=10 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=11)]
captain=captain2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=4)]
#---------START of adding IS------------
captain=cbind(captain, pearl[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
#captain=captain[,-as.numeric(ncol(captain))]
#---------END of adding IS------------
which(colSums(is.na(jack[which(jack$group=="CS"),]))>=5 & colSums(is.na(jack[which(jack$group=="PA"),]))>=13 & colSums(is.na(jack[which(jack$group=="PHT"),]))>=13 & colSums(is.na(jack[which(jack$group=="PPGL"),]))>=12)

test2=test[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=11 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=10 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=12 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=11)]
test2=test2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=4)]
#test2=test2[,-as.numeric(ncol(test2))]
test2=cbind(test2, test[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
#captain=captain[,-which(colSums(is.na(pearl[which(jack$group=="CS"),]))>=10 & colSums(is.na(pearl[which(jack$group=="PA"),]))>=26 & colSums(is.na(pearl[which(jack$group=="PPGL"),]))>=25 & colSums(is.na(pearl[which(jack$group=="PHT"),]))>=26)]
# calculate MV per sample, remove samples with >20% of total number of peaks
MV=rowSums(is.na(captain))
MV=MV/ncol(captain)
#captain=captain[-as.numeric(which(MV>0.3)),]
sparrow=jack#[-as.numeric(which(MV>0.3)),]
#TINT
captain[is.na(captain)]=0
sumcap=function(x){
  sum(captain[x,])
}
sums=sapply(1:as.numeric(nrow(captain)), sumcap)
captain[,as.numeric(ncol(captain))+1]=sums
tincap=captain/captain[,as.numeric(ncol(captain))]
tincap$V6346=NULL
captain=tincap
# PQN using HVs and using only peaks with no MV for calculating quotients
captain[captain==0]=NA
dauntless=captain
black=dauntless[which(sparrow$group=="HV" & sparrow$center=="Dresden (2)"),]
check1=apply(black, 2, median, na.rm=T)
check2=t(t(dauntless)/check1)
check3=apply(check2, 1, median, na.rm=T)
checkmate=captain/check3
checkmate[is.na(checkmate)]=0
#norrington=checkmate
#newcheck=checkmate*100000000
#checkmate=as.matrix(checkmate)
#test=sprintf("%.20f", checkmate)
#write.csv(test, "test.csv")
#write.csv(checkmate2, "test4.csv")
#test <- read_csv("Z:/NickB/all_neg2/BLANKS/results_XCMS_offline/test.csv")
#test4=as.matrix(test3)
#write.csv(newt, "test4.csv")
#test5 <- read_csv("Z:/NickB/all_neg2/BLANKS/results_XCMS_offline/test4.csv")
#Batch correction
checkmate.m=as.matrix(checkmate)
checkmate.m[which(checkmate.m==0)]=min(checkmate.m[which(checkmate.m>0)])
library(sva)
dat=t(checkmate.m)
batch=as.numeric(as.character(sparrow$BATCH))

combat_par=t(ComBat(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE, ref.batch = NULL))
combat_par[checkmate==0]=0
combat_par[combat_par<0]=0

checkmate=combat_par

#commodore=checkmate
#RSD filter:30% threshold
hay=function(x) {
  sd(checkmate[which(sparrow$group=="Qc Nick"),x], na.rm = T)
}
sdQc_Nick=sapply((1:(as.numeric(ncol(checkmate)))), hay)
bay=function(x) {
  mean(checkmate[which(sparrow$group=="Qc Nick"),x], na.rm = T)
}
meanQc_Nick=sapply((1:(as.numeric(ncol(checkmate)))), bay)
cvQc_Nick=sdQc_Nick/meanQc_Nick
cvQc_Nick[is.na(cvQc_Nick)]=0
hay=function(x) {
  sd(checkmate[which(sparrow$group=="QC plasmamonster"),x], na.rm = T)
}
sdQC_plasmamonster=sapply((1:(as.numeric(ncol(checkmate)))), hay)
bay=function(x) {
  mean(checkmate[which(sparrow$group=="QC plasmamonster"),x], na.rm = T)
}
meanQC_plasmamonster=sapply((1:(as.numeric(ncol(checkmate)))), bay)
cvQC_plasmamonster=sdQC_plasmamonster/meanQC_plasmamonster
#cvQC_plasmamonster[is.na(cvQC_plasmamonster)]=0.31
norrington=checkmate[,which(cvQC_plasmamonster<0.3)]

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=ncheckO.imputed$data
#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(sparrow$group=="QC plasmamonster"),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2)), 1))
#QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2))))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(sparrow[,1:20], final)
FEATURE_NAMES=test2[1:9,which(cvQC_plasmamonster<0.3)]
setwd("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT")
write.csv(final, file="EQP_neg_proc_newHV.csv")
#write.csv(cvQc_Nick, file="cvQc_Nick_nEQP_corr3.csv")
#write.csv(cvQC_plasmamonster, file="cvQC_plasmamonster_nEQP_corr3.csv")
final=final[,-(1:20)]
final=final[-c(6,9,27,33,36,38,66,84,85,184,187,213),]
library(readr)
EQP_neg_proc_METADATA3 <- read_csv("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/EQP_neg_proc_METADATA4.csv")

library(mixOmics)


X=final
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = EQP_neg_proc_METADATA3$GROUP, title = "ALL SAMPLES QTOF", ind.names = EQP_neg_proc_METADATA3$`ENSAT-HT ID`,legend = T,cex=4)


X=final[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST"),]
Y=EQP_neg_proc_METADATA3[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PATIENT SAMPLES QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
gr=Y$GROUP
gr[gr=="PA"]=gr[gr=="PPGL"]=gr[gr=="CS"]="AHT"
plotIndiv(pca.EQP, group = gr, title = "PATIENT SAMPLES QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plsda.EQP <- plsda(X, gr, ncomp = 20, scale=F)                  #number of components
perf.plsda.EQP <- perf(plsda.EQP, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
plsda.EQP <- plsda(X, gr, ncomp = perf.plsda.EQP$choice.ncomp[2,1], scale = F)


X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT"),]
gr=Y$`CENTER`
gr[gr=="FRPA"]=gr[gr=="GBGL"]=1
gr[gr=="ITTU"]=gr[gr=="GYDR"]=2
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = gr, title = "PHT PATIENT SAMPLES QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plsda.EQP <- plsda(X, gr, ncomp = 20, scale=F)                  #number of components
perf.plsda.EQP <- perf(plsda.EQP, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
plsda.EQP <- plsda(X, gr, ncomp = perf.plsda.EQP$choice.ncomp[2,1], scale = F)

X=final[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$`CENTER`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
Y=EQP_neg_proc_METADATA3[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$`CENTER`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
gr=Y$GROUP
gr[gr=="PA"]=gr[gr=="PPGL"]=gr[gr=="CS"]=1
gr[gr=="PHT"]=2
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = gr, title = "NO PHT OUTLIERS PATIENT SAMPLES QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plsda.EQP <- plsda(X, gr, ncomp = 20, scale=F)                  #number of components
perf.plsda.EQP <- perf(plsda.EQP, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
plsda.EQP <- plsda(X, gr, ncomp = perf.plsda.EQP$choice.ncomp[2,1], scale = F)


X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PA"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PA"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS PA QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS PA QTOF", legend = T,pch="sphere", style="3d")
WHYF=Y
ex=X

X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PPGL"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PPGL"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS PPGL QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS PPGL QTOF", legend = T,pch="sphere", style="3d")
WHYF=Y
ex=X

X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="CS"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="CS"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS CS QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS CS QTOF", legend = T,pch="sphere", style="3d")
WHYF=Y
ex=X

X=final[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PPGL" | EQP_neg_proc_METADATA3$GROUP=="CS" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
Y=EQP_neg_proc_METADATA3[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PPGL" | EQP_neg_proc_METADATA3$GROUP=="CS" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "NO PHT OUTLIERS PHT VS PA QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)

X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PPGL"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="PPGL"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS PPGL QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)

X=final[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PA" | EQP_neg_proc_METADATA3$GROUP=="CS" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
Y=EQP_neg_proc_METADATA3[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PA" | EQP_neg_proc_METADATA3$GROUP=="CS" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "NO PHT OUTLIERS PHT VS PPGL QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)

X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="CS"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT" | EQP_neg_proc_METADATA3$GROUP=="CS"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "PHT VS CS QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)

X=final[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PA" | EQP_neg_proc_METADATA3$GROUP=="PPGL" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
Y=EQP_neg_proc_METADATA3[-which(is.na(EQP_neg_proc_METADATA3$GROUP) | EQP_neg_proc_METADATA3$GROUP=="HV" | EQP_neg_proc_METADATA3$GROUP=="POST" | EQP_neg_proc_METADATA3$GROUP=="PA" | EQP_neg_proc_METADATA3$GROUP=="PPGL" | EQP_neg_proc_METADATA3$`CENTER ID`=="GBGL"  | (EQP_neg_proc_METADATA3$`CENTER ID`=="FRPA" & EQP_neg_proc_METADATA3$GROUP=="PHT")),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$GROUP, title = "NO PHT OUTLIERS PHT VS CS QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="CS")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="CS"),])))
    folds2=split(sample(which((WHYF$GROUP=="PHT"))), rep(1:f, length=nrow(X[which((WHYF$GROUP=="PHT")),])))
    folds=list(folds1, folds2)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$GROUP)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


allconf=Reduce('+', allconf)
((allconf[2,2]/(allconf[2,2]+allconf[2,1]))+(allconf[1,1]/(allconf[1,1]+allconf[1,2])))/2


vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

each.vip=function(x) {
  e.vip=function(y) {
    f.vip=function(z){
      allvip[[x]][[y]][z,which(allvip[[x]][[y]][z,min(as.numeric(1:ncol(allvip[[x]][[y]])))]>0)]
    }
    avip=sapply((1:as.numeric(ncol(X))), f.vip)
    avip=matrix(avip)
  }
  favip=lapply((1:f), e.vip)
}

allncvip=lapply((1:50), each.vip)
allv=function(x){
  do.call(cbind, allncvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checke2[order(checke2)]
checker=rowSums(is.na(allvvv))
names(checker)=colnames(X)
checker[order(checker)]
FEATURE_NAMES[,colnames(FEATURE_NAMES) %in% names(checker[which(checker<400)])]

MC.mitch=function(x) {
  fold.MC=function(y){
    mitch[[x]][[y]][[8]]
  }
  fold.MCs=lapply((1:f), fold.MC)
}

allMC=lapply((1:50), MC.mitch)

CC.mitch=function(x) {
  fold.CC=function(y){
    mitch[[x]][[y]][[9]]
  }
  fold.CCs=lapply((1:f), fold.CC)
}

allCC=lapply((1:50), CC.mitch)

cen=unmap(WHYF$`CENTER ID`)
enp.fac=as.factor(WHYF$GENDER)
levels(enp.fac)=c("0", "1")
levels(enp.fac)=as.numeric(levels(enp.fac))
gender=as.matrix(enp.fac)
gender=as.numeric(gender)
gender=as.matrix(gender)
colnames(gender)=c("gender")
enp.fac=as.factor(WHYF$GROUP)
levels(enp.fac)=c("0", "1")
levels(enp.fac)=as.numeric(levels(enp.fac))
group=as.matrix(enp.fac)
group=as.numeric(group)
group=as.matrix(group)
colnames(group)=c("GROUP")
enp.fac=as.factor(WHYF$`SAMPLE AGE`)
meisio=as.numeric(levels(enp.fac))
meisio[which(meisio<median(meisio, na.rm = T))]=0
meisio[which(meisio>=median(meisio, na.rm = T))]=1
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SA=as.matrix(enp.fac)
SA=as.numeric(SA)
SA=as.matrix(SA)
colnames(SA)=c("SA")
enp.fac=as.factor(WHYF$`PATIENT AGE`)
meisio=as.numeric(levels(enp.fac))
meisio[which(meisio<median(meisio, na.rm = T))]=0
meisio[which(meisio>=median(meisio, na.rm = T))]=1
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
AGE=as.matrix(enp.fac)
AGE=as.numeric(AGE)
AGE=as.matrix(AGE)
colnames(AGE)=c("AGE")
enp.fac=as.factor(WHYF$`QTOF No`)
meisio=as.numeric(levels(enp.fac))
meisio[which(meisio<median(meisio, na.rm = T))]=0
meisio[which(meisio>=median(meisio, na.rm = T))]=1
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SEQ=as.matrix(enp.fac)
SEQ=as.numeric(SEQ)
SEQ=as.matrix(SEQ)
colnames(SEQ)=c("SEQ")
enp.fac=as.factor(WHYF$BATCH)
meisio=as.numeric(levels(enp.fac))
meisio[which(meisio<median(meisio, na.rm = T))]=0
meisio[which(meisio>=median(meisio, na.rm = T))]=1
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
BATCH=as.matrix(enp.fac)
BATCH=as.numeric(BATCH)
BATCH=as.matrix(BATCH)
colnames(BATCH)=c("BATCH")
enp.fac=as.factor(WHYF$SAMPLEPREP)
meisio=as.numeric(levels(enp.fac))
meisio[which(meisio<median(meisio, na.rm = T))]=0
meisio[which(meisio>=median(meisio, na.rm = T))]=1
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SAMPLEPREP=as.matrix(enp.fac)
SAMPLEPREP=as.numeric(SAMPLEPREP)
SAMPLEPREP=as.matrix(SAMPLEPREP)
colnames(SAMPLEPREP)=c("SAMPLEPREP")
enp.fac=as.factor(WHYF$`CENTER ID`)
levels(enp.fac)=c("1", "2", "3", "4", "5", "6")
levels(enp.fac)=as.numeric(levels(enp.fac))
center=as.matrix(enp.fac)
center=as.numeric(center)
center=as.matrix(center)
colnames(center)=c("center")
clus=rep(1,as.numeric(nrow(WHYF)))
clus[which(WHYF$GROUP=="PHT"&WHYF$`CENTER ID`=="FRPA")]=0
clus[which(WHYF$`CENTER ID`=="GBGL")]=0
clus=as.factor(clus)
levels(clus)=as.numeric(levels(clus))
clus=as.matrix(clus)
clus=as.numeric(clus)
clus=as.matrix(clus)
colnames(clus)=c("clus")

X=as.matrix(X)
WHY=cbind(center,group, SA, gender, AGE, BATCH, SAMPLEPREP, SEQ, clus)
set.seed(1)
system.time(perm1<-PermuteTest(X, WHY, Interact = c("21", "23", "24", "25", "26", "27", "28", "29"), N=1000))

rm(mitch, par.output)
CSVPHT_Q1=list(all.BER, checker, checke2, allCC, allMC, perm1)

save.image("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/EQP_neg_fm_CSVPHT_newHV.RData")

X=final[which(EQP_neg_proc_METADATA3$GROUP=="PHT"),]
Y=EQP_neg_proc_METADATA3[which(EQP_neg_proc_METADATA3$GROUP=="PHT"),]
pca.EQP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = Y$`CENTER ID`, title = "PHT QTOF", ind.names = Y$`ENSAT-HT ID`,legend = T,cex=4)
plotIndiv(pca.EQP, group = Y$`CENTER ID`, title = "PHT QTOF", legend = T,pch="sphere", style="3d")
WHYF=Y
ex=X

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$`CENTER ID`=="FRPA")), rep(1:f, length=nrow(X[which(WHYF$`CENTER ID`=="FRPA"),])))
    folds2=split(sample(which(WHYF$`CENTER ID`=="GBGL")), rep(1:f, length=nrow(X[which(WHYF$`CENTER ID`=="GBGL"),])))
    folds3=split(sample(which((WHYF$`CENTER ID`=="GYDR"))), rep(1:f, length=nrow(X[which((WHYF$`CENTER ID`=="GYDR")),])))
    folds4=split(sample(which((WHYF$`CENTER ID`=="ITTU"))), rep(1:f, length=nrow(X[which((WHYF$`CENTER ID`=="ITTU")),])))
    folds=list(folds1, folds2, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=WHYF[train,]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack$`CENTER ID`)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "0", "1", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=WHYF[test,]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack$`CENTER ID`)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    enp.fac=as.factor(Y)
    levels(enp.fac)=c("0", "0", "1", "1")
    levels(enp.fac)=as.numeric(levels(enp.fac))
    group=as.matrix(enp.fac)
    group=as.numeric(group)
    group=as.matrix(group)
    colnames(group)=c("GROUP")
    Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    MCs=jack[-which(Prediction==Y),(1:14)]
    CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)


allconf=Reduce('+', allconf)
((allconf[2,2]/(allconf[2,2]+allconf[2,1]))+(allconf[1,1]/(allconf[1,1]+allconf[1,2])))/2


vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

each.vip=function(x) {
  e.vip=function(y) {
    f.vip=function(z){
      allvip[[x]][[y]][z,which(allvip[[x]][[y]][z,min(as.numeric(1:ncol(allvip[[x]][[y]])))]>0)]
    }
    avip=sapply((1:as.numeric(ncol(X))), f.vip)
    avip=matrix(avip)
  }
  favip=lapply((1:f), e.vip)
}

allncvip=lapply((1:50), each.vip)
allv=function(x){
  do.call(cbind, allncvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checke2[order(checke2)]
checker=rowSums(is.na(allvvv))
names(checker)=colnames(X)
checker[order(checker)]
FEATURE_NAMES[,colnames(FEATURE_NAMES) %in% names(checker[which(checker<300)])]

PHT_Q1=list(all.BER, checker, checke2)

save.image("D:/Office 2010 Toolkit and EZ-Activator v2.2.3/metabolomics hypertension/QTOF_EQP_ENSAT-HT/EQP_neg_fm_PHT_newHV.RData")

