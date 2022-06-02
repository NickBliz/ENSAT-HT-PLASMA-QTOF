setwd("D:/MEGA/Publication 4/QTOF/ALL_NEG/BLANKS")
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
EQP_av_m[c(56,57,173,174,294,295,406,407,518,519,635,636),]=tzelis[c(56,57,173,174,294,295,406,407,518,519,635,636),1:ncol(tzelis)]
#EQP_av_m=cbind(EQPc[,(1:20)], EQP_av_m)
#EQP_av_m_m=EQP_av_m[-which(EQP_av_m$`F/B`=="backwards"),]

EQP_av_m_m=cbind(EQPc[,(1:20)], EQP_av_m)
test=jel[which(row.names(EQP_av_m_m) %in% row.names(jel)),]
EQP_bp=t(EQP_av_m_m) #REMEMBER! <jel> has your variable names!

jack=t(EQP_bp)
jack=as.data.frame(jack)
pearl=jack[,(21:ncol(jack))]
indx <- sapply(pearl, is.factor)
pearl[indx] <- lapply(pearl[indx], function(x) as.numeric(as.character(x)))
pearl[pearl==0]=NA
#REMOVE SAMPLES SO THAT CORRECTION IS FEASIBLE: 355:370, 443:458  (RUN2) Qc Nick - Patient plasma 1PPLWW-1024
jack=jack[-c(355:370,372,373,440, 441, 443:458),]
pearl=pearl[-c(355:370,372,373,440, 441, 443:458),]
#remove peaks not present in at least 80% of samples in each group
require("plyr")
count(jack, c("RUN"))
count(jack, c("group"))
#captain=pearl[,-which(colSums(is.na(pearl[which(jack$group=="QC plasmamonster"),]))>=3)]
#captain=pearl[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=112 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=72 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=122 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=120 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=120 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=114)]
#jel=jel[,-(which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=112 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=72 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=122 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=120 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=120 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=114))]
#captain=captain[,-which(colSums(is.na(pearl[which(jack$group=="CS"),]))>=10 & colSums(is.na(pearl[which(jack$group=="PA"),]))>=26 & colSums(is.na(pearl[which(jack$group=="PPGL"),]))>=25 & colSums(is.na(pearl[which(jack$group=="PHT"),]))>=26)]
captain2=pearl[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=22 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=13 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=22)]
captain=captain2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=8)]
captain=cbind(captain, pearl[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
test2=test[,-which(colSums(is.na(pearl[which(jack$RUN=="run1"),]))>=22 | colSums(is.na(pearl[which(jack$RUN=="run2"),]))>=13 | colSums(is.na(pearl[which(jack$RUN=="run3"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run4"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run5"),]))>=24 | colSums(is.na(pearl[which(jack$RUN=="run6"),]))>=22)]
test2=test2[,-which(colSums(is.na(captain2[which(jack$group=="QC plasmamonster"),]))>=8)]
test2=cbind(test2, test[,c(as.numeric(ncol(pearl)), as.numeric(ncol(pearl))-1, as.numeric(ncol(pearl))-2)])
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
tincap$V6399=NULL
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


