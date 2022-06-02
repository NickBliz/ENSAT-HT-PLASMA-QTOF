X=EQP_NEG
row.names(X)=EQP_METADATA_NEG$`sample name 2`
X=X[which(EQP_METADATA_NEG$GROUP=="PHT"),]
WHYF=data.frame(EQP_METADATA_NEG[which(EQP_METADATA_NEG$GROUP=="PHT"),])
pca.EQP = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = WHYF$CENTER.ID, title = "NEG QTOF", ind.names = F,legend = T,cex=1, pch = "sphere", style = "3d")
ex=X
clus=Y=WHYF$CENTER.ID
#clus[which(clus=="FRPA1" | clus=="GBGL2")]=1
#clus[-which(clus=="1")]=0

plotIndiv(pca.EQP, group = clus, title = "NEG QTOF", ind.names = F,legend = T,cex=1, pch = "sphere", style = "3d")

f=6
n=50

hay=function(x) {
  {
    #folds1=split(sample(which(clus=="0")), rep(1:f, length=nrow(X[which(clus=="0"),])))
    #folds2=split(sample(which(clus=="1")), rep(1:f, length=nrow(X[which(clus=="1"),])))
    folds1=split(sample(which(clus=="FRPA1")), rep(1:f, length=nrow(X[which(clus=="FRPA1"),])))
    folds2=split(sample(which(clus=="GBGL2")), rep(1:f, length=nrow(X[which(clus=="GBGL2"),])))
    folds3=split(sample(which(clus=="GYDR")), rep(1:f, length=nrow(X[which(clus=="GYDR"),])))
    folds4=split(sample(which(clus=="ITTU3")), rep(1:f, length=nrow(X[which(clus=="ITTU3"),])))
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
    jack=clus[train]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    #enp.fac=as.factor(Y)
    #levels(enp.fac)=as.numeric(levels(enp.fac))
    #group=as.matrix(enp.fac)
    #group=as.numeric(group)
    #group=as.matrix(group)
    #colnames(group)=c("GROUP")
    #Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    plsda.train <- plsda(X, Y, ncomp = ncomp, scale = F)
    jack=clus[test]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    #enp.fac=as.factor(Y)
    #levels(enp.fac)=as.numeric(levels(enp.fac))
    #group=as.matrix(enp.fac)
    #group=as.numeric(group)
    #group=as.matrix(group)
    #colnames(group)=c("GROUP")
    #Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
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

each.0=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.0=lapply((1:n), each.0)
all.0=unlist(all.0)

each.1=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.1=lapply((1:n), each.1)
all.1=unlist(all.1)

each.2=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.2=lapply((1:n), each.2)
all.2=unlist(all.2)

each.3=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.3=lapply((1:n), each.3)
all.3=unlist(all.3)

vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]][,as.numeric(ncol(mitch[[x]][[y]][[3]]))]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

allv=function(x){
  do.call(cbind, allvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checker=checke2[order(checke2)]
checker[which(checker>1)]

#C1VC21=list(all.BER, checker, checke2)
C1VC22=list(allconf, checker, checke2)


hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]], folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=clus[train]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=clus[test]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(splsda.train)
    #MCs=jack[-which(Prediction==Y),(1:14)]
    #CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
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


each.FRPA1=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.FRPA1=lapply((1:n), each.FRPA1)
all.FRPA1=unlist(all.FRPA1)

each.GBGL2=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.GBGL2=lapply((1:n), each.GBGL2)
all.GBGL2=unlist(all.GBGL2)

each.GYDR=function(x){
  allconf[[x]][3,3]/sum(allconf[[x]][3,])
}

all.GYDR=lapply((1:n), each.GYDR)
all.GYDR=unlist(all.GYDR)

each.ITTU3=function(x){
  allconf[[x]][4,4]/sum(allconf[[x]][4,])
}

all.ITTU3=lapply((1:n), each.ITTU3)
all.ITTU3=unlist(all.ITTU3)


plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
system.time(perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = f, auc = F, progressBar = T, nrepeat = n, dist="mahalanobis.dist"))
#600sec
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
system.time(tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = f, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = n, scale = F, cpus = 6))
#80sec
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_C1VC221=test.predict$B.hat[,,ncomp]
rich=B.hat_C1VC221#[,1]
rich=rich[which(rowSums(rich)>0 | rowSums(rich)<0),]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)

C1VC221=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
C1VC222=list(allconf, all.BER, all.FRPA1, all.GBGL2 , all.GYDR , all.ITTU3)








X=EQP_NEG
row.names(X)=EQP_METADATA_NEG$`sample name 2`
X=X[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),]
WHYF=data.frame(EQP_METADATA_NEG[-which(EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),])
pca.EQP = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = WHYF$CENTER.ID, title = "NEG QTOF", ind.names = F,legend = T,cex=1, pch = "sphere", style = "3d")
ex=X
clus=as.numeric(as.character(WHYF$SAMPLE.AGE))
clus[which(clus<median(clus))]=1
clus[-which(clus==1)]=0
Y=clus
plotIndiv(pca.EQP, group = clus, title = "NEG QTOF", ind.names = F,legend = T,cex=1, pch = "sphere", style = "3d")

f=6
n=50

hay=function(x) {
  {
    folds1=split(sample(which(clus=="0")), rep(1:f, length=nrow(X[which(clus=="0"),])))
    folds2=split(sample(which(clus=="1")), rep(1:f, length=nrow(X[which(clus=="1"),])))
    #folds1=split(sample(which(clus=="FRPA1")), rep(1:f, length=nrow(X[which(clus=="FRPA1"),])))
    #folds2=split(sample(which(clus=="GBGL2")), rep(1:f, length=nrow(X[which(clus=="GBGL2"),])))
    #folds3=split(sample(which(clus=="GYDR")), rep(1:f, length=nrow(X[which(clus=="GYDR"),])))
    #folds4=split(sample(which(clus=="ITTU3")), rep(1:f, length=nrow(X[which(clus=="ITTU3"),])))
    folds=list(folds1, folds2)#, folds3, folds4)
  }
  return(folds)
}

set.seed(1)
newF=lapply(1:50, hay)



hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])#, folds[[3]][[i]], folds[[4]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=clus[train]
    ENP=ex[train,]
    dataset_train=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    #enp.fac=as.factor(Y)
    #levels(enp.fac)=as.numeric(levels(enp.fac))
    #group=as.matrix(enp.fac)
    #group=as.numeric(group)
    #group=as.matrix(group)
    #colnames(group)=c("GROUP")
    #Y=as.factor(group)
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    ncomp=min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))
    plsda.train <- plsda(X, Y, ncomp = ncomp, scale = F)
    jack=clus[test]
    ENP=ex[test,]
    dataset_test=ENP
    group=as.factor(jack)
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    #enp.fac=as.factor(Y)
    #levels(enp.fac)=as.numeric(levels(enp.fac))
    #group=as.matrix(enp.fac)
    #group=as.numeric(group)
    #group=as.matrix(group)
    #colnames(group)=c("GROUP")
    #Y=as.factor(group)
    X=as.matrix(X)
    test.predict <- predict(plsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
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

each.0=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.0=lapply((1:n), each.0)
all.0=unlist(all.0)

each.1=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.1=lapply((1:n), each.1)
all.1=unlist(all.1)


vip.mitch=function(x) {
  fold.vip=function(y){
    mitch[[x]][[y]][[3]][,as.numeric(ncol(mitch[[x]][[y]][[3]]))]
  }
  fold.vips=lapply((1:f), fold.vip)
}

allvip=lapply((1:50), vip.mitch)

allv=function(x){
  do.call(cbind, allvip[[x]])
}

allvv=lapply(1:50, allv)
allvvv=do.call(cbind, allvv)
mode(allvvv)="numeric"
library(matrixStats)
checke2=rowMedians(allvvv, na.rm = T)
names(checke2)=colnames(X)
checker=checke2[order(checke2)]
checker[which(checker>1)]
checkk=checker[which(checker>1)]

C1VC23=list(allconf, checker, checke2)

hay=function(x) {
  folds=newF[[x]]
  say=function(i){
    test=c(folds[[1]][[i]], folds[[2]][[i]])
    train <- setdiff(1:nrow(WHYF), test)
    jack=clus[train]
    ENP=ex[train,]
    dataset_train=ENP
    group=jack
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    plsda.train <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 5, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 5, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.train$choice.keepX)
    jack=clus[test]
    ENP=ex[test,]
    dataset_test=ENP
    group=jack
    ex1=as.matrix(ENP)
    X=ex1
    Y=group
    X=as.matrix(X)
    test.predict <- predict(splsda.train, X, dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, ncomp]                         #number of components
    confusion.mat = get.confusion_matrix(truth = Y, predicted = Prediction)
    BER.res = get.BER(confusion.mat)
    rich=vip(plsda.train)
    #MCs=jack[-which(Prediction==Y),(1:14)]
    #CCs=jack[which(Prediction==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, Prediction, Y)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
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

each.0=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.0=lapply((1:n), each.0)
all.0=unlist(all.0)

each.1=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.1=lapply((1:n), each.1)
all.1=unlist(all.1)






plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
system.time(perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = f, auc = F, progressBar = T, nrepeat = n, dist="mahalanobis.dist"))
#600sec
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
system.time(tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = f, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = n, scale = F, cpus = 6))
#80sec
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_C1VC2311=test.predict$B.hat[,,ncomp]
rich=B.hat_C1VC2311[,1]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)



C1VC2311=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
C1VC2312=list(allconf, all.BER, all.0, all.1)
