X=EQP_NEG[-which(EQP_METADATA_NEG$GROUP=="CS" | EQP_METADATA_NEG$GROUP=="PPGL" | EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),]
pca.EQP = pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.EQP, group = EQP_METADATA_NEG$GROUP[-which(EQP_METADATA_NEG$GROUP=="CS" | EQP_METADATA_NEG$GROUP=="PPGL" | EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV")], title = "NEG QTOF", ind.names = F,legend = T,cex=4, pch = "sphere", style = "3d")
plotIndiv(pca.EQP, group = EQP_METADATA_NEG$`CENTER ID`[-which(EQP_METADATA_NEG$GROUP=="CS" | EQP_METADATA_NEG$GROUP=="PPGL" | EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV")], title = "NEG QTOF", ind.names = F,legend = T,cex=4, pch = "sphere", style = "3d")

WHYF=data.frame(EQP_METADATA_NEG[-which(EQP_METADATA_NEG$GROUP=="CS" | EQP_METADATA_NEG$GROUP=="PPGL" | EQP_METADATA_NEG$GROUP=="X" | EQP_METADATA_NEG$GROUP=="Qc Nick" | EQP_METADATA_NEG$GROUP=="QC plasmamonster" | EQP_METADATA_NEG$GROUP=="HV"),])
sa=as.numeric(WHYF$`SAMPLE.AGE`)
kruskal.test(sa~WHYF$CENTER.ID)
#p-value < 2.2e-16
sa[which(sa<median(sa))]="<median"
sa[-which(sa=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = sa, legend = T, pch=c(19), cex = 5, title="(c): PCA SAMPLE AGE")
cor.test(pI$df$x, sa, method = "kendall")
#p-value = 2.406e-14
cor.test(pI$df$x, sa, method = "spearman")
#p-value < 2.2e-16
cor.test(pI$df$y, sa, method = "kendall")
#p-value = 0.5576
cor.test(pI$df$y, sa, method = "spearman")
#p-value = 0.5997
SEQ=as.numeric(WHYF$QTOF.No)
kruskal.test(SEQ~WHYF$CENTER.ID)
#p-value = 0.02776
SEQ[which(SEQ<median(SEQ))]="<median"
SEQ[-which(SEQ=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = SEQ, legend = T, pch=c(19), cex = 5, title="(c): PCA RUN ORDER")
cor.test(pI$df$x, SEQ, method = "kendall")
#p-value = 0.8754
cor.test(pI$df$x, SEQ, method = "spearman")
#p-value = 0.9223
cor.test(pI$df$y, SEQ, method = "kendall")
#p-value = 0.8722
cor.test(pI$df$y, SEQ, method = "spearman")
#p-value = 0.7872
RUN=as.numeric(WHYF$BATCH)
kruskal.test(RUN~WHYF$CENTER.ID)
#p-value = 0.8465
RUN[which(RUN<median(RUN))]="<median"
RUN[-which(RUN=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = RUN, legend = T, pch=c(19), cex = 5, title="(c): PCA BATCH")
cor.test(pI$df$x, RUN, method = "kendall")
#p-value = 0.9045
cor.test(pI$df$x, RUN, method = "spearman")
#p-value = 0.9266
cor.test(pI$df$y, RUN, method = "kendall")
#p-value = 0.7003
cor.test(pI$df$y, RUN, method = "spearman")
#p-value = 0.6896


X1=X
colnames(X1)=FEATURE_NAMES[2,]
pca.EQP = pca(X1, ncomp = 10, center = TRUE, scale = FALSE)
plotVar(pca.EQP, cex=3)

#final=cbind(sparrow[,1:20], final)

WHYF=data.frame(EQP_METADATA_NEG)
WHYF=WHYF[which(WHYF$GROUP=="PA" | WHYF$GROUP=="PHT"),]
ex=X
Y=WHYF$GROUP
length(which(Y=="PA"))
length(which(Y=="PHT"))

f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
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
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = F)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
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

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)





plsda.select <- plsda(X, Y, ncomp = 10, scale=F)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = F)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = F, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_PAVPHT1=test.predict$B.hat[,,ncomp]
rich=B.hat_PAVPHT1[,1]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)



PAVPHT1=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
PAVPHT2=list(allconf, all.BER, all.PA, all.PHT)


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
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 7, type.measure = "class", standardize = F)
    #fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F, lambda = cvfit1$lambda.min)
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
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)





X=ex
Y[which(Y=="PA")]="1"
Y[which(Y=="PHT")]="0"
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=FALSE)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)

PAVPHT3=list(cvfit1, rich)
PAVPHT4=list(allconf, all.BER, all.PA, all.PHT)


sex=as.factor(WHYF$GENDER)
levels(sex)=c("0", "1")
levels(sex)=as.numeric(levels(sex))
sex=as.numeric(sex)
sex[sex==1]=0
sex[sex==2]=1
age=as.numeric(WHYF$PATIENT.AGE)
ex=cbind(ex, sex, age)

X=ex
Y=WHYF$GROUP
#colnames(X)=round(as.numeric(colnames(X)),3)
pca.EQP=mixOmics::pca(X, scale = T, center = T, ncomp = 10)
plotIndiv(pca.EQP, group = Y, legend = T, style = "3d", pch="sphere")
plotIndiv(pca.EQP, group = WHYF$CENTER.ID, legend = T, style = "3d", pch="sphere")
length(which(Y=="PA"))
length(which(Y=="PHT"))

sa=as.numeric(WHYF$`SAMPLE.AGE`)
kruskal.test(sa~WHYF$CENTER.ID)
#p-value < 2.2e-16
sa[which(sa<median(sa))]="<median"
sa[-which(sa=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = sa, legend = T, pch=c(19), cex = 5, title="(c): PCA SAMPLE AGE")
cor.test(pI$df$x, sa, method = "kendall")
#p-value = 7.049e-10
cor.test(pI$df$x, sa, method = "spearman")
#p-value = 1.392e-10
cor.test(pI$df$y, sa, method = "kendall")
#p-value = 5.2e-06
cor.test(pI$df$y, sa, method = "spearman")
#p-value = 1.081e-06
SEQ=as.numeric(WHYF$QTOF.No)
kruskal.test(SEQ~WHYF$CENTER.ID)
#p-value = 0.02776
SEQ[which(SEQ<median(SEQ))]="<median"
SEQ[-which(SEQ=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = SEQ, legend = T, pch=c(19), cex = 5, title="(c): PCA RUN ORDER")
cor.test(pI$df$x, SEQ, method = "kendall")
#p-value = 0.2999
cor.test(pI$df$x, SEQ, method = "spearman")
#p-value = 0.3238
cor.test(pI$df$y, SEQ, method = "kendall")
#p-value = 0.5372
cor.test(pI$df$y, SEQ, method = "spearman")
#p-value = 0.5514
RUN=as.numeric(WHYF$BATCH)
kruskal.test(RUN~WHYF$CENTER.ID)
#p-value = 0.8465
RUN[which(RUN<median(RUN))]="<median"
RUN[-which(RUN=="<median")]=">=median"
pI=plotIndiv(pca.EQP, group = RUN, legend = T, pch=c(19), cex = 5, title="(c): PCA BATCH")
cor.test(pI$df$x, RUN, method = "kendall")
#p-value = 0.4999
cor.test(pI$df$x, RUN, method = "spearman")
#p-value = 0.495
cor.test(pI$df$y, RUN, method = "kendall")
#p-value = 0.9934
cor.test(pI$df$y, RUN, method = "spearman")
#p-value = 0.9583


f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
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
    plsda.train <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    perf.plsda.train <- perf(plsda.train, validation = "Mfold", folds = 7, auc = F, progressBar = T, nrepeat = 1, dist="mahalanobis.dist")
    set.seed(1)
    tune.splsda.train <- tune.splsda(X, Y, ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), validation = "Mfold", folds = 7, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 1, scale = T)
    ncomp=ceiling(min(which(tune.splsda.train$error.rate==min(tune.splsda.train$error.rate)))/nrow(tune.splsda.train$error.rate))
    splsda.train <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.train$choice.keepX)
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

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)




X=ex
plsda.select <- plsda(X, Y, ncomp = 10, scale=T)                  #number of components
set.seed(1)  #for keeping perf the same as with lapply
perf.plsda.select <- perf(plsda.select, validation = "Mfold", folds = 8, auc = F, progressBar = T, nrepeat = 50, dist="mahalanobis.dist")
ncomp=perf.plsda.select$choice.ncomp[2,]
set.seed(1)
tune.splsda.select <- tune.splsda(X, Y, ncomp = ncomp, validation = "Mfold", folds = 8, progressBar = T, dist = "mahalanobis.dist", measure = "BER", nrepeat = 50, scale = T)
ncomp=tune.splsda.select$choice.ncomp$ncomp
splsda.select <- splsda(X, Y, ncomp = ncomp, scale = T, keepX = tune.splsda.select$choice.keepX)
rich=vip(splsda.select)
rich[which(rowSums(rich)>0),]
test.predict <- predict(splsda.select, X, dist = "mahalanobis.dist")
B.hat_PAVPHT1=test.predict$B.hat[,,ncomp]
rich=B.hat_PAVPHT1[,1]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)



PAVPHT5=list(plsda.select, splsda.select, perf.plsda.select, tune.splsda.select, rich)
PAVPHT6=list(allconf, all.BER, all.PA, all.PHT)


f=8
n=50

hay=function(x) {
  {
    folds1=split(sample(which(WHYF$GROUP=="PA")), rep(1:f, length=nrow(X[which(WHYF$GROUP=="PA"),])))
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
    fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = T)                  #number of components
    set.seed(1)  #for keeping perf the same as with lapply
    cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 7, type.measure = "class")
    #fit = glmnet(X, Y, family = "binomial", alpha=1, standardize = F, lambda = cvfit1$lambda.min)
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
    test.predict <- predict(cvfit1, X, type = "class", s = "lambda.min")     #change distance
    confusion.mat = get.confusion_matrix(truth = Y, predicted = test.predict)
    BER.res = get.BER(confusion.mat)
    sensitivity=confusion.mat[1,1]/sum(confusion.mat[1,])
    specificity=confusion.mat[2,2]/sum(confusion.mat[2,])
    rich=coef(cvfit1, s = "lambda.min")
    MCs=jack[-which(test.predict==Y),(1:14)]
    CCs=jack[which(test.predict==Y),(1:14)]
    conrich=list(confusion.mat, BER.res, rich, dataset_train, dataset_test, test.predict, Y, MCs, CCs)
    return(conrich)
  }
  check=lapply(1:f, say)
}

#set.seed(1)
system.time( par.output <- mclapply.hack( 1:n,hay))

mitch=par.output

conf.mitch=function(x) {
  mitch[[x]][[8]][[1]]+mitch[[x]][[7]][[1]]+mitch[[x]][[6]][[1]]+mitch[[x]][[5]][[1]]+mitch[[x]][[4]][[1]]+mitch[[x]][[3]][[1]]+mitch[[x]][[2]][[1]]+mitch[[x]][[1]][[1]]
}

allconf=lapply((1:n), conf.mitch)

each.BER=function(x){
  get.BER(allconf[[x]])
}

all.BER=lapply((1:n), each.BER)
all.BER=unlist(all.BER)

each.PA=function(x){
  allconf[[x]][1,1]/sum(allconf[[x]][1,])
}

all.PA=lapply((1:n), each.PA)
all.PA=unlist(all.PA)

each.PHT=function(x){
  allconf[[x]][2,2]/sum(allconf[[x]][2,])
}

all.PHT=lapply((1:n), each.PHT)
all.PHT=unlist(all.PHT)





X=as.matrix(ex)
Y[which(Y=="PA")]="1"
Y[which(Y=="PHT")]="0"
set.seed(1)  #for keeping perf the same as with lapply
cvfit1 = cv.glmnet(X, Y, family = "binomial", alpha=1, nfolds = 8, type.measure = "class", standardize=T)
rich=coef(cvfit1, s = "lambda.min")
rich=rich[-1,]
rich=rich[which(rich>0 | rich<0)]
rich=as.matrix(rich)
rich
rich_names=FEATURE_NAMES[,which(FEATURE_NAMES[1,] %in% row.names(rich))]
rich=cbind(t(rich_names), rich)

PAVPHT7=list(cvfit1, rich)
PAVPHT8=list(allconf, all.BER, all.PA, all.PHT)

