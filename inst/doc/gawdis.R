## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(gawdis)

## -----------------------------------------------------------------------------
dummy$trait
ex1 <- gowdis(dummy$trait)
round(ex1, 3)##just to see only 3 decimals, enough
class(ex1)

## -----------------------------------------------------------------------------
distance.num2<-dist(dummy$trait[, "num2", drop=F])/max(dist(dummy$trait$num2))#note the drop=F not to loose the species names in the process#
round(distance.num2, 3)
distance.bin2<-dist(dummy$trait$bin2)
distance.bin2
distance.twotraits.byhand<-(distance.num2+distance.bin2)/2
distance.twotraits.gowdis<-gowdis(dummy$trait[, c("num2", "bin2")])
plot(distance.twotraits.gowdis, distance.twotraits.byhand)
abline(0, 1)

## -----------------------------------------------------------------------------
dummy$trait$fac2
gowdis(dummy$trait[, "fac2", drop=F])#notice that gowdis wants a matrix, not a vector, here a simple solution to solve this and keep species names

## -----------------------------------------------------------------------------
cor(distance.twotraits.gowdis, distance.num2)
cor(distance.twotraits.gowdis, distance.bin2)

## -----------------------------------------------------------------------------
distance.twotraits.byhand.weighted<-0.67*(distance.num2)+0.33*(distance.bin2)
distance.twotraits.gowdis.weighted<-gowdis(dummy$trait[, c("num2", "bin2")], w=c(2, 1))
plot(distance.twotraits.byhand.weighted, distance.twotraits.gowdis.weighted)
abline(0, 1)

## -----------------------------------------------------------------------------
cor(distance.twotraits.gowdis.weighted, distance.num2)
cor(distance.twotraits.gowdis.weighted, distance.bin2)

## -----------------------------------------------------------------------------
equalcont <- gawdis(dummy$trait[,c("num2", "bin2")])
equalcont
attr(equalcont,"correls")
attr(equalcont,"weights")

## ---- warning=FALSE-----------------------------------------------------------
ex1 <- gowdis(dummy$trait)
ex1.gaw1 <- gawdis(dummy$trait, w.type ="equal")
ex1.gaw2 <- gawdis(dummy$trait, w.type ="user", W=rep(1, 8))
plot(ex1, ex1.gaw1)
abline(0, 1)
plot(ex1, ex1.gaw2)
abline(0, 1)

## -----------------------------------------------------------------------------
analytical<-gawdis(dummy$trait[,c(2,4,6,8)], w.type ="analytic")# it is not needed to add the argument w.type, this is the approach used by default if not defined#
attr(analytical, "correls")
attr(analytical, "weights")
iterations<-gawdis(dummy$trait[,c(2,4,6,8)], w.type ="optimized", opti.maxiter=100)# here we used 'only' 100 iterations, to speed up the process and because with such few species this is likely more than enough#
attr(iterations, "correls")
attr(iterations, "weights")
plot(analytical, iterations)

## -----------------------------------------------------------------------------
dim(tussock$trait)
head(tussock$trait)
head(tussock$trait[, 3:7])
cor(tussock$trait[, 3:7], use = "complete")

## -----------------------------------------------------------------------------
tussock.trait<-tussock$trait[, c("height", "LDMC", "leafN","leafS", "leafP", "SLA", "seedmass", "raunkiaer", "pollination", "clonality", "growthform")]
tussock.trait.log<-tussock.trait#some traits needed log-tranformation, just creating a matrix to store the new data
tussock.trait.log$height<-log(tussock.trait$height)
tussock.trait.log$seedmass<-log(tussock.trait$seedmass)
tussock.trait.log$leafS<-log(tussock.trait$leafS)
colnames(tussock.trait.log)

## -----------------------------------------------------------------------------
#straightgowdis<-gowdis(tussock.trait.log)
straightgowdis.2<-gawdis(tussock.trait.log, w.type = "equal", silent = T)#we compute 'normal' gower with the new function because it provides more results
#plot(straightgowdis, straightgowdis.2)# if you want to check that the results are the same
cors.gow<-attr(straightgowdis.2,"correls")
cors.gow[12]<-mantel(straightgowdis.2, gowdis(tussock.trait.log[, 2:6]), na.rm=T)$statistic
names(cors.gow)[12]<-"leaves"
cors.gow

## -----------------------------------------------------------------------------
testpcoa<-dbFD(cailliez(straightgowdis.2), tussock$abun)###checking how many PCoA axes are retained
pcoaaxes<-dudi.pco(cailliez(straightgowdis.2), scannf = FALSE, nf = 11)
gowdis.PCoA<-dist(pcoaaxes$li)
sum(pcoaaxes$eig[1:11])/sum(pcoaaxes$eig)#how much variability the axes explain
#contribution of traits on the combined dissimilarity, done by hand#
cors.pcoa<-vector()
for(i in 1:dim(tussock.trait.log)[2]){
  cors.pcoa[i]<-mantel(gowdis.PCoA, gowdis(as.data.frame(tussock.trait.log[, i])), na.rm=T)$statistic
}
cors.pcoa[12]<-mantel(gowdis.PCoA, gowdis(tussock.trait.log[, 2:6]), na.rm=T)$statistic
names(cors.pcoa)<-c(colnames(tussock.trait.log), "leaves")#contributions of traits to the overall multi-trait dissimilarity
cors.pcoa

## -----------------------------------------------------------------------------
gaw.nogroups<-gawdis(tussock.trait.log, w.type = "optimized", opti.maxiter = 200)#there are NAs so the iteration approach is the only possible
cors.gaw<-attr(gaw.nogroups,"correls")
cors.gaw[12]<-mantel(gaw.nogroups, gowdis(as.data.frame(tussock.trait.log[, 2:6])), na.rm=T)$statistic
names(cors.gaw)[12]<-"leaves"
cors.gaw

## -----------------------------------------------------------------------------
colnames(tussock.trait.log)
gaw.groups<-gawdis(tussock.trait.log, w.type = "optimized", opti.maxiter = 200, groups.weight=T, groups = c(1,2, 2, 2, 2, 2, 3, 4, 5, 6, 7))#there are NAs so the iteration approach is the only possible
cors.gaw.gr<-attr(gaw.groups,"correls")
cors.gaw.gr[12]<-attr(gaw.groups,"group.correls")[2]
names(cors.gaw.gr)[12]<-"leaves"
cors.gaw.gr

## -----------------------------------------------------------------------------
bodysize<-c(10, 20, 30, 40, 50, NA, 70)
carnivory<-c(1, 1, 0, 1, 0,1, 0)
red<-c(1, 0, 0.5, 0, 0.2, 0, 1)
yellow<-c(0, 1, 0, 0, 0.3, 1, 0)
blue<-c(0, 0, 0.5,1, 0.5, 0, 0)
colors.fuzzy<-cbind(red, yellow, blue)
names(bodysize)<-paste("sp", 1:7, sep="")
names(carnivory)<-paste("sp", 1:7, sep="")
rownames(colors.fuzzy)<-paste("sp", 1:7, sep="")
tall<-as.data.frame(cbind(bodysize, carnivory, colors.fuzzy))
tall

## -----------------------------------------------------------------------------
round(gowdis(tall[, 3:5]), 3)

## -----------------------------------------------------------------------------
round(gowdis(tall[, 3:5])/max(gowdis(tall[, 3:5])), 3)

## -----------------------------------------------------------------------------
dissim.bodysize<-gowdis(tall[, "bodysize", drop=F])
dissim.carnivory<-gowdis(tall[, "carnivory", drop=F])
dissim.colour<-gowdis(tall[, 3:5])/max(gowdis(tall[, 3:5]))
dall<-list(as.matrix(dissim.bodysize), as.matrix(dissim.carnivory), as.matrix(dissim.colour))
mean.dissim.all<-as.dist(apply(simplify2array(dall), c(1, 2), mean, na.rm=T), 2)
round(mean.dissim.all, 3)

## -----------------------------------------------------------------------------
gawdis(tall, w.type="equal", groups =c(1, 2, 3, 3, 3), fuzzy=TRUE)

