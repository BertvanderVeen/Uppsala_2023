library(gllvm)
library(Momocs)
library(shapes)
#mouse molers
dat <- array(unlist(apodemus[[1]]),dim=c(nrow(apodemus[[1]][[1]]),2,length(apodemus[[1]])))
dat<-sapply(seq(dim(dat)[3]), function(x) matrix(unlist(apply(dat[,,x],1,c,simplify = F)),ncol=1))
#mod<-gllvm(dat,num.lv=2,family="gamma",n.init=5,row.eff="random")
mod<-gllvm(dat,num.lv=2,family="gamma",n.init=5,row.eff="fixed")
save(mod,file="/home/bertv/Dropbox/Uppsala/workshop/mice.RData")
#mod3<-gllvm(dat,num.lv=2,family="gamma",n.init=5)
#AIC(mod,mod2,mod3)
eta<-predict(mod,type="response")
coordPred<-array(unlist(sapply(seq(ncol(dat)), function(x) matrix((eta)[,x],ncol=2,byrow=T),simplify = F)),dim=c(nrow(dat)/2,2,ncol(dat)))

plot(NA,type="n",xlim=c(range(coordPred[,1,])),ylim=range(coordPred[,2,]), ylab="y",xlab="x") 
cols<-viridis::magma(ncol(dat), alpha = 0.5)
for(i in 1:ncol(dat)){
  lines(coordPred[,1,i],coordPred[,2,i],  col = cols)
}


#mosquito
dat <- array(unlist(mosquito[[1]]),dim=c(nrow(mosquito[[1]][[1]]),2,length(mosquito[[1]])))
dat<-sapply(seq(dim(dat)[3]), function(x) matrix(unlist(apply(dat[,,x],1,c,simplify = F)),ncol=1))
mod<-gllvm(dat,num.lv=2,family="gaussian",n.init=5,row.eff="random")
#save(mod,file="/home/bertv/Dropbox/Uppsala/workshop/mosquito.RData")
load("/home/bertv/Dropbox/Uppsala/workshop/mosquito.RData")
eta<-predict(mod)
coordPred<-array(unlist(sapply(seq(ncol(dat)), function(x) matrix((eta)[,x],ncol=2,byrow=T),simplify = F)),dim=c(nrow(dat)/2,2,ncol(dat)))

plot(NA,type="n",xlim=c(range(coordPred[,1,])),ylim=range(coordPred[,2,]), ylab="y",xlab="x") 
cols<-viridis::magma(ncol(dat), alpha = 0.5)
for(i in 1:ncol(dat)){
  lines(coordPred[,1,i],coordPred[,2,i],  col = cols)
}


#hearts
dat<-sapply(seq(length(hearts[[1]])), function(x) matrix(unlist(apply(hearts[[1]][[x]],1,c,simplify = F)),ncol=1))
mod<-gllvm(dat[,1:50],num.lv=2,family="gaussian",starting.val="random",n.init=5,row.eff="random")
eta<-predict(mod)

coordPred<-array(unlist(sapply(seq(ncol(mod$y)), function(x) matrix((eta)[,x],ncol=2,byrow=T),simplify = F)),dim=c(nrow(mod$y)/2,2,ncol(mod$y)))

plot(NA,type="n",xlim=c(range(coordPred[,1,])),ylim=range(coordPred[,2,])) 
for(i in 1:ncol(mod$y)){
  lines(coordPred[,1,i],coordPred[,2,i])
}
