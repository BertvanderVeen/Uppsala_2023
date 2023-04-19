ggEllips <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2), predict.region = FALSE, level =0.95,
jitter = FALSE, jitter.amount = 0.2, s.colors = 1, s.cex = 1.2, symbols = FALSE, cex.spp = 0.7, spp.colors = "blue", arrow.scale = 0.8, arrow.spp.scale = 0.8, arrow.ci = TRUE, arrow.lty = "solid", spp.arrows = NULL, spp.arrows.lty = "dashed", cex.env = 0.7, lab.dist = 0.1, lwd.ellips = 0.5, col.ellips = 4, lty.ellips = 1, type = NULL, rotate  = FALSE, ...){

if (!any(class(object) %in% "gllvm"))
  stop("Class of the object isn't 'gllvm'.")
if(is.null(spp.arrows)){
  if(object$quadratic!=FALSE){
    spp.arrows <- TRUE
  }else{
    spp.arrows <- FALSE
  }
}

arrow.scale <- abs(arrow.scale)
a <- jitter.amount
n <- NROW(object$y)
p <- NCOL(object$y)
num.lv <- object$num.lv
num.lv.c <- object$num.lv.c
num.RR <- object$num.RR
quadratic <- object$quadratic
gr_par_list <- list(...)
# If both scales are not given, use MASS::eqscplot
if(("ylim" %in% names(gr_par_list)) & ("xlim" %in% names(gr_par_list))){
  plotfun <- plot
} else {
  plotfun <- MASS::eqscplot
}

if (!is.null(ind.spp)) {
  ind.spp <- min(c(p, ind.spp))
} else {
  ind.spp <- p
}
if(length(spp.colors)==1){
  spp.colors <- rep(spp.colors,p)
}else if(length(spp.colors)!=p){
  stop("spp.colors needs to be of length p or 1.")
}
if ((num.lv+(num.lv.c+num.RR)) == 0)
  stop("No latent variables to plot.")

if (is.null(rownames(object$params$theta)))
  rownames(object$params$theta) = paste("V", 1:p)

if((num.lv.c+num.RR)==0&is.null(type)){
  type <- "residual"
}else if(is.null(type)){
  if(num.lv.c==0){
    type <- "marginal"
  }else{
    type <- "conditional"
  }
}
if(!is.null(type)){

}
# This must be done, otherwise the scale of the ordination is non informative if the scale of params$theta () differ drastically:
if(type == "residual"|num.lv>0){
  # First create an index for the columns with unconstrained LVs
  which.scl <- NULL
  if(num.lv>0){
    which.scl <- (num.lv.c+num.RR+1):(num.lv.c+num.RR+num.lv)
  }

  if(quadratic!=FALSE){
    which.scl <- c(which.scl, which.scl+num.lv.c+num.RR+num.lv)
  }
  # Add indices of constrained LVs if present
  if(num.lv.c>0&type=="residual"){
    which.scl <- c(which.scl, 1:num.lv.c)
    if(quadratic!=FALSE){
      which.scl <- c(which.scl, 1:num.lv.c+num.lv.c+num.RR+num.lv)
    }
  }
  which.scl <- sort(which.scl)
  # Do the scaling
  if(quadratic==FALSE){
    object$params$theta[,which.scl] <- object$params$theta[,which.scl, drop=FALSE]%*%diag(object$params$sigma.lv, nrow = length(object$params$sigma.lv), ncol = length(object$params$sigma.lv))
  }else{
    sig <- diag(c(object$params$sigma.lv,object$params$sigma.lv^2))
    object$params$theta[,which.scl] <- object$params$theta[,which.scl,drop=FALSE]%*%sig
  }

}

lv <- gllvm::getLV(object, type = type)

  do_svd <- svd(lv)
  if(num.lv.c>0|(num.RR+num.lv)>0){
    do_svd$v <- svd(gllvm::getLV(object))$v
  }
  if(!rotate){
    do_svd$v <- diag(ncol(do_svd$v))
  }
  svd_rotmat_sites <- do_svd$v
  svd_rotmat_species <- do_svd$v



  choose.lvs <- lv
  if(quadratic == FALSE){choose.lv.coefs <- object$params$theta}else{choose.lv.coefs<-gllvm::optima(object,sd.errors=F)}

  if(spp.arrows){
    lvth <- max(abs(choose.lvs))
    idx <- choose.lv.coefs>(-lvth)&choose.lv.coefs<lvth
  }else{
    idx <- matrix(TRUE,ncol=num.lv+num.lv.c+num.RR,nrow=p)
  }

  bothnorms <- vector("numeric",ncol(choose.lv.coefs))
  for(i in 1:ncol(choose.lv.coefs)){
    bothnorms[i] <- sqrt(sum(choose.lvs[,i]^2)) * sqrt(sum(choose.lv.coefs[idx[,i],i]^2))
  }

  scaled_cw_sites <- t(t(choose.lvs) / sqrt(colSums(choose.lvs^2)) * (bothnorms^alpha))
  scaled_cw_species <- choose.lv.coefs
  for(i in 1:ncol(scaled_cw_species)){
    scaled_cw_species[,i] <- choose.lv.coefs[,i] / sqrt(sum(choose.lv.coefs[idx[,i],i]^2)) * (bothnorms[i]^(1-alpha))
  }

  choose.lvs <- scaled_cw_sites%*%svd_rotmat_sites
  choose.lv.coefs <- scaled_cw_species%*%svd_rotmat_species

  B<-(diag((bothnorms^alpha)/sqrt(colSums(gllvm::getLV(object,type = type)^2)))%*%svd_rotmat_sites)


  if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}

    if(object$row.eff=="random" && dim(object$A)[2]>dim(object$lvs)[2]){
      object$A<- object$A[,-1,-1]
    }

    if(object$randomB==FALSE|(num.lv.c+num.RR)==0|object$randomB!=FALSE&type=="residual"){
      sdb<-gllvm:::CMSEPf(object, type = type)$A

      if(num.RR>0){
        #variational covariances but add 0s for RRR
        A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
        A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
      }else{A<-object$A}
      if(type%in%c("conditional")){
        if((num.lv.c+num.lv)==1){
          A<-A*object$params$sigma.lv
        }else{
          for(i in 1:n){
            A[i,,]<-diag(object$params$sigma.lv)%*%A[i,,]%*%diag(object$params$sigma.lv)
          }
        }

      }
      if(type!="marginal"){
        object$A<-sdb+A
      }else{
        object$A<-sdb
      }

    }else{
      sdb<-CMSEPf(object, type = "residual")$A

      if(num.RR>0){
        #variational covariances but add 0s for RRR
        A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
        if((num.lv+num.lv.c)>0)A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
      }else{A<-object$A}
      if(type%in%c("conditional")){
        if((num.lv.c+num.lv)==1){
          A<-A*object$params$sigma.lv
        }else{
          for(i in 1:n){
            A[i,,]<-diag(object$params$sigma.lv)%*%A[i,,]%*%diag(object$params$sigma.lv)
          }
        }

        A<-sdb+A
      }
      #calculate covariance sfor randomB
      covsB <- CMSEPf(object, type = "residual",return.covb = T)
      covsB <- covsB[colnames(covsB)=="b_lv",colnames(covsB)=="b_lv"]

      if(object$randomB=="P"|object$randomB=="single"){
        covsB <- covsB + as.matrix(Matrix::bdiag(lapply(seq(dim(object$Ab.lv)[1]), function(k) object$Ab.lv[k , ,])))
      }else if(object$randomB=="LV"){
        covsB <- covsB + as.matrix(Matrix::bdiag(lapply(seq(dim(object$Ab.lv)[1]), function(q) object$Ab.lv[q , ,])))
      }

      if(type=="marginal"){
        #no prediction errors for the residual
        A[,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] <- array(0,dim = c(n,num.lv.c+num.RR,num.lv.c+num.RR))
      }
      S <- diag(1:(num.lv.c+num.RR))
      if(num.lv.c>0&type=="conditional"){
        diag(S)[1:num.lv.c] <- object$params$sigma.lv[1:num.lv.c]
      }
      for(i in 1:n){
        Q <- as.matrix(Matrix::bdiag(replicate(num.RR+num.lv.c,object$lv.X[i,,drop=F],simplify=F)))
        A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] <- S%*%A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)]%*%S
        temp <- Q%*%covsB%*%t(Q) #variances and single dose of covariances
        temp[col(temp)!=row(temp)] <- 2*temp[col(temp)!=row(temp)] ##should be double the covariance
        A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] <- A[i,1:(num.RR+num.lv.c),1:(num.RR+num.lv.c)] + temp
      }
      object$A <- A

    }

    r=0
    ell <- NULL
    for (i in 1:n) {
      if(!object$TMB && object$Lambda.struc == "diagonal"){
        covm <- (t(B)%*%diag(object$A[i,1:num.lv+r])%*%B)[which.lvs,which.lvs];
        # covm <- diag(object$A[i,which.lvs+r]);
      } else {
        covm <- (t(B)%*%object$A[i,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.RR+num.lv.c)+r]%*%B)[which.lvs,which.lvs];
        # covm <- object$A[i,which.lvs+r,which.lvs+r];
      }
      ell<-rbind(ell,cbind(ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=num.lv+num.RR+num.lv.c)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips),i))


    }
    ell <- data.frame(ell)
    colnames(ell)<-c("x","y","group")
    return(ell)
    }



  ellipse<-function(center, covM, rad, col=4, lwd = 1, lty = 1){
    seg <- 51
    Qc <- chol(covM, pivot = TRUE)
    angles <- (0:seg) * 2 * pi / seg
    unit.circ <- cbind(cos(angles), sin(angles))
    order <- order(attr(Qc, "pivot"))
    ellips <- t(center + rad * t(unit.circ %*% Qc[, order]))
    return(ellips)
  }
