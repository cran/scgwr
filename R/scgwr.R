scgwr      <- function (coords, y, x=NULL, knn = 100, kernel = "gau", p = 4, approach="CV", nsamp=NULL ){

  ic_score <- function(par, X, y) {
    par2 <- par^2
    R0   <- rep(par2[1], p + 1)
    for (r0 in 2:(p + 1)) R0[r0] <- par2[1]^r0
    R0   <-R0/sum( R0 )

    MMM  <- R0[id_p] * MM
    mmm  <- R0[id_py]* mm
    MMM2 <- 2*par2[2]*( R0[id_p] * MM ) + ( R0[id_p]^2 * MM2 )

    MMMsum <- NULL
    MMMsum2<- NULL
    mmmsum <- NULL
    for (ll1 in 1:nx) {
      for (ll2 in 1:nx) {
        MMM_sub <- colSums(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
        MMMsum <- rbind(MMMsum, MMM_sub)

        MMM_sub2<- colSums(MMM2[(id_x1 == ll1) & (id_x2 == ll2), ])
        MMMsum2 <- rbind(MMMsum2, MMM_sub2)
      }
      mmm_sub <- colSums(mmm[id_xy == ll1, ])
      mmmsum <- rbind(mmmsum, mmm_sub)
    }
    MMMlist <- list(NULL)
    MMMlist2<- list(NULL)
    mmmlist <- list(NULL)
    xxxlist <- list(NULL)# code is inefficient (it must not needed)
    for (ll3 in 1:nsamp) {
      MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx) + par2[2] * XX0
      MMMlist2[[ll3]]<- matrix(MMMsum2[, ll3], nx, nx) + par2[2] ^ 2 * XX0
      mmmlist[[ll3]] <- mmmsum[, ll3] + par2[2] * Xy0
      xxxlist[[ll3]] <- X[ll3,]
    }
    tryres <- try(MMMinv <- lapply(MMMlist, solve), silent = TRUE)
    if (class(tryres) == "try-error") {
      obj <- 10^10
    } else {
      beta <- t(matrix(unlist(mapply("%*%", MMMinv, mmmlist,
                                     SIMPLIFY = FALSE)), nrow = nx, ncol = nsamp))

      trS00  <- mapply("%*%", MMMinv , xxxlist,SIMPLIFY = FALSE)
      trS0   <- mapply("%*%", xxxlist, trS00   ,SIMPLIFY = FALSE)
      trS    <- Reduce("+", trS0)#*g00

      trSS00  <- mapply("%*%", MMMlist2 , trS00,   SIMPLIFY = FALSE)
      trSS0   <- mapply(function(x,y) sum(x*y), trS00, trSS00, SIMPLIFY = FALSE)
      trSS    <- Reduce("+", trSS0)

      pred    <- rowSums(X * beta)#[nlist,]
      sse     <- sum((y - pred)^2)#[nlist]
      sig     <- sqrt( sse/nsamp )
      AICc    <- 2*nsamp*log(sig)+nsamp*log(2*pi)+nsamp*(nsamp+trS)/(nsamp-2-trS)
    }

    return(AICc)
  }

  cv_score <- function(par,X,y) {
    par2 <- par^2
    R0   <- rep(par2[1], p + 1)
    for (r0 in 2:(p + 1)) R0[r0] <- par2[1]^r0
    R0   <-R0/sum( R0 )

    MMM  <- R0[id_p] * MM#[,nlist]
    mmm  <- R0[id_py] * mm#[,nlist]
    MMMsum <- NULL
    mmmsum <- NULL
    for (ll1 in 1:nx) {
      for (ll2 in 1:nx) {
        MMM_sub<- colSums(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
        MMMsum <- rbind(MMMsum, MMM_sub)
      }
      mmm_sub <- colSums(mmm[id_xy == ll1, ])
      mmmsum <- rbind(mmmsum, mmm_sub)
    }
    MMMlist <- list(NULL)
    mmmlist <- list(NULL)

    for (ll3 in 1:nsamp) {
      MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx) + par2[2] * XX0
      mmmlist[[ll3]] <- mmmsum[, ll3] + par2[2] * Xy0
    }
    tryres <- try(MMMinv <- lapply(MMMlist, solve), silent = TRUE)
    if (class(tryres) == "try-error") {
      error <- 10^10
    }
    else {
      beta <- t(matrix(unlist(mapply("%*%", MMMinv, mmmlist,
                                     SIMPLIFY = FALSE)), nrow = nx, ncol = nsamp))
      pred <- rowSums(as.matrix(X) * beta)#[nlist,]
      error <- sum((y - pred)^2)#[nlist]
    }
    return(error)
  }

  scgwr_summary <- function(par) {

    par2 <- par^2
    #n   <- length(y)

    if( 100001 >= n){#nsamp == n
      niter<-1
      ilist<-list( 1:n )
    } else {
      niter<-round(n/100000)
      ilist<-list( 1:100000 )
      nmax <-0
      ii   <-1
      while(nmax<=n){
        ilist0     <- c(nmax + 1:100000)
        ilist0_id  <- ilist0 <= n
        if( sum( ilist0_id ) >= 1 ){
          ilist[[ii]]<- ilist0[ilist0_id]
        }
        nmax<-nmax+100000
        ii  <-ii+1
      }
    }

    BETA       <-matrix(0,nrow=n,ncol=nx)
    BSE_no_sig <-BETA
    ENP        <-0
    PRED       <-rep(0,n)
    for(iii in 1:length(ilist)){
      jd  <- ilist[[iii]]
      n_iii<-length(jd)
      Dnn <- get.knnx(coords,coords[jd,], knn)#not knn -1 because of using get.knnx
      Did  <- Dnn$nn.index
      if (kernel == "gau") {
        G0    <- exp(-(Dnn$nn.dist/ban0)^2)
      } else if (kernel == "exp") {
        G0    <- exp(-Dnn$nn.dist/ban0)
      }

      MM <- matrix(0, ncol = n_iii, nrow = length(id_p))
      MM2<- matrix(0, ncol = n_iii, nrow = length(id_p))
      mm <- matrix(0, ncol = n_iii, nrow = sum(id_x1 == 1))
      for (i in 1:n_iii) {
        G <- matrix(1, nrow = knn, ncol = p + 1)
        for (p0 in 1:p) G[, p0 + 1] <- G0[i, ]^(2^(p/2)/2^p0)
        X_sub <- as.matrix(X[Did[i, ],])
        y_sub <- y[Did[i, ]]
        kk0   <- 1
        for (k0 in 1:nx) {
          GXp <- G   * X_sub[, k0]
          GXp2<- G^2 * X_sub[, k0]
          for (k1 in 1:(p + 1)) {
            MM [(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp [,k1] * X_sub)
            MM2[(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp2[,k1] * X_sub)
            mm [(id_xy == k0) & (id_py == k1), i] <- sum(GXp[,k1] * y_sub)
          }
        }
      }

      R0   <- rep(par2[1], p + 1)
      for (r0 in 2:(p + 1)) R0[r0] <- par2[1]^r0
      R0   <-R0/sum( R0 )

      MMM  <- R0[id_p] * MM
      mmm  <- R0[id_py] * mm
      MMM2 <- 2*par2[2]*MMM + R0[id_p]^2 * MM2

      MMMsum <- NULL
      MMMsum2<- NULL
      mmmsum <- NULL
      for (ll1 in 1:nx) {
        for (ll2 in 1:nx) {
          MMM_sub <- colSums(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
          MMMsum  <- rbind(MMMsum, MMM_sub)

          MMM_sub2<- colSums(MMM2[(id_x1 == ll1) & (id_x2 == ll2), ])
          MMMsum2 <- rbind(MMMsum2, MMM_sub2)
        }
        mmm_sub <- colSums(mmm[id_xy == ll1, ])
        mmmsum <- rbind(mmmsum, mmm_sub)
      }
      MMMlist <- list(NULL)
      MMMlist2<- list(NULL)
      mmmlist <- list(NULL)
      xxxlist <- list(NULL)# code is inefficient (it must not needed)
      for (ll3 in 1:n_iii) {
        MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx)  + par2[2] * XX0
        MMMlist2[[ll3]]<- matrix(MMMsum2[, ll3], nx, nx) + par2[2] ^ 2 * XX0
        mmmlist[[ll3]] <- mmmsum[, ll3] + par2[2] * Xy0
        xxxlist[[ll3]] <- as.matrix( X[ilist[[iii]],] )[ll3,]
      }
      tryres <- try(MMMinv <- lapply(MMMlist, solve), silent = TRUE)
      beta   <- t(matrix(unlist(mapply("%*%", MMMinv, mmmlist,
                                       SIMPLIFY = FALSE)), nrow = nx, ncol = n_iii))

      trS00  <- mapply("%*%", MMMinv , xxxlist,SIMPLIFY = FALSE)
      trS0   <- mapply("%*%", xxxlist, trS00   ,SIMPLIFY = FALSE)
      trS    <- Reduce("+", trS0)#*g00

      trSS00  <- mapply("%*%", MMMlist2 , trS00,   SIMPLIFY = FALSE)
      trSS0   <- mapply(function(x,y) sum(x*y), trS00, trSS00, SIMPLIFY = FALSE)
      trSS    <- Reduce("+", trSS0)
      enp     <- 2*trS - trSS

      pred    <- rowSums(X[ilist[[iii]],] * beta)
      bse00   <- mapply("%*%", MMMlist2, MMMinv,SIMPLIFY = FALSE)
      bse0    <- mapply("%*%", MMMinv  , bse00 ,SIMPLIFY = FALSE)
      bse0    <- lapply(bse0, function(x) sqrt(diag( x )))
      bse     <- matrix(unlist(bse0), nrow = n_iii, ncol = nx, byrow=TRUE)

      BETA[ilist[[iii]],]       <- beta
      BSE_no_sig[ilist[[iii]],] <- bse
      ENP     <-ENP + enp
      PRED[ilist[[iii]]]        <- pred
    }

    pred    <- PRED
    sse     <- sum((y - pred)^2)
    sig     <- sqrt( sse/( n - enp ) )
    bse     <- c(sig)*BSE_no_sig#sqrt(sig) is replaced with sig
    R2      <- cor(pred, y)^2
    adjR2   <- 1 - (1 - R2) * (n - 1)/(n - enp - 1)

    beta    <- BETA
    bt      <- beta/bse
    bp0     <- data.frame(2 - 2 * pt(abs(bt), df = c(n) - enp))
    bp      <- bp0    #bp<- data.frame(bp0 * c(abs(enp/nx))) # this part is changed
    bp[bp > 1] <- 1

    beta    <- data.frame(beta)
    bse     <- data.frame(bse)
    bt      <- data.frame(beta/bse)
    names(beta) <- names(bse) <- names(bt) <- names(bp) <-names(bp0) <- xname

    r       <- par2[1]
    pen     <- par2[2]
    spar    <- data.frame(c(r, pen))
    names(spar)    <- "Estimate"
    rownames(spar) <- c("scale", "penalty")

    sig0    <- sqrt( sse/n )
    loglik  <-  -n*log(sig0)- n/2*log(2*pi)
    AICc    <- -2*loglik + n*(n+trS)/(n-2-trS)

    if(approach=="CV"){
      e_stat <- data.frame(stat = c( sig, R2, adjR2, loglik,AICc,cvscore))#sqrt(sig) is replaced with sig
      rownames(e_stat) <- c("resid_SE", "R2", "adjR2", "logLik","AICc","CV_score(RMSE)")
    } else {
      e_stat <- data.frame(stat = c( sig, R2, adjR2, loglik,AICc))   #sqrt(sig) is replaced with sig
      rownames(e_stat) <- c("resid_SE", "R2", "adjR2", "logLik","AICc")
    }

    return(list(beta=beta, bse=bse, bt=bt, bp=bp, bp0=bp0, e_stat=e_stat, spar=spar, pred=pred,
                resid=resid, enp = enp ))
  }

  proc  <-function(coords, X, y, coords_samp, X_samp, y_samp, knn, approach, kernel,p){
    n_samp <- length(y_samp)
    Dnn    <- get.knnx(coords,coords_samp, knn-1)
    if( approach != "CV"){
      Dnn$nn.index<-as.matrix( cbind( 1:n_samp, Dnn$nn.index) )
      Dnn$nn.dist <-as.matrix( cbind( 0  , Dnn$nn.dist ) )
    }

    bmax <- sqrt((max(coords[,1])-min(coords[,1]))^2+(max(coords[,2])-min(coords[,2]))^2)
    Dband<- seq(0,bmax,len=10000)
    Did  <- Dnn$nn.index
    if (kernel == "gau") {
      ban0  <- median(Dnn$nn.dist[, min(50,knn)])/sqrt(3)
      G0    <- exp(-(Dnn$nn.dist/ban0)^2)
      G0band<- exp(-(Dband/ban0)^2)
    } else if (kernel == "exp") {
      ban0  <- median(Dnn$nn.dist[, min(50,knn)])/3
      G0    <- exp(-Dnn$nn.dist/ban0)
      G0band<- exp(-Dband/ban0)
    }

    XX0   <- crossprod(X)
    Xy0   <- crossprod(X,y)
    id_p  <- rep(rep(1:(p + 1), each = nx), nx)
    id_x1 <- rep(1:nx, each = nx * (p + 1))
    id_x2 <- rep(1:nx, nx * (p + 1))
    id_py <- id_p [id_x1 == 1]
    id_xy <- id_x2[id_x1 == 1]
    MM    <- matrix(0, ncol = n_samp, nrow = length(id_p))
    mm    <- matrix(0, ncol = n_samp, nrow = sum(id_x1 == 1))
    if(approach=="AICc"){
      MM2 <- matrix(0, ncol = n_samp, nrow = length(id_p))
      G   <- matrix(1, nrow = knn, ncol = p + 1)
    } else {
      MM2 <- NULL
      G   <- matrix(1, nrow = knn - 1, ncol = p + 1)
    }
    for (i in 1:n_samp){
      for (p0 in 1:p) G[, p0 + 1] <- G0[i, ]^(2^(p/2)/2^p0)
      X_sub <- as.matrix(X[Did[i, ],])
      y_sub <- y[Did[i, ]]
      kk0   <- 1
      for (k0 in 1:nx) {
        GXp   <- G * X_sub[, k0]
        if(approach=="AICc"){
          GXp2<- G^2 * X_sub[, k0]
        }
        for (k1 in 1:(p + 1)) {
          MM [(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp [,k1] * X_sub)
          mm [(id_xy == k0) & (id_py == k1), i] <- sum(GXp[,k1] * y_sub)
          if(approach=="AICc"){
            MM2[(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp2[,k1] * X_sub)
          }
        }
      }
    }
    return(list(MM=MM, MM2=MM2, mm=mm, XX0=XX0,Xy0=Xy0,id_p=id_p,id_x1=id_x1,id_x2=id_x2,
                id_py=id_py,id_xy=id_xy,G0=G0,Dnn=Dnn, ban0=ban0))
  }

  ##########################################################
  n <- length(y)
  if (knn >= n) knn   <- n - 1
  if (is.null(x)) {
    X     <- as.matrix(rep(1, n))
    xname <- "(Intercept)"
    x_id <- NULL
  } else {
    X00 <- as.matrix(x)
    if (is.numeric(X00) == F) {
      mode(X00) <- "numeric"
    }

    if( dim( as.matrix(x) )[2]==1 ){
      x_id <- (sd(x) !=0)
    } else {
      x_id <- (apply(x,2,sd) !=0)#(x %>% summarise_if(is.numeric, funs(sd)) !=0)
    }
    x_nm <- colnames(x)[x_id] #added

    if( all(x_id == 0)){
      X <- as.matrix(rep(1, n))
      xname <- "(Intercept)"
      x_id <- NULL
    } else {
      X <- as.matrix(cbind(1, X00[, x_id]))
      xname <- c("(Intercept)", names(as.data.frame( X00[, x_id])))
    }
  }
  nx  <- dim(X)[2]

  if(approach=="AICc"){
    if( !is.null( nsamp ) ){
      message("nsamp is ignored because it is currently not available for the AICc-based estimation")
    }

    nlist<-1:n
    nsamp<-n
  } else if( !is.null(nsamp) ){
    if(nsamp < n){
      nlist<-sample(n,nsamp)
    } else {
      nlist<-1:n
      nsamp<-n
    }
  } else {
    nlist<-1:n
    nsamp<-n
  }

  coords_samp<- coords[nlist,]
  X_samp     <- X[nlist,]
  y_samp     <- y[nlist]

  misc       <- proc(coords=coords,coords_samp=coords_samp, X=X, y=y, X_samp=X_samp, y_samp=y_samp, knn=knn, approach=approach, kernel=kernel, p=p)
  MM         <- misc$MM
  MM2        <- misc$MM2
  mm         <- misc$mm
  XX0        <- misc$XX0
  Xy0        <- misc$Xy0
  G0         <- misc$G0
  ban0       <- misc$ban0
  id_p       <- misc$id_p
  id_x1      <- misc$id_x1
  id_x2      <- misc$id_x2
  id_py      <- misc$id_py
  id_xy      <- misc$id_xy

  par0       <- c(1, 0.01)
  if(approach == "CV"){
    res0    <- optim(par0,fn = cv_score, X=X_samp,y=y_samp,method ="L-BFGS-B",lower=c(-10,-Inf),upper=c(10,Inf))
    cvscore <- sqrt(res0$value/n)
  } else {
    res0    <- optim(par0,fn = ic_score, X=X_samp,y=y_samp, method ="L-BFGS-B",lower=c(-10,-Inf),upper=c(10,Inf))
    cvscore <-NULL
  }
  par   <- res0$par
  summ  <- scgwr_summary(par)
  beta  <- summ$beta
  bse   <- summ$bse
  bt    <- summ$bt
  bp    <- summ$bp
  bp0   <- summ$bp0
  e_stat<- summ$e_stat
  pred  <- summ$pred
  resid <- summ$resid
  spar  <- summ$spar
  enp   <- summ$enp

  other <- list( kernel = kernel, knn = knn, x_id = x_id, xname = xname, ban0 = ban0, nx=nx,
                 p=p, X=X, y=y, par = par, XX0 = XX0, Xy0 = Xy0, enp = enp, coords=coords )

  result<- list(b = beta, bse = bse, t = bt, p = bp0, par = spar,
                e = e_stat, pred = pred, resid = resid, other = other, call = match.call() )
  class( result ) <- "scgwr"
  return( result )
}

print.scgwr <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Spatially varying coefficients on x (summary)----\n")
  cat("\nCoefficient estimates:\n")
  print( summary( x$b ) )
  cat("\nStatistical significance:\n")
  p01<-apply(x$p,2,function(x) sum(x<0.01))
  p05<-apply(x$p,2,function(x) sum(x<0.05)) - p01
  p10<-apply(x$p,2,function(x) sum(x<0.10)) - p01 - p05
  p90<-length(x$p[,1]) - p01 - p05 - p10
  pv <-data.frame(rbind( p90, p10, p05, p01))
  names(pv)[1]  <- "Intercept"
  row.names(pv) <- c("Not significant", "Significant (10% level)",
                     "Significant ( 5% level)","Significant ( 1% level)")
  print(pv)
  cat("\n----Variance parameters----------------------------------\n")
  print( x$par )
  cat("\n----Error statistics-------------------------------------\n")
  print(x$e)
  invisible(x)
}
