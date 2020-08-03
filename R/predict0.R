predict0<-function( mod, coords0, x0 = NULL ){

  coords<- mod$other$coords
  kernel<- mod$other$kernel
  knn   <- mod$other$knn
  x_id  <- mod$other$x_id
  nx    <- mod$other$nx
  ban0  <- mod$other$ban0
  X     <- mod$other$X
  y     <- mod$other$y
  XX0   <- mod$other$XX0
  Xy0   <- mod$other$Xy0
  par   <- mod$other$par
  sig   <- mod$e$stat[1]
  enp   <- mod$other$enp
  xname <- mod$other$xname
  p     <- mod$other$p
  n     <- length( y )

  n0    <- dim(coords0)[1]
  if( !is.null(x0 ) ){
    X00 <- as.matrix(x0)
    if (is.numeric(X00) == F) {
      mode(X00) <- "numeric"
    }

    if( all(x_id == 0)){
      X0 <- as.matrix(rep(1, n))
    } else {
      if( n0 ==1 ){
        X0 <- as.matrix(t(c(1, c(X00[, x_id]))))
      } else {
        X0 <- as.matrix(cbind(1, X00[, x_id]))
      }
    }
  }


  Dnn0<- get.knnx(coords, coords0, knn)

  Did0<- Dnn0$nn.index
  if (kernel == "gau") {
    G0    <- exp(-(Dnn0$nn.dist/ban0)^2)
    #G0band<- exp(-(Dband/ban0)^2)
  } else if (kernel == "exp") {
    G0    <- exp(-Dnn0$nn.dist/ban0)
    #G0band<- exp(-Dband/ban0)
  }

  id_p  <- rep(rep(1:(p + 1), each = nx), nx)
  id_x1 <- rep(1:nx, each = nx * (p + 1))
  id_x2 <- rep(1:nx, nx * (p + 1))
  id_py <- id_p [id_x1 == 1]
  id_xy <- id_x2[id_x1 == 1]
  MM    <- matrix(0, ncol = n0, nrow = length(id_p))
  MM2   <- matrix(0, ncol = n0, nrow = length(id_p))
  mm    <- matrix(0, ncol = n0, nrow = sum(id_x1 == 1))
  G     <- matrix(1, nrow = knn, ncol = p + 1)

  for (i in 1:n0){
    for (p0 in 1:p) G[, p0 + 1] <- G0[i, ]^(2^(p/2)/2^p0)
    X_sub <- as.matrix( X[Did0[i, ],] )
    y_sub <- y[Did0[i, ]]
    kk0   <- 1
    for (k0 in 1:nx) {
      GXp   <- G * X_sub[, k0]
      GXp2  <- G^2 * X_sub[, k0]
      for (k1 in 1:(p + 1)) {
        MM [(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp [,k1] * X_sub)
        MM2[(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp2[,k1] * X_sub)
        mm [(id_xy == k0) & (id_py == k1), i] <- sum(GXp[,k1] * y_sub)
      }
    }
  }

  par2 <- par^2
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
      #if( n0==1 ){
      #  MMM_sub <- sum(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
      #  MMM_sub2<- sum(MMM2[(id_x1 == ll1) & (id_x2 == ll2), ])
      #} else {
        MMM_sub <- colSums(as.matrix(MMM[(id_x1 == ll1) & (id_x2 == ll2), ]))
        MMM_sub2<- colSums(as.matrix(MMM2[(id_x1 == ll1) & (id_x2 == ll2), ]))
      #}

      MMMsum  <- rbind(MMMsum,  MMM_sub)
      MMMsum2 <- rbind(MMMsum2, MMM_sub2)
    }
    mmm_sub <- colSums(as.matrix(mmm[id_xy == ll1, ]))
    mmmsum <- rbind(mmmsum, mmm_sub)
  }
  MMMlist <- list(NULL)
  MMMlist2<- list(NULL)
  mmmlist <- list(NULL)
  #xxxlist <- list(NULL)# code is inefficient (it must not needed)
  for (ll3 in 1:n0) {
    MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx)  + par2[2] * XX0
    MMMlist2[[ll3]]<- matrix(MMMsum2[, ll3], nx, nx) + par2[2] ^ 2 * XX0
    mmmlist[[ll3]] <- mmmsum[, ll3] + par2[2] * Xy0
    #xxxlist[[ll3]] <- X[ll3,]
  }
  tryres <- try(MMMinv <- lapply(MMMlist, solve), silent = TRUE)
  beta   <- t(matrix(unlist(mapply("%*%", MMMinv, mmmlist,
                                   SIMPLIFY = FALSE)), nrow = nx, ncol = n0))

  if( !is.null(x0)){
    pred    <- rowSums(X0 * beta)
  } else {
    pred    <- beta
  }

  bse00   <- mapply("%*%", MMMlist2, MMMinv,SIMPLIFY = FALSE)
  bse0    <- mapply("%*%", MMMinv  , bse00 ,SIMPLIFY = FALSE)
  bse0    <- lapply(bse0, function(x) c( sig )*sqrt(diag( x )))
  bse     <- matrix(unlist(bse0), nrow = n0, ncol = nx, byrow=TRUE)

  bt      <- beta/bse
  bp0     <- data.frame(2 - 2 * pt(abs(bt), df = c(n) - enp))
  bp      <- bp0
  bp[bp > 1] <- 1

  beta <- data.frame(beta)
  bse <- data.frame(bse)
  bt <- data.frame(beta/bse)
  names(beta) <- names(bse) <- names(bt) <- names(bp) <-names(bp0) <- xname

  return(list( pred=pred, b=beta, bse=bse, t=bt, p=bp ))
}
