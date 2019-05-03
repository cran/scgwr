scgwr <- function (coords, y, x, knn = 100, kernel = "gau", p = 4){

  cv_score2 <- function(par) {
    par2 <- par^2
    R0   <- rep(par2[1], p + 1)
    for (r0 in 2:(p + 1)) R0[r0] <- par2[1]^r0
    MMM  <- R0[id_p] * MM
    mmm  <- R0[id_py] * mm
    MMMsum <- NULL
    mmmsum <- NULL
    for (ll1 in 1:nx) {
      for (ll2 in 1:nx) {
        MMM_sub <- colSums(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
        MMMsum <- rbind(MMMsum, MMM_sub)
      }
      mmm_sub <- colSums(mmm[id_xy == ll1, ])
      mmmsum <- rbind(mmmsum, mmm_sub)
    }
    MMMlist <- list(NULL)
    mmmlist <- list(NULL)
    for (ll3 in 1:n) {
      MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx) + par2[2] * XX0
      mmmlist[[ll3]] <- mmmsum[, ll3] + par2[2] * Xy0
    }
    tryres <- try(MMMinv <- lapply(MMMlist, solve), silent = TRUE)
    if (class(tryres) == "try-error") {
      error <- 10^10
    }
    else {
      beta <- t(matrix(unlist(mapply("%*%", MMMinv, mmmlist,
                                     SIMPLIFY = FALSE)), nrow = nx, ncol = n))
      pred <- rowSums(X * beta)
      error <- sum((y - pred)^2)
    }
    return(error)
  }

  n <- length(y)
  if (knn >= n)
    knn <- n - 1
  if (is.null(x)) {
    X <- as.matrix(rep(1, n))
    xname <- "(Intercept)"
    x_id <- NULL
  } else {
    X00 <- as.matrix(x)
    if (is.numeric(X00) == F) {
      mode(X00) <- "numeric"
    }
    x_id <- (x %>% summarise_if(is.numeric, funs(sd)) !=0)
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
  nx <- dim(X)[2]
  Dnn <- get.knn(coords, knn)
  Did <- Dnn$nn.index
  if (kernel == "gau") {
     ban0 <- median(Dnn$nn.dist[, min(50,knn)])/sqrt(3)
     G0 <- exp(-(Dnn$nn.dist/ban0)^2)
  } else if (kernel == "exp") {
     ban0 <- median(Dnn$nn.dist[, min(50,knn)])/3
     G0 <- exp(-Dnn$nn.dist/ban0)
  }

  XX0   <- crossprod(X)
  Xy0   <- crossprod(X,y)
  id_p  <- rep(rep(1:(p + 1), each = nx), nx)
  id_x1 <- rep(1:nx, each = nx * (p + 1))
  id_x2 <- rep(1:nx, nx * (p + 1))
  id_py <- id_p[id_x1 == 1]
  id_xy <- id_x2[id_x1 == 1]
  MM <- matrix(0, ncol = n, nrow = length(id_p))
  mm <- matrix(0, ncol = n, nrow = sum(id_x1 == 1))
  for (i in 1:n) {
    G <- matrix(1, nrow = knn + 1, ncol = p + 1)
    for (p0 in 1:p) G[, p0 + 1] <- c(1, G0[i, ])^(2^(p/2)/2^p0)
    X_sub <- X[Did[i, ],]
    y_sub <- y[Did[i, ]]
    kk0   <- 1
    for (k0 in 1:nx) {
      GXp <- G[-1, ] * X_sub[, k0]
      for (k1 in 1:(p + 1)) {
        MM[(id_x1 == k0) & (id_p == k1), i]  <- colSums(GXp[,k1] * X_sub)
        mm[(id_xy == k0) & (id_py == k1), i] <- sum(GXp[,k1] * y_sub)
      }
    }
  }
  par <- c(1, 0.01)
  res0 <- optim(fn = cv_score2, par)
  res0$value
  res0$par
  r <- res0$par[1]^2
  pen <- res0$par[2]^2
  cvscore <-res0$value
  R0 <- rep(r, p + 1)
  for (r0 in 2:(p + 1)) R0[r0] <- r^r0
  Xp <- matrix(0, ncol = nx * (p + 1), nrow = knn + 1)
  RGXsum <- matrix(0, nrow = n, ncol = (p + 1) * nx)
  gX <- list(NULL)
  for (i in 1:n) {
    G <- matrix(R0[1], nrow = knn + 1, ncol = p + 1)
    for (p0 in 1:p) G[, p0 + 1] <- R0[p0] * c(1, G0[i, ])^(2^(p/2)/2^p0)
    X_sub <- X[c(i, Did[i, ]), ]
    gX[[i]] <- rowSums(G) * X_sub
    kk0 <- 1
    for (k0 in 1:nx) {
      Xp[, kk0:(kk0 + p)] <- G * X_sub[, k0]
      kk0 <- kk0 + p + 1
    }
    RGXsum[i, ] <- colSums(Xp)
  }
  for (i in 1:n) {
    kk0 <- 1
    for (k0 in 1:nx) {
      for (k1 in 1:(p + 1)) {
        MM[(id_x1 == k0) & (id_p == k1), i] <- MM[(id_x1 == k0) & (id_p == k1), i] + X[i, k0] * X[i, ]
        mm[(id_xy == k0) & (id_py == k1), i] <- mm[(id_xy ==k0) & (id_py == k1), i] + X[i, k0] * y[i]
      }
      kk0 <- kk0 + p + 1
    }
  }
  MMM <- R0[id_p] * MM
  mmm <- R0[id_py] * mm
  MMMsum <- NULL
  mmmsum <- NULL
  for (ll1 in 1:nx) {
    for (ll2 in 1:nx) {
      MMM_sub <- colSums(MMM[(id_x1 == ll1) & (id_x2 == ll2), ])
      MMMsum <- rbind(MMMsum, MMM_sub)
    }
    mmm_sub <- colSums(mmm[id_xy == ll1, ])
    mmmsum <- rbind(mmmsum, mmm_sub)
  }
  MMMlist <- list(NULL)
  MMMlist2 <- list(NULL)
  mmmlist <- list(NULL)
  for (ll3 in 1:n) {
    MMMlist[[ll3]] <- matrix(MMMsum[, ll3], nx, nx) + pen * XX0
    MMMlist2[[ll3]] <- matrix(MMMsum[, ll3], nx, nx)
    mmmlist[[ll3]] <- mmmsum[, ll3] + pen * Xy0
  }
  MMMinv <- lapply(MMMlist, solve)
  beta0 <- mapply("%*%", MMMinv, mmmlist, SIMPLIFY = FALSE)
  beta <- matrix(NA, nrow = n, ncol = nx)
  for (l8 in 1:n) {
    beta[l8, ] <- beta0[[l8]]
  }
  pred <- rowSums(X * beta)
  resid <- y - pred
  SSE <- c(sum(resid^2))
  x <- as.list(as.data.frame(t(X)))
  trHH0 <- as.list(as.data.frame(matrix(mapply(function(x,y) {
    return(x %*% c(y))
  }, MMMinv, x), nrow = nx)))
  trH <- 0
  for (i in 1:n) {
    G <- matrix(R0[1], nrow = knn + 1, ncol = p + 1)
    for (p0 in 1:p) G[, p0 + 1] <- R0[p0] * c(1, G0[i, ])^(2^(p/2)/2^p0)
    g <- rowSums(G)
    X_sub <- X[c(i, Did[i, ]), ]
    trH <- trH + (trHH0[[i]] %*% c(X_sub[1, ] * g[1]))
  }

  G2 <- matrix(0, nrow = n, ncol = knn + 1)
  for (p0 in 1:p) G2 <- G2 + as.matrix(cbind(1, G0))^(2^(p/2)/2^p0) * R0[p]
  G22 <- G2
  XG2X <- list(NULL)
  for (i in 1:n) {
    XG2X_0 <- matrix(0, nx, nx)
    for (k01 in 1:nx) {
      for (k02 in 1:nx) {
        XG2X_0[k01, k02] <- sum(X[c(i, Did[i, ]), k01] *
                                  (G2[i, ]^2 + 2 * pen * G2[i, ]) * X[c(i, Did[i, ]), k02]) + pen^2 * XX0[k01, k02]
      }
    }
    XG2X[[i]] <- XG2X_0
  }
  bse_0 <- matrix(0, ncol = nx, nrow = n)
  trHH_0<- rep(0, n)
  for (l9 in 1:n) {
    XXMat      <-MMMinv[[l9]] %*% XG2X[[l9]] %*% MMMinv[[l9]]
    bse_0 [l9,]<- sqrt(diag(XXMat))
    trHH_0[l9 ]<- c(X[l9,]) %*% XXMat %*% c(X[l9,])
  }

  trHH<-sum( trHH_0 )
  enp <- 2 * trH - trHH
  sig <- c(sqrt(SSE/(n - enp)))
  bse <- sig*bse_0

  bt <- beta/bse
  bp0 <- data.frame(2 - 2 * pt(abs(bt), df = c(n) - enp))
  bp <- data.frame(bp0 * c(abs(enp/nx)))
  bp[bp > 1] <- 1
  beta <- data.frame(beta)
  bse <- data.frame(bse)
  bt <- data.frame(beta/bse)
  names(beta) <- names(bse) <- names(bt) <- names(bp0) <- names(bp) <- xname
  spar <- data.frame(c(r, pen))
  names(spar) <- "Estimate"
  rownames(spar) <- c("scale", "penalty")
  loglik <- -n * log(sig) - n/2 * log(2 * pi)
  AICc <- -2 * loglik + n * (n + trH)/(n - 2 - trH)
  R2 <- cor(pred, y)^2
  adjR2 <- 1 - (1 - R2) * (n - 1)/(n - enp - 1)
  cvscore2<- sqrt(cvscore/n)
  e_stat <- data.frame(stat = c(SSE, sig, R2, adjR2, loglik,
                                AICc,cvscore2))
  rownames(e_stat) <- c("SSE", "resid_SE", "R2", "adjR2", "logLik",
                        "AICc","CV_score(RMSE)")
  return(list(b = beta, bse = bse, t = bt, p = bp0, pa = bp, par = spar,
              e = e_stat, pred = pred, resid = resid))
}
