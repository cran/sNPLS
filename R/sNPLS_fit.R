#' Fit a sNPLS model
#'
#' @description Fits a N-PLS regression model imposing sparsity on \code{wj} and \code{wk} matrices
#' @param XN A three-way array containing the predictors.
#' @param Y A matrix containing the response.
#' @param ncomp Number of components in the projection
#' @param threshold_j Threshold value on Wj. Scaled between [0, 1)
#' @param threshold_k Threshold value on Wk. scaled between [0, 1)
#' @param keepJ Number of variables to keep for each component, ignored if threshold_j is provided
#' @param keepK Number of 'times' to keep for each component, ignored if threshold_k is provided
#' @param scale.X Perform unit variance scaling on X?
#' @param center.X Perform mean centering on X?
#' @param scale.Y Perform unit variance scaling on Y?
#' @param center.Y Perform mean centering on Y?
#' @param conver Convergence criterion
#' @param max.iteration Maximum number of iterations
#' @param silent Show output?
#' @param method Select between L1 penalization (sNPLS), variable selection with Selectivity Ratio (sNPLS-SR) or variable selection with VIP (sNPLS-VIP)
#' @return A fitted sNPLS model
#' @references C. A. Andersson and R. Bro. The N-way Toolbox for MATLAB Chemometrics & Intelligent Laboratory Systems. 52 (1):1-4, 2000.
#' @references Hervas, D. Prats-Montalban, J. M., Garcia-Cañaveras, J. C., Lahoz, A., & Ferrer, A. (2019). Sparse N-way partial least squares by L1-penalization. Chemometrics and Intelligent Laboratory Systems, 185, 85-91.
#' @examples
#' X_npls<-array(rpois(7500, 10), dim=c(50, 50, 3))
#'
#' Y_npls <- matrix(2+0.4*X_npls[,5,1]+0.7*X_npls[,10,1]-0.9*X_npls[,15,1]+
#' 0.6*X_npls[,20,1]- 0.5*X_npls[,25,1]+rnorm(50), ncol=1)
#' #Discrete thresholding
#' fit <- sNPLS(X_npls, Y_npls, ncomp=3, keepJ = rep(2,3) , keepK = rep(1,3))
#' #Continuous thresholding
#' fit2 <- sNPLS(X_npls, Y_npls, ncomp=3, threshold_j=0.5, threshold_k=0.5)
#' #USe sNPLS-SR method
#' fit3 <- sNPLS(X_npls, Y_npls, ncomp=3, threshold_j=0.5, threshold_k=0.5, method="sNPLS-SR")
#' @importFrom stats sd
#' @export
sNPLS <- function(XN, Y, ncomp = 2, threshold_j=0.5, threshold_k=0.5, keepJ = NULL, keepK = NULL, scale.X=TRUE, center.X=TRUE,
         scale.Y=TRUE, center.Y=TRUE, conver = 1e-16, max.iteration = 10000, silent = F, method="sNPLS"){

  mynorm <- function(x) sqrt(sum(diag(crossprod(x))))
  thresholding <- function(x, nj) {
    ifelse(abs(x) > abs(x[order(abs(x))][nj]),
           (abs(x) - abs(x[order(abs(x))][nj])) *
             sign(x), 0)
  }
  rel_thresholding <- function(x, j_rel){
    ifelse(abs(x)-max(abs(x))*j_rel <= 0, 0, sign(x)*(abs(x)-max(abs(x))*j_rel))
  }
  if(!method %in% c("sNPLS", "sNPLS-SR", "sNPLS-VIP")) stop("'method' not recognized")
  if (length(dim(Y)) == 3) Y <- unfold3w(Y)
  if (length(dim(XN)) != 3) stop("'XN' is not a three-way array")
  if (!is.null(rownames(XN))){
    y.names <- x.names <- rownames(XN)
  } else {
    y.names <- x.names <- 1:dim(XN)[1]
  }
  if (!is.null(colnames(XN))){
    var.names <- colnames(XN)
  } else {
    var.names <- paste("X.", 1:dim(XN)[2], sep = "")
  }
  if (!is.null(dimnames(XN)[[3]])){
    x3d.names <- dimnames(XN)[[3]]
  } else {
    x3d.names <- paste("Z.", 1:dim(XN)[3], sep = "")
  }
  if (!is.null(colnames(Y))){
    yvar.names <- colnames(Y)
  } else {
    yvar.names <- paste("Y.", 1:dim(Y)[2], sep = "")
  }
  if(!center.X) center.X <- rep(0, ncol(XN)*dim(XN)[3])
  if(!center.Y) center.Y <- rep(0, ncol(Y))
  if(!scale.X) scale.X <- rep(1, ncol(XN)*dim(XN)[3])
  if(!scale.Y) scale.Y <- rep(1, ncol(Y))
  if(is.null(keepJ) | is.null(keepK) | method!="sNPLS"){
    cont_thresholding <- TRUE
    message(paste("Using continuous thresholding (", method, ")", sep=""))
  } else {
    cont_thresholding <- FALSE
    message("Using discrete L1-thresholding")
  }
  if(length(threshold_j) == 1 & ncomp > 1) threshold_j <- rep(threshold_j, ncomp)
  if(length(threshold_k) == 1 & ncomp > 1) threshold_k <- rep(threshold_k, ncomp)

  # Matrices initialization
  U <- Q <- X <- P <- NULL
  if(method == "sNPLS"){
    WsupraJ <- WsupraK <- Tm <- NULL
  } else{
    WsupraJ <- matrix(nrow=ncol(XN), ncol=ncomp)
    WsupraK <- matrix(nrow=dim(XN)[3], ncol=ncomp)
    Tm <- matrix(nrow=nrow(XN), ncol=ncomp)
  }
  Yorig <- Y
  Y <- scale(Y, center = center.Y, scale = scale.Y)
  y_center <- attr(Y, "scaled:center")
  y_scale <- attr(Y, "scaled:scale")
  B <- matrix(0, ncol = ncomp, nrow = ncomp)
  Gu <- vector("list", ncomp)
  S <- svd(Y)$d
  u <- Y[, S == max(S)]  #Column with the highest variance

  # Unfolding of XN en 2-D
  X <- unfold3w(XN)
  #Check for zero variance columns and fix them with some noise
  if(any(apply(X, 2, sd)==0)){
    X[,apply(X, 2, sd)==0] <- apply(X[,apply(X, 2, sd)==0, drop=FALSE], 2, function(x) jitter(x))
  }

  # Center and scale
  Xd <- scale(X, center = center.X, scale = scale.X)
  x_center <- attr(Xd, "scaled:center")
  x_scale <- attr(Xd, "scaled:scale")

  # Main loop for each component
  for (f in 1:ncomp) {
    if(!cont_thresholding){
      nj <- ncol(XN) - keepJ[f]
      nk <- dim(XN)[3] - keepK[f]
    }
    if(method %in% c("sNPLS-VIP", "sNPLS-SR")){
      nj <- ncol(XN)
      nk <- dim(XN)[3]
      wselj <- rep(1, nj)
      wselk <- rep(1, nk)
    }
    it = 1
    while (it < max.iteration) {
      Zrow <- crossprod(u, Xd)
      Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
      svd.z <- svd(Z)
      wsupraj <- svd.z$u[, 1]
      # L1 penalization for wsupraj
      if(method == "sNPLS"){
        if(cont_thresholding){
          wsupraj <- rel_thresholding(wsupraj, threshold_j[f])
        } else {
          if (nj != 0) {
            wsupraj <- thresholding(wsupraj, nj)
          }
        }
      }

      ##########
      wsuprak <- svd.z$v[, 1]
      # L1 penalization for wsuprak
      if(method == "sNPLS"){
        if(cont_thresholding){
          wsuprak <- rel_thresholding(wsuprak, threshold_k[f])
        } else {
          if (nk != 0) {
            wsuprak <- thresholding(wsuprak, nk)
          }
        }
      }

      if(method %in% c("sNPLS-VIP", "sNPLS-SR")){
        W <- kronecker(wsuprak, wsupraj) *  kronecker(wselk, wselj)
        tf <- Xd %*% W
        qf <- crossprod(Y, tf)/mynorm(crossprod(Y, tf))
        uf <- Y %*% qf
        Tm[,f] <- tf
        Q <- cbind(Q, qf)
        WsupraJ[,f] <- wsupraj
        WsupraK[,f] <- wsuprak
      }

      if(method == "sNPLS-VIP"){
        #VIPj
        SCE <- sum(sum(Tm[,1:f, drop=FALSE] %*% t(Q[,1:f, drop=FALSE])^2))
        VIPj <- sqrt(nrow(WsupraJ)*((WsupraJ[,f, drop=FALSE]^2)*SCE)/sum(SCE))
        VIPj_01 <- (VIPj-min(VIPj))/(max(VIPj)-min(VIPj))

        #VIPk
        SCE <- sum(sum(Tm[,1:f, drop=FALSE] %*% t(Q[,1:f, drop=FALSE])^2))
        VIPk <- sqrt(nrow(WsupraK)*((WsupraK[,f, drop=FALSE]^2)*SCE)/sum(SCE))
        VIPk_01 <- (VIPk-min(VIPk))/(max(VIPk)-min(VIPk))

        wselk <- as.numeric(VIPk_01 > threshold_k[f])
        wselj <- as.numeric(VIPj_01 > threshold_j[f])
      }

      if(method == "sNPLS-SR"){
        TM <- MASS::ginv(crossprod(Tm[,1:f, drop=FALSE])) %*% t(Tm[,1:f, drop=FALSE])
        WkM <- MASS::ginv(crossprod(WsupraK[,1:f, drop=FALSE])) %*% t(WsupraK[,1:f, drop=FALSE])
        WjM <- MASS::ginv(crossprod(WsupraJ[,1:f, drop=FALSE])) %*% t(WsupraJ[,1:f, drop=FALSE])
        Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
        #SR
        P[[f]] = t(as.matrix(Gu[[f]]) %*% t(kronecker(WsupraK[,1:f, drop=FALSE], WsupraJ[,1:f, drop=FALSE])))
        Xres <- Xd - Tm[,1:f, drop=FALSE] %*% t(P[[f]])
        xpred <- Tm[,1:f, drop=FALSE] %*% t(P[[f]])
        SSexp <- xpred^2
        SSres <- (Xd-xpred)^2
        SSexp_cube <- array(SSexp, dim=c(nrow(SSexp), nj, nk))
        SSres_cube <- array(SSres, dim=c(nrow(SSexp), nj, nk))
        SR_k <- numeric(nk)
        for(k in 1:nk){
          SR_k[k] <- sum(SSexp_cube[,,k])/sum(SSres_cube[,,k])
        }
        SR_k_01 <- (SR_k-min(SR_k))/(max(SR_k)-min(SR_k))
        SR_j <- numeric(nj)
        for(j in 1:nj){
          SR_j[j] <- sum(SSexp_cube[,j,])/sum(SSres_cube[,j,])
        }
        SR_j_01 <- (SR_j-min(SR_j))/(max(SR_j)-min(SR_j))
        wselj <- rep(1, nj)
        wselk <- rep(1, nk)
        PFcalck <- SR_k_01
        wselk <- as.numeric(PFcalck > threshold_k[f])
        PFcalcj <- SR_j_01
        wselj <- as.numeric(PFcalcj > threshold_j[f])
      }

      ##########
      if(method == "sNPLS"){
        tf <- Xd %*% kronecker(wsuprak, wsupraj)
        qf <- crossprod(Y, tf)/mynorm(crossprod(Y, tf))
        uf <- Y %*% qf
      }
      if (sum((uf - u)^2) < conver) {
        if (!silent) {
          cat(paste("Component number ", f, "\n"))
          cat(paste("Number of iterations: ", it, "\n"))
        }
        it <- max.iteration
        if(method == "sNPLS"){
          Tm <- cbind(Tm, tf)
          WsupraJ <- cbind(WsupraJ, wsupraj)
          WsupraK <- cbind(WsupraK, wsuprak)
          bf <- MASS::ginv(crossprod(Tm)) %*% t(Tm) %*% uf
          TM <- MASS::ginv(crossprod(Tm)) %*% t(Tm)
          Q <- cbind(Q, qf)
        } else {
          Tm[,f] <- tf
          WsupraJ[,f] <- wsupraj*wselj
          WsupraK[,f] <- wsuprak*wselk
          bf <- MASS::ginv(crossprod(Tm[,1:f, drop=FALSE])) %*% t(Tm[,1:f, drop=FALSE]) %*% uf
        }
        B[1:length(bf), f] <- bf
        U <- cbind(U, uf)
        if(method == "sNPLS-VIP"){
          TM <- MASS::ginv(crossprod(Tm[,1:f, drop=FALSE])) %*% t(Tm[,1:f, drop=FALSE])
          WkM <- MASS::ginv(crossprod(WsupraK[,1:f, drop=FALSE])) %*% t(WsupraK[,1:f, drop=FALSE])
          WjM <- MASS::ginv(crossprod(WsupraJ[,1:f, drop=FALSE])) %*% t(WsupraJ[,1:f, drop=FALSE])
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          P[[f]] = t(as.matrix(Gu[[f]]) %*% t(kronecker(WsupraK[,1:f, drop=FALSE], WsupraJ[,1:f, drop=FALSE])))
        }
        if(method == "sNPLS"){
          TM <- MASS::ginv(crossprod(Tm)) %*% t(Tm)
          WkM <- MASS::ginv(crossprod(WsupraK)) %*% t(WsupraK)
          WjM <- MASS::ginv(crossprod(WsupraJ)) %*% t(WsupraJ)
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          P[[f]] = t(as.matrix(Gu[[f]]) %*% t(kronecker(WsupraK, WsupraJ)))
        }

        if(method == "sNPLS"){
          Y <- Y - Tm %*% bf %*% t(qf)
        } else {
          Y <- Y - Tm[,1:f, drop=FALSE] %*% bf %*% t(qf)
        }
        S <- svd(Y)$d
        u <- Y[, S == max(S)]
      } else {
        u <- uf
        it <- it + 1
      }
    }
  }
  Yadjsc <- Tm %*% B %*% t(Q)
  Yadj <- Yadjsc * y_scale + y_center
  SqrdE <- sum((Yorig - Yadj)^2)
  rownames(WsupraJ) <- var.names
  rownames(WsupraK) <- x3d.names
  rownames(Q) <- yvar.names
  rownames(Tm) <- rownames(U) <- x.names
  colnames(Tm) <- colnames(WsupraJ) <- colnames(WsupraK) <- colnames(B) <-
    colnames(U) <- colnames(Q) <- names(Gu) <- names(P) <- paste("Comp.", 1:ncomp)
  output <- list(T = Tm, Wj = WsupraJ, Wk = WsupraK, B = B, U = U, Q = Q, P = P,
                 Gu = Gu, ncomp = ncomp, Xd=Xd, Yadj = Yadj, SqrdE = SqrdE,
                 Standarization = list(ScaleX = x_scale, CenterX = x_center,
                                       ScaleY = y_scale, CenterY = y_center),
                 Method = method)
  class(output)<-"sNPLS"
  return(output)
}

#' R-matrix from a sNPLS model fit
#'
#' @description Builds the R-matrix from a sNPLS model fit
#' @param x A sNPLS model obtained from \code{sNPLS}
#' @return Returns the R-matrix of the model, needed to compute the coefficients
Rmatrix<-function(x) {
  WsupraK <- x$Wk
  WsupraJ <- x$Wj
  R <- matrix(nrow = dim(x$Wj)[1] * dim(x$Wk)[1], ncol = x$ncomp)
  ncomp <- x$ncomp
  kroneckers<-sapply(1:x$ncomp, function(x) kronecker(WsupraK[, x], WsupraJ[, x]))
  tkroneckers<-apply(kroneckers, 2, function(x) t(x))
  R[,1] <- kroneckers[,1]
  if(ncomp>1){
    for(i in 2:ncomp){
      pi <- pi0 <- Matrix::Matrix(diag(dim(R)[1]), sparse=TRUE)
      for (j in 1:(i - 1)) {
        pi <- Matrix::Matrix(pi %*% pi0 - kroneckers[,j] %*% t(tkroneckers[,j]), sparse=TRUE)
      }
      w <- kroneckers[, i]
      pi <- pi %*% w
      R[, i] <- Matrix::as.matrix(pi)
    }
  }
  return(R)
}

#' Unfolding of three-way arrays
#'
#' @description Unfolds a three-way array into a matrix
#' @param x A three-way array
#' @return Returns a matrix with dimensions \code{dim(x)[1] x dim(x)[2]*dim(x([3]))}
unfold3w <- function(x) {
    dim(x) <- c(dim(x)[1], dim(x)[2] * dim(x)[3])
    return(x)
}

#' Cross-validation for a sNPLS model
#'
#' @description Performs cross-validation for a sNPLS model
#' @param X_npls A three-way array containing the predictors.
#' @param Y_npls A matrix containing the response.
#' @param ncomp A vector with the different number of components to test
#' @param samples Number of samples for performing random search in continuous thresholding
#' @param keepJ A vector with the different number of selected variables to test for discrete thresholding
#' @param keepK A vector with the different number of selected 'times' to test for discrete thresholding
#' @param nfold Number of folds for the cross-validation
#' @param parallel Should the computations be performed in parallel? Set up strategy first with \code{future::plan()}
#' @param method Select between sNPLS, sNPLS-SR or sNPLS-VIP
#' @param ... Further arguments passed to sNPLS
#' @return A list with the best parameters for the model and the CV error
#' @examples
#' \dontrun{
#' X_npls<-array(rpois(7500, 10), dim=c(50, 50, 3))
#'
#' Y_npls<-matrix(2+0.4*X_npls[,5,1]+0.7*X_npls[,10,1]-0.9*X_npls[,15,1]+
#' 0.6*X_npls[,20,1]- 0.5*X_npls[,25,1]+rnorm(50), ncol=1)
#' #Grid search for discrete thresholding
#' cv1<- cv_snpls(X_npls, Y_npls, ncomp=1:2, keepJ = 1:3, keepK = 1:2, parallel = FALSE)
#' #Random search for continuous thresholding
#' cv2<- cv_snpls(X_npls, Y_npls, ncomp=1:2, samples=20, parallel = FALSE)
#' }
#' @importFrom stats runif
#' @export
cv_snpls <- function(X_npls, Y_npls, ncomp = 1:3, samples=20,
                     keepJ = NULL, keepK = NULL, nfold = 10, parallel = TRUE,  method="sNPLS", ...) {

  if(parallel) message("Your parallel configuration is ", attr(future::plan(), "class")[3])
  if(!method %in% c("sNPLS", "sNPLS-SR", "sNPLS-VIP")) stop("'method' not recognized")
  if(length(dim(Y_npls)) == 3) Y_npls <- unfold3w(Y_npls)
  top <- ceiling(dim(X_npls)[1]/nfold)
  foldid <- sample(rep(1:nfold, top), dim(X_npls)[1], replace = F)
  if(is.null(keepJ) | is.null(keepK)){
    cont_thresholding <- TRUE
    message("Using continuous thresholding")
  } else {
    cont_thresholding <- FALSE
    message("Using discrete thresholding")
  }
  if(cont_thresholding){
    search.grid <- expand.grid(list(ncomp = ncomp, threshold_j=runif(samples),
                                    threshold_k=runif(samples)))
  } else {
    search.grid <- expand.grid(list(ncomp = ncomp, keepJ = keepJ,
                                    keepK = keepK))
  }
  SqrdE <- numeric()
  applied_fun <- function(y) {
    sapply(1:nfold, function(x) {
      suppressMessages(tryCatch(do.call(cv_fit, c(list(xtrain = X_npls[x != foldid, , ],
                                      ytrain = Y_npls[x != foldid, , drop = FALSE],
                                      xval = X_npls[x == foldid, , ],
                                      yval = Y_npls[x == foldid, , drop = FALSE],
                                      ncomp = y["ncomp"], method=method),
                                 list(threshold_j=y["threshold_j"])[cont_thresholding],
                                 list(threshold_k=y["threshold_k"])[cont_thresholding],
                                 list(keepJ=rep(y["keepJ"], y["ncomp"]))[!cont_thresholding],
                                 list(keepK=rep(y["keepK"], y["ncomp"]))[!cont_thresholding], ...)),
               error=function(x) NA))
    })
  }
  if (parallel) {
    cv_res <- future.apply::future_apply(search.grid, 1, applied_fun, future.seed=TRUE)
  } else {
    cv_res <- pbapply::pbapply(search.grid, 1, applied_fun)
  }
  cv_mean <- apply(cv_res, 2, function(x) mean(x, na.rm = TRUE))
  cv_se <- apply(cv_res, 2, function(x) sd(x, na.rm=TRUE)/sqrt(nfold))
  best_model <- search.grid[which.min(cv_mean), ]
  output <- list(best_parameters = best_model, cv_mean = cv_mean,
                 cv_se = cv_se, cv_grid = search.grid)
  class(output)<-"cvsNPLS"
  return(output)
}

#' Internal function for \code{cv_snpls}
#'
#' @param xtrain A three-way training array
#' @param ytrain A response training matrix
#' @param xval A three-way test array
#' @param yval A response test matrix
#' @param ncomp Number of components for the sNPLS model
#' @param threshold_j Threshold value on Wj. Scaled between [0, 1)
#' @param threshold_k Threshold value on Wk. Scaled between [0, 1)
#' @param keepJ Number of variables to keep for each component, ignored if threshold_j is provided
#' @param keepK Number of 'times' to keep for each component, ignored if threshold_k is provided
#' @param method Select between sNPLS, sNPLS-SR or sNPLS-VIP
#' @param ... Further arguments passed to sNPLS
#' @return Returns the CV mean squared error
#' @importFrom stats predict
#' @export
cv_fit <- function(xtrain, ytrain, xval, yval, ncomp, threshold_j=NULL, threshold_k=NULL, keepJ=NULL, keepK=NULL,  method, ...) {
  fit <- sNPLS(XN = xtrain, Y = ytrain, ncomp = ncomp, keepJ=keepJ, keepK=keepK, threshold_j = threshold_j,
                 threshold_k = threshold_k, silent = TRUE, method=method, ...)
  Y_pred <- predict(fit, xval)
  CVE <- sqrt(mean((Y_pred - yval)^2))
  return(CVE)
}

#' Plots for sNPLS model fits
#'
#' @description Different plots for sNPLS model fits
#' @param x A sNPLS model fit
#' @param comps Vector with the components to plot. It can be of length \code{ncomp} for types "time" and "variables" and of length 2 otherwise.
#' @param type The type of plot. One of those: "T", "U", "Wj", "Wk", "time" or "variables"
#' @param labels Should rownames be added as labels to the plot?
#' @param group Vector with categorical variable defining groups (optional)
#' @param ... Not used
#' @return A plot of the type specified in the \code{type} parameter
#' @export
plot.sNPLS <- function(x, type = "T", comps = c(1, 2), labels=TRUE, group=NULL, ...) {
  if (type == "T")
    p<-plot_T(x, comps = comps, labels=labels, group=group)
  if (type == "U")
    p<-plot_U(x, comps = comps, labels=labels, group=group)
  if (type == "Wj")
    p<-plot_Wj(x, comps = comps, labels=labels)
  if (type == "Wk")
    p<-plot_Wk(x, comps = comps, labels=labels)
  if (type == "time")
    p<-plot_time(x, comps = comps)
  if (type == "variables")
    p<-plot_variables(x, comps = comps)
  p
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param labels Should rownames be added as labels to the plot?
#' @param group Vector with categorical variable defining groups
#' @return A plot of the T matrix of a sNPLS model fit
plot_T <- function(x, comps, labels, group=NULL){
  df <- data.frame(x$T)[comps]
  if(!is.null(group)) df <- data.frame(df, group=as.factor(group))
  names(df)[1:2] <- paste("Comp.", comps, sep="")
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=names(df)[1], y=names(df)[2])) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = 0,
                                                                      lty=2) +
    ggplot2::geom_hline(yintercept = 0, lty=2) + ggplot2::xlab(names(df)[1]) +
    ggplot2::ylab(names(df)[2]) + ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                                                 axis.title=ggplot2::element_text(size=12,face="bold"))
  if(labels) p1 <- p1 + ggrepel::geom_text_repel(ggplot2::aes(label=rownames(df)))
  if(!is.null(group)) p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color=group))
  p1
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param labels Should rownames be added as labels to the plot?
#' @param group Vector with categorical variable defining groups
#' @return A plot of the U matrix of a sNPLS model fit
plot_U <- function(x, comps, labels, group=NULL){
  df <- data.frame(x$U)[comps]
  if(!is.null(group)) df <- data.frame(df, group=as.factor(group))
  names(df)[1:2] <- paste("Comp.", comps, sep="")
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=names(df)[1], y=names(df)[2])) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = 0,
                                                                      lty=2) +
    ggplot2::geom_hline(yintercept = 0, lty=2) + ggplot2::xlab(names(df)[1]) +
    ggplot2::ylab(names(df)[2]) + ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                                                 axis.title=ggplot2::element_text(size=12,face="bold"))
  if(labels) p1 <- p1 + ggrepel::geom_text_repel(ggplot2::aes(label=rownames(df)))
  if(!is.null(group)) p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color=group))
  p1
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param labels Should rownames be added as labels to the plot?
#' @return A plot of Wj coefficients
plot_Wj <- function(x, comps, labels){
  df <- data.frame(x$Wj)[comps]
  names(df) <- paste("Comp.", comps, sep="")
  var_names_zero <- rownames(df)
  var_names_zero[rowSums(df)==0]<- ""
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=names(df)[1], y=names(df)[2])) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = 0,
                                                                      lty=2) +
    ggplot2::geom_hline(yintercept = 0, lty=2) + ggplot2::xlab(names(df)[1]) +
    ggplot2::ylab(names(df)[2]) + ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                                                 axis.title=ggplot2::element_text(size=12,face="bold"))
  if(labels) p1 <- p1 + ggrepel::geom_text_repel(ggplot2::aes(label=var_names_zero), size=3)
  p1
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param labels Should rownames be added as labels to the plot?
#' @return A plot of the Wk coefficients
plot_Wk <- function(x, comps, labels){
  df <- data.frame(x$Wk)[comps]
  names(df) <- paste("Comp.", comps, sep="")
  var_names_zero <- rownames(df)
  var_names_zero[rowSums(df)==0]<- ""
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=names(df)[1], y=names(df)[2])) +
    ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = 0,
                                                                      lty=2) +
    ggplot2::geom_hline(yintercept = 0, lty=2) + ggplot2::xlab(names(df)[1]) +
    ggplot2::ylab(names(df)[2]) + ggplot2::theme(axis.text=ggplot2::element_text(size=12),
                                                 axis.title=ggplot2::element_text(size=12,face="bold"))
  if(labels) p1 <- p1 + ggrepel::geom_text_repel(ggplot2::aes(label=var_names_zero), size=4)
  p1
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector with the components to plot
#' @return A plot of Wk coefficients for each component
plot_time <- function(x, comps){
  df <- clickR::forge(data.frame(row=1:nrow(x$Wk), x$Wk), affixes=as.character(comps),
                      var.name="Component")
  ggplot2::ggplot(df, ggplot2::aes_string(x="row", y="Comp..", color="Component")) +
    ggplot2::geom_line(lwd=1.05) + ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0, alpha=0.2) +
    ggplot2::xlab("Time (index)") + ggplot2::ylab("Wk") +
    ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                   axis.title = ggplot2::element_text(size=12,face="bold"))
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector with the components to plot
#' @return A plot of Wj coefficients for each component
plot_variables <- function(x, comps){
  df <- clickR::forge(data.frame(row=1:nrow(x$Wj), x$Wj), affixes=as.character(comps),
                      var.name="Component")
  ggplot2::ggplot(df, ggplot2::aes_string(x="row", y="Comp..", color="Component")) +
    ggplot2::geom_line(lwd=1.05) + ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0, alpha=0.2) +
    ggplot2::xlab("Variable (index)") + ggplot2::ylab("Wj") +
    ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                   axis.title=ggplot2::element_text(size=12,face="bold"))
}

#' Predict for sNPLS models
#'
#' @description Predict function for sNPLS models
#' @param object A sNPLS model fit
#' @param newX A three-way array containing the new data
#' @param rescale Should the prediction be rescaled to the original scale?
#' @param ... Further arguments passed to \code{predict}
#' @return A matrix with the predictions
#' @export
predict.sNPLS <- function(object, newX, rescale = TRUE, ...) {
  newX <- unfold3w(newX)
  # Centrado y escalado
  #Xstd <- t((t(newX) - object$Standarization$CenterX)/object$Standarization$ScaleX)
  #Xstd <- sweep(sweep(newX, 2, object$Standarization$CenterX), 2, object$Standarization$ScaleX, "/")
  Xstd <- scale(newX, center=object$Standarization$CenterX, scale=object$Standarization$ScaleX)
  R <- Rmatrix(object)
  Bnpls <- R %*% object$B %*% t(object$Q)
  Yval <- Xstd %*% Bnpls
  if (rescale) {
    Yval <- Yval * object$Standarization$ScaleY + object$Standarization$CenterY
  }
  return(Yval)
}

#' Fitted method for sNPLS models
#'
#' @description Fitted method for sNPLS models
#' @param object A sNPLS model fit
#' @param ... Further arguments passed to \code{fitted}
#' @return Fitted values for the sNPLS model
#' @export
fitted.sNPLS <- function(object, ...){
  return(object$Yadj)
}

#' Plot cross validation results for sNPLS objects
#'
#' @description Plot function for visualization of cross validation results for sNPLS models
#' @param x A cv_sNPLS object
#' @param ... Not used
#' @return A facet plot with the results of the cross validation
#' @export
plot.cvsNPLS <- function(x, ...) {
  cont_thresholding <- names(x$cv_grid)[2] == "threshold_j"
  df_grid <- data.frame(threshold_j=x$cv_grid[,2], threshold_k=x$cv_grid[,3], CVE=x$cv_mean, Ncomp=paste("Ncomp =", x$cv_grid$ncomp, sep=" "))
  if(!cont_thresholding) names(df_grid)[c(1, 2)] <- c("KeepJ", "KeepK")
  if(cont_thresholding){
    ggplot2::ggplot(df_grid, ggplot2::aes_string(x="threshold_j", y="CVE"))+ggplot2::geom_point()+ggplot2::geom_smooth()+ggplot2::facet_grid(cut(threshold_k,10) ~ Ncomp)+
      ggplot2::theme_bw()
  } else {
      ggplot2::ggplot(df_grid, ggplot2::aes_string(x="KeepJ", y="CVE"))+ggplot2::geom_point()+ggplot2::geom_line()+ggplot2::facet_grid(KeepK ~ Ncomp)+
        ggplot2::theme_bw()
  }
}

#' Coefficients from a sNPLS model
#'
#' @description Extract coefficients from a sNPLS model
#' @param object A sNPLS model fit
#' @param as.matrix Should the coefficients be presented as matrix or vector?
#' @param ... Further arguments passed to \code{coef}
#' @return A matrix (or vector) of coefficients
#' @export
coef.sNPLS <- function(object, as.matrix = FALSE, ...) {
  R <- Rmatrix(object)
  Bnpls <- R %*% object$B %*% t(object$Q)
  colnames(Bnpls) <- paste("Estimate", colnames(Bnpls))
  if (as.matrix){
    dim(Bnpls) <- c(dim(object$Wj)[1], dim(object$Wk)[1], dim(Bnpls)[2])
    rownames(Bnpls) <- rownames(object$Wj)
    colnames(Bnpls) <- rownames(object$Wk)
  }
  return(Bnpls)
}

#' Repeated cross-validation for sNPLS models
#'
#' @description Performs repeated cross-validatiodn and represents results in a plot
#' @param X_npls A three-way array containing the predictors.
#' @param Y_npls A matrix containing the response.
#' @param ncomp A vector with the different number of components to test
#' @param samples Number of samples for performing random search in continuous thresholding
#' @param keepJ A vector with the different number of selected variables to test in discrete thresholding
#' @param keepK A vector with the different number of selected 'times' to test in discrete thresholding
#' @param nfold Number of folds for the cross-validation
#' @param times Number of repetitions of the cross-validation
#' @param parallel Should the computations be performed in parallel? Set up strategy first with \code{future::plan()}
#' @param method Select between sNPLS, sNPLS-SR or sNPLS-VIP
#' @param ... Further arguments passed to cv_snpls
#' @return A density plot with the results of the cross-validation and an (invisible) \code{data.frame} with these results
#' @importFrom stats var
#' @export
repeat_cv <- function(X_npls, Y_npls, ncomp = 1:3, samples=20, keepJ=NULL, keepK=NULL, nfold = 10, times=30, parallel = TRUE, method="sNPLS", ...){
  if(!method %in% c("sNPLS", "sNPLS-SR", "sNPLS-VIP")) stop("'method' not recognized")
  if(parallel) message("Your parallel configuration is ", attr(future::plan(), "class")[3])
  if(is.null(keepJ) | is.null(keepK)){
    cont_thresholding <- TRUE
    message("Using continuous thresholding")
  } else {
    cont_thresholding <- FALSE
    message("Using discrete thresholding")
  }
  if(parallel){
    rep_cv<-future.apply::future_sapply(1:times, function(x) suppressMessages(cv_snpls(X_npls, Y_npls, ncomp=ncomp, parallel = FALSE, nfold = nfold, samples=samples, keepJ=keepJ, keepK=keepK, method=method, ...)), future.seed=TRUE)
  } else {
    rep_cv<-pbapply::pbreplicate(times, suppressMessages(cv_snpls(X_npls, Y_npls, ncomp=ncomp, parallel = FALSE, nfold = nfold, samples=samples, keepJ=keepJ, keepK=keepK, method=method, ...)))
  }
  resdata<-data.frame(ncomp=sapply(rep_cv[1,], function(x) x[[1]]), threshold_j=sapply(rep_cv[1,], function(x) x[[2]]),
                      threshold_k=sapply(rep_cv[1,], function(x) x[[3]]))
  if(!cont_thresholding) names(resdata)[c(2, 3)] <- c("keepJ", "keepK")
  class(resdata)<-c("repeatcv", "data.frame")
  return(resdata)
}

#' Density plot for repeat_cv results
#'
#' @description Plots a grid of slices from the 3-D kernel denity estimates of the repeat_cv function
#' @param x A repeatcv object
#' @param ... Further arguments passed to plot
#' @return A grid of slices from a 3-D density plot of the results of the repeated cross-validation
#' @importFrom grDevices colorRampPalette
#' @importFrom stats ftable density setNames
#' @export
plot.repeatcv <- function(x, ...){
  x.old <- x
  x<-x[,sapply(x, function(x) var(x)>0), drop=FALSE]
  if(ncol(x) < ncol(x.old)) warning(paste("\n", colnames(x.old)[!colnames(x.old) %in% colnames(x)], "is constant at", x.old[1,colnames(x.old)[!colnames(x.old) %in% colnames(x)]]))
  if(ncol(x) == 1){
    densities <- density(x[,1])
    df_grid <- setNames(data.frame(densities$x, densities$y), c(colnames(x), "density"))
    p <- ggplot2::ggplot(x, ggplot2::aes_string(x=colnames(x)))+ggplot2::geom_density(color="gray", fill="gray", alpha=0.3)+ggplot2::theme_classic()
  } else{
    H.pi <- ks::Hpi(x)
    fhat <- ks::kde(x, H=H.pi, compute.cont=TRUE, gridsize = rep(151, ncol(x)))
    if(ncol(x) == 3){
      ncomp_values <- sapply(sort(unique(fhat$x[,1])), function(x) which.min(abs(fhat$eval.points[[1]]-x)))
      positions <- as.data.frame(fhat$x)
      positions$Ncomp <- factor(positions$ncomp)
      df_grid <- setNames(expand.grid(fhat$eval.points[[2]], fhat$eval.points[[3]]), names(positions)[c(2,3)])
      combl <- nrow(df_grid)
      df_grid <- df_grid[rep(1:nrow(df_grid), length(ncomp_values)),]
      df_grid$density <- unlist(lapply(ncomp_values, function(x) as.numeric(matrix(fhat$estimate[x,,], ncol=1))))
      df_grid$Ncomp <- factor(rep(sort(unique(positions$ncomp)), each=combl))
      p <- ggplot2::ggplot(df_grid, ggplot2::aes_string(names(df_grid)[1], names(df_grid)[2], fill="density"))+ggplot2::geom_raster()+
        ggplot2::scale_fill_gradientn(colours =colorRampPalette(c("white", "blue", "red"))(10))+ggplot2::theme_classic()+
        ggplot2::geom_count(inherit.aes = FALSE, ggplot2::aes_string(x=names(df_grid)[1], y=names(df_grid)[2]), data=positions) +ggplot2::facet_grid(~Ncomp)
      if(names(df_grid)[1] == "threshold_j"){
        p <- p + ggplot2::scale_x_continuous(limits=c(0, 1)) + ggplot2::scale_y_continuous(limits=c(0, 1))
      } else {
        p <- p + ggplot2::scale_x_continuous(breaks=if(round(diff(range(df_grid$keepJ)))<=10) round(max(0, min(df_grid$keepJ)):max(df_grid$keepJ)) else round(seq(max(0, min(df_grid$keepJ)), max(df_grid$keepJ), by= ceiling(max(df_grid$keepJ)/20)*2)))+
          ggplot2::scale_y_continuous(breaks=if(round(diff(range(df_grid$keepK)))<=10) round(max(0, min(df_grid$keepK)):max(df_grid$keepK)) else round(seq(max(0, min(df_grid$keepK)), max(df_grid$keepK), by= ceiling(max(df_grid$keepK)/20)*2)))+
          ggplot2::scale_size_area(breaks=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))],
                                   labels=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))])
      }
    } else {
      positions <- as.data.frame(fhat$x)
      df_grid <- expand.grid(V1=fhat$eval.points[[1]], V2=fhat$eval.points[[2]])
      names(df_grid)<-colnames(positions)
      df_grid$density <- as.numeric(matrix(fhat$estimate, ncol=1))
      p <- ggplot2::ggplot(df_grid, ggplot2::aes_string(colnames(df_grid)[1], colnames(df_grid)[2], fill="density"))+ggplot2::geom_raster()+
        ggplot2::scale_fill_gradientn(colours =colorRampPalette(c("white", "blue", "red"))(10))+ggplot2::theme_classic()+
        ggplot2::geom_count(inherit.aes = FALSE, ggplot2::aes_string(x=colnames(df_grid)[1], y=colnames(df_grid)[2]), data=positions)
      if("threshold_j" %in% names(df_grid) | "threshold_k" %in% names(df_grid)){
        p <- p + ggplot2::scale_x_continuous(limits=c(0, 1)) + ggplot2::scale_y_continuous(limits=c(0, 1))
      } else {
        p <- p + ggplot2::scale_x_continuous(breaks=if(round(diff(range(df_grid[,1])))<=10) round(max(0, min(df_grid[,1])):max(df_grid[,1])) else round(seq(max(0, min(df_grid[,1])), max(df_grid[,1]), by= ceiling(max(df_grid[,1])/20)*2)))+
          ggplot2::scale_y_continuous(breaks=if(round(diff(range(df_grid[,2])))<=10) round(max(0, min(df_grid[,2])):max(df_grid[,2])) else round(seq(max(0, min(df_grid[,2])), max(df_grid[,2]), by= ceiling(max(df_grid[,2])/20)*2)))+
          ggplot2::scale_size_continuous(breaks=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))],
                                         labels=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))])
      }
    }
  }
  print(p)
  data.frame(x.old[1, !colnames(x.old) %in% colnames(x), drop=FALSE], df_grid[which.max(df_grid$density),])
}

#' Summary for sNPLS models
#'
#' @description Summary of a sNPLS model fit
#' @param object A sNPLS object
#' @param ... Further arguments passed to summary.default
#' @return A summary inclunding number of components, squared error and coefficients of the fitted model
#' @importFrom stats coef
#' @export
summary.sNPLS <- function(object, ...){
  cat(object$Method, "model with", object$ncomp, "components and squared error of", round(object$SqrdE,3), "\n", "\n")
  cat("Coefficients: \n")
  round(coef(object, as.matrix=TRUE), 3)
}

#' Compute Selectivity Ratio for a sNPLS model
#'
#' @description Estimates Selectivity Ratio for the different components of a sNPLS model fit
#' @param model A sNPLS model
#' @return A list of data.frames, each of them including the computed Selectivity Ratios for each variable
#' @export
SR <- function(model){
  output <- lapply(1:model$ncomp, function(f){
    Xres <- model$Xd - model$T[,1:f, drop=FALSE] %*% t(model$P[[f]])
    Xpred <- model$T[,1:f, drop=FALSE] %*% t(model$P[[f]])
    SSexp <- Xpred^2
    SSres <- (model$Xd-Xpred)^2
    SSexp_cube <- array(SSexp, dim=c(nrow(SSexp), nrow(model$Wj), nrow(model$Wk)))
    SSres_cube <- array(SSres, dim=c(nrow(SSexp), nrow(model$Wj), nrow(model$Wk)))
    SR_k <- sapply(1:nrow(model$Wk), function(x) sum(SSexp_cube[,,x])/sum(SSres_cube[,,x]))
    SR_j <- sapply(1:nrow(model$Wj), function(x) sum(SSexp_cube[,x,])/sum(SSres_cube[,x,]))
    list(SR_k=round(SR_k,3),
         SR_j=round(SR_j,3))
  })
  SR_j <- do.call("cbind", sapply(output, function(x) x[2]))
  SR_k <- do.call("cbind", sapply(output, function(x) x[1]))
  rownames(SR_j) <- rownames(model$Wj)
  rownames(SR_k) <- rownames(model$Wk)
  colnames(SR_j) <- colnames(SR_k) <- paste("Component ", 1:model$ncomp, sep="")
  list(SR_j=SR_j, SR_k=SR_k)
}
