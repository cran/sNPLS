#' Fit a sNPLS model
#'
#' @description Fits a N-PLS regression model imposing a L1 penalization on \code{wj} and \code{wk} matrices
#' @param XN A three-way array containing the predictors.
#' @param Y A matrix containing the response.
#' @param ncomp Number of components in the projection
#' @param conver Convergence criterion
#' @param max.iteration Maximum number of iterations
#' @param keepJ Number of variables to keep for each component
#' @param keepK Number of 'times' to keep for each component
#' @param silent Show output?
#' @return A fitted sNPLS model
#' @references C. A. Andersson and R. Bro. The N-way Toolbox for MATLAB Chemometrics & Intelligent Laboratory Systems. 52 (1):1-4, 2000.
#' @references Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation. Journal of Multivariate Analysis 99, 1015-1034
#' @examples
#' X_npls<-array(rpois(7500, 10), dim=c(50, 50, 3))
#'
#' Y_npls<-matrix(2+0.4*X_npls[,5,1]+0.7*X_npls[,10,1]-0.9*X_npls[,15,1]+
#' 0.6*X_npls[,20,1]- 0.5*X_npls[,25,1]+rnorm(50), ncol=1)
#'
#' fit<-sNPLS(X_npls, Y_npls, ncomp=3, keepJ = rep(2,3) , keepK = rep(1,3))
#' @importFrom stats predict sd
#' @export
sNPLS <- function(XN, Y, ncomp = 2, conver = 1e-16, max.iteration = 10000,
                  keepJ = rep(ncol(XN), ncomp), keepK = rep(rev(dim(XN))[1], ncomp),
                  silent = F) {

    mynorm <- function(x) sqrt(sum(diag(crossprod(x))))
    if (length(dim(XN)) != 3)
        stop("'XN' is not a three-way array")
    if (!is.null(rownames(XN)))
        y.names <- x.names <- rownames(XN) else y.names <- x.names <- 1:dim(XN)[1]
    if (!is.null(colnames(XN)))
        var.names <- colnames(XN) else var.names <- paste("X.", 1:dim(XN)[2], sep = "")
    if (!is.null(dimnames(XN)[[3]]))
        x3d.names <- dimnames(XN)[[3]] else x3d.names <- paste("Z.", 1:dim(XN)[3], sep = "")
    if (!is.null(colnames(Y)))
        yvar.names <- colnames(Y) else yvar.names <- paste("Y.", 1:dim(Y)[2], sep = "")

    # Matrices initialization
    Tm <- U <- Q <- WsupraJ <- WsupraK <- X <- P <- NULL
    Yorig <- Y
    Y <- scale(Y)
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
    Xd <- scale(X)
    x_center <- attr(Xd, "scaled:center")
    x_scale <- attr(Xd, "scaled:scale")

    # Main loop for each component
    for (f in 1:ncomp) {
        nj <- ncol(XN) - keepJ[f]
        nk <- dim(XN)[3] - keepK[f]
        it = 1
        while (it < max.iteration) {
            Zrow <- crossprod(u, Xd)
            Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
            svd.z <- svd(Z)
            wsupraj <- svd.z$u[, 1]
            # L1 penalization for wsupraj
            if (nj != 0) {
                wsupraj <- ifelse(abs(wsupraj) > abs(wsupraj[order(abs(wsupraj))][nj]),
                                  (abs(wsupraj) - abs(wsupraj[order(abs(wsupraj))][nj])) *
                                    sign(wsupraj), 0)
            }
            ##########
            wsuprak <- svd.z$v[, 1]
            # L1 penalization for wsuprak
            if (nk != 0) {
                wsuprak <- ifelse(abs(wsuprak) > abs(wsuprak[order(abs(wsuprak))][nk]),
                                  (abs(wsuprak) - abs(wsuprak[order(abs(wsuprak))][nk])) *
                                    sign(wsuprak), 0)
            }
            ##########
            tf <- Xd %*% kronecker(wsuprak, wsupraj)
            qf <- crossprod(Y, tf)/mynorm(crossprod(Y, tf))
            uf <- Y %*% qf
            if (sum((uf - u)^2) < conver) {
                if (!silent) {
                  print(paste("Component number ", f))
                  print(paste("Number of iterations: ", it))
                }
                it <- max.iteration
                Tm <- cbind(Tm, tf)
                WsupraJ <- cbind(WsupraJ, wsupraj)
                WsupraK <- cbind(WsupraK, wsuprak)
                bf <- MASS::ginv(crossprod(Tm)) %*% t(Tm) %*% uf
                B[1:length(bf), f] <- bf
                Q <- cbind(Q, qf)
                U <- cbind(U, uf)
                TM <- MASS::ginv(crossprod(Tm)) %*% t(Tm)
                WkM <- MASS::ginv(crossprod(WsupraK)) %*% t(WsupraK)
                WjM <- MASS::ginv(crossprod(WsupraJ)) %*% t(WsupraJ)
                Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
                # Xd <- Xd - Tm%*%Gu[[f]] %*% t(kronecker(WsupraK,WsupraJ))
                P[[f]] = t(as.matrix(Gu[[f]]) %*% t(kronecker(WsupraK, WsupraJ)))
                Y <- Y - Tm %*% bf %*% t(qf)
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
                   Gu = Gu, ncomp = ncomp, Yadj = Yadj, SqrdE = SqrdE,
                   Standarization = list(ScaleX = x_scale, CenterX = x_center,
                                         ScaleY = y_scale, CenterY = y_center))
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
#' @param keepJ A vector with the different number of selected variables to test
#' @param keepK A vector with the different number of selected 'times' to test
#' @param nfold Number of folds for the cross-validation
#' @param parallel Should the computations be performed in parallel?
#' @param free_cores If parallel computations are performed how many cores are left unused
#' @return A list with the best parameters for the model and the CV error
#' @examples
#' \dontrun{
#' X_npls<-array(rpois(7500, 10), dim=c(50, 50, 3))
#'
#' Y_npls<-matrix(2+0.4*X_npls[,5,1]+0.7*X_npls[,10,1]-0.9*X_npls[,15,1]+
#' 0.6*X_npls[,20,1]- 0.5*X_npls[,25,1]+rnorm(50), ncol=1)
#'
#' cv1<- cv_snpls(X_npls, Y_npls, ncomp=1:2, keepJ = 1:3, keepK = 1:2, parallel = FALSE)
#' }
#' @export
cv_snpls <- function(X_npls, Y_npls, ncomp = 1:3, keepJ = 1:ncol(X_npls),
                     keepK = 1:dim(X_npls)[3], nfold = 10, parallel = TRUE, free_cores = 2) {
    if (parallel & (parallel::detectCores()>1)) {
        cl <- parallel::makeCluster(max(2, parallel::detectCores() - free_cores))
        parallel::clusterExport(cl, list(deparse(substitute(X_npls)),
                                         deparse(substitute(Y_npls))))
        parallel::clusterCall(cl, function() require(sNPLS))
    }
    top <- ceiling(dim(X_npls)[1]/nfold)
    foldid <- sample(rep(1:nfold, top), dim(X_npls)[1], replace = F)
    search.grid <- expand.grid(list(ncomp = ncomp, keepJ = keepJ, keepK = keepK))
    SqrdE <- numeric()
    applied_fun <- function(y) {
        sapply(1:nfold, function(x) {
            tryCatch(cv_fit(xtrain = X_npls[x != foldid, , ],
                            ytrain = Y_npls[x != foldid, , drop = FALSE],
                            xval = X_npls[x == foldid, , ],
                            yval = Y_npls[x == foldid, , drop = FALSE],
                            ncomp = y["ncomp"],
                            keepJ = rep(y["keepJ"], y["ncomp"]),
                            keepK = rep(y["keepK"], y["ncomp"])),
                     error=function(x) NA)
          })
    }
    if (parallel) {
        cv_res <- parallel::parApply(cl, search.grid, 1, applied_fun)
        parallel::stopCluster(cl)
    } else cv_res <- pbapply::pbapply(search.grid, 1, applied_fun)
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
#' @param keepJ Number of variables to keep for each component
#' @param keepK Number of 'times' to keep for each component
#' @return Returns the CV mean squared error
#' @export
cv_fit <- function(xtrain, ytrain, xval, yval, ncomp, keepJ, keepK) {
  fit <- sNPLS(XN = xtrain, Y = ytrain, ncomp = ncomp, keepJ = keepJ,
               keepK = keepK, silent = TRUE)
  Y_pred <- predict(fit, xval)
  CVE <- sqrt(mean((Y_pred - yval)^2))
  return(CVE)
}

#' Plots for sNPLS model fits
#'
#' @description Different plots for sNPLS model fits
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param type The type of plot. One of those: "T", "U", "Wj", "Wk", "time" or "variables"
#' @param ... Options passed to \code{plot}
#' @return A plot of the type specified in the \code{type} parameter
#' @importFrom graphics abline matplot plot text layout par plot.new
#' @export
plot.sNPLS <- function(x, type = "T", comps = c(1, 2), ...) {
  old.mar<- par()$mar
  if (type == "T")
    plot_T(x, comps = comps, ...)
  if (type == "U")
    plot_U(x, comps = comps, ...)
  if (type == "Wj")
    plot_Wj(x, comps = comps, ...)
  if (type == "Wk")
    plot_Wk(x, comps = comps, ...)
  if (type == "time")
    plot_time(x, comps = comps, ...)
  if (type == "variables")
    plot_variables(x, comps = comps, ...)
  par(mar=old.mar)
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param xlim Limits of the X axis
#' @param ylim Limits of the Y axis
#' @param ... Options passed to \code{plot}
#' @return A plot of the T matrix of a sNPLS model fit
plot_T <- function(x,
                   comps,
                   xlim = c(min(x$T[, comps[1]])-diff(range(x$T[, comps[1]]))/10,
                            max(x$T[, comps[1]])+diff(range(x$T[, comps[1]]))/10),
                   ylim = c(min(x$T[, comps[2]])-diff(range(x$T[, comps[2]]))/10,
                            max(x$T[, comps[2]])+diff(range(x$T[, comps[2]]))/10),
                   ...){
  plot(x$T[, comps[1]], x$T[, comps[2]],
       pch  = 16,
       xlab = colnames(x$T)[comps[1]],
       ylab = colnames(x$T)[comps[2]],
       ylim = ylim,
       xlim = xlim, ...)
  abline(h = 0, v = 0, lty = 2)
  text(x$T[, comps[1]], x$T[, comps[2]], labels = rownames(x$T),
       pos = plotrix::thigmophobe(x$T[, comps[1]], x$T[, comps[2]]))
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param xlim Limits of the X axis
#' @param ylim Limits of the Y axis
#' @param ... Options passed to \code{plot}
#' @return A plot of the U matrix of a sNPLS model fit
plot_U <- function(x,
                   comps,
                   ylim = c(min(x$U[, comps[2]])-diff(range(x$U[, comps[2]]))/10,
                            max(x$U[, comps[2]])+diff(range(x$U[, comps[2]]))/10),
                   xlim = c(min(x$U[, comps[1]])-diff(range(x$U[, comps[1]]))/10,
                            max(x$U[, comps[1]])+diff(range(x$U[, comps[1]]))/10),
                   ...){
  plot(x$U[, comps[1]], x$U[, comps[2]],
       pch  = 16,
       xlab = colnames(x$U)[comps[1]], ylab = colnames(x$U)[comps[2]],
       ylim = ylim,
       xlim = xlim, ...)
  abline(h = 0, v = 0, lty = 2)
  text(x$U[, comps[1]], x$U[, comps[2]], labels = rownames(x$U),
       pos = plotrix::thigmophobe(x$U[, comps[1]], x$U[, comps[2]]))
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param xlim Limits of the X axis
#' @param ylim Limits of the Y axis
#' @param ... Options passed to \code{plot}
#' @return A plot of Wj coefficients
plot_Wj <- function(x,
                    comps,
                    xlim = c(min(x$Wj[, comps[1]]) - diff(range(x$Wj[, comps[1]]))/10,
                            max(x$Wj[, comps[1]]) + diff(range(x$Wj[, comps[1]]))/10),
                    ylim= c(min(x$Wj[, comps[2]]) - diff(range(x$Wj[, comps[2]]))/10,
                            max(x$Wj[, comps[2]]) + diff(range(x$Wj[, comps[2]]))/10),
                    ...){
  plot(x$Wj[, comps[1]], x$Wj[, comps[2]],
       pch  = 16,
       xlab = colnames(x$Wj)[comps[1]],
       ylab = colnames(x$Wj)[comps[2]],
       ylim = ylim,
       xlim = xlim, ...)
  abline(h = 0, v = 0, lty = 2)
  text(x$Wj[, comps[1]][x$Wj[, comps[1]] != 0 | x$Wj[, comps[2]] != 0],
       x$Wj[, comps[2]][x$Wj[,comps[1]] != 0 | x$Wj[, comps[2]] != 0],
       pos = tryCatch(plotrix::thigmophobe(x$Wj[, comps[1]][x$Wj[, comps[1]] != 0 |
                                                              x$Wj[, comps[2]] != 0],
                                  x$Wj[, comps[2]][x$Wj[, comps[1]] != 0 |
                                                     x$Wj[, comps[2]] != 0]),
                      error=function(x) 1),
       labels = which(x$Wj[, comps[1]] != 0 | x$Wj[, comps[2]] != 0))
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector of length two with the components to plot
#' @param xlim Limits of the X axis
#' @param ylim Limits of the Y axis
#' @param ... Options passed to \code{plot}
#' @return A plot of the Wk coefficients
plot_Wk <- function(x,
                    comps,
                    xlim = c(min(x$Wk[, comps[1]]) - diff(range(x$Wk[, comps[1]]))/10,
                                     max(x$Wk[, comps[1]])+diff(range(x$Wk[, comps[1]]))/10),
                    ylim = c(min(x$Wk[, comps[2]]) - diff(range(x$Wk[, comps[2]]))/10,
                             max(x$Wk[, comps[2]]) + diff(range(x$Wk[, comps[2]]))/10),
                    ...){
  plot(x$Wk[, comps[1]], x$Wk[, comps[2]],
       pch  = 16,
       xlab = colnames(x$Wk)[comps[1]],
       ylab = colnames(x$Wk)[comps[2]],
       ylim = ylim,
       xlim = xlim, ...)
  abline(h = 0, v = 0, lty = 2)
  text(x$Wk[, comps[1]][x$Wk[, comps[1]] != 0 | x$Wk[, comps[2]] != 0],
       x$Wk[, comps[2]][x$Wk[,comps[1]] != 0 | x$Wk[, comps[2]] != 0],
       pos = tryCatch(plotrix::thigmophobe(x$Wk[, comps[1]][x$Wk[, comps[1]] != 0 |
                                                              x$Wk[, comps[2]] != 0],
                                  x$Wk[, comps[2]][x$Wk[, comps[1]] != 0 |
                                                     x$Wk[, comps[2]] != 0]),
                      error=function(x) 1),
       labels = which(x$Wk[, comps[1]] != 0 | x$Wk[, comps[2]] != 0))
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector with the components to plot
#' @param xlab X-axis label
#' @param ... Options passed to \code{plot}
#' @return A plot of Wk coefficients for each component
#' @importFrom graphics legend
plot_time <- function(x, comps, xlab="Time", ...) {
  layout(cbind(1,2), widths = c(5,1))
  par(mar=c(5.1, 5.1, 4.1, 1))
  matplot(1:dim(x$Wk[,comps])[1], x$Wk[,comps], type = "l", xlab = xlab,
          ylab = "Wk", col = 1:5, lty=1:5, ...)
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("left", colnames(x$Wk)[comps], col=1:5, lty=1:5, lwd=2)
  layout(1)
}

#' Internal function for \code{plot.sNPLS}
#'
#' @param x A sNPLS model fit
#' @param comps A vector with the components to plot
#' @param xlab X-axis label
#' @param ... Options passed to \code{plot}
#' @return A plot of Wj coefficients for each component
#' @importFrom graphics legend
plot_variables <- function(x, comps, xlab="Variables", ...) {
  layout(cbind(1,2), widths = c(5,1))
  par(mar=c(5.1, 5.1, 4.1, 1))
  matplot(1:dim(x$Wj[,comps])[1], x$Wj[,comps], type = "l", xlab = xlab,
          ylab = "Wj", col = 1:5, lty=1:5, ...)
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("left", colnames(x$Wk)[comps], col=1:5, lty=1:5, lwd=2)
  layout(1)
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
  Xstd <- t((t(newX) - object$Standarization$CenterX)/object$Standarization$ScaleX)
  R <- Rmatrix(object)
  Bnpls <- R %*% object$B %*% t(object$Q)
  Yval <- Xstd %*% Bnpls
  if (rescale) {
    Yval <- Yval * object$Standarization$ScaleY + object$Standarization$CenterY
  }
  return(Yval)
}

#' Plot cross validation results for sNPLS objects
#'
#' @description Plot function for visualization of cross validation results for sNPLS models
#' @param x A cv_sNPLS object
#' @param facets Chose between a facet plot or a 3-D scatter plot
#' @param ... Arguments passed to \code{car::scatter3d}
#' @return A 3D scatter plot with the results of the cross validation
#' @export
plot.cvsNPLS <- function(x, facets=TRUE, ...) {
  df_grid <- data.frame(KeepJ=x$cv_grid$keepJ, KeepK=paste("KeepK =", x$cv_grid$keepK, sep=" "), CVE=x$cv_mean, Ncomp=paste("Ncomp =", x$cv_grid$ncomp, sep=" "))
  if(facets){
    ggplot2::ggplot(df_grid, ggplot2::aes_string(x="KeepJ", y="CVE"))+ggplot2::geom_line()+ggplot2::facet_grid(KeepK ~ Ncomp)+
      ggplot2::scale_x_continuous(breaks=if(round(diff(range(df_grid$KeepJ)))<=10) round(max(0, min(df_grid$KeepJ)):max(df_grid$KeepJ)) else round(seq(max(0, min(df_grid$KeepJ)), max(df_grid$KeepJ), by= ceiling(max(df_grid$KeepJ)/20)*2)))+
      ggplot2::theme_bw()
  } else{
    car::scatter3d(x$cv_grid$keepJ, x$cv_mean, x$cv_grid$keepK, groups = if (length(unique(x$cv_grid$ncomp)) > 1) factor(x$cv_grid$ncomp) else NULL,
                   surface = TRUE, fit = "smooth", axis.col = c("black", "black", "black"), xlab = "KeepJ", ylab = "CVE", zlab = "KeepK", parallel = FALSE, ...)
    rgl::grid3d(c("x", "y", "z"), col = "gray", lwd = 1, lty = 1, n = 5)
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
  colnames(Bnpls) <- "Estimate"
  if (as.matrix){
    dim(Bnpls) <- c(dim(object$Wj)[1], dim(object$Wk)[1])
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
#' @param keepJ A vector with the different number of selected variables to test
#' @param keepK A vector with the different number of selected 'times' to test
#' @param nfold Number of folds for the cross-validation
#' @param parallel Should the computations be performed in parallel?
#' @param free_cores If parallel computations are performed how many cores are left unused
#' @param ... Currently not used
#' @param times Number of repetitions of the cross-validation
#' @return A density plot with the results of the cross-validation and an (invisible) \code{data.frame} with these results
#' @importFrom stats var
#' @export
repeat_cv<-function(X_npls, Y_npls, ncomp = 1:3, keepJ = 1:ncol(X_npls), keepK = 1:dim(X_npls)[3],
                    nfold = 10, parallel = TRUE, free_cores = 2, times=30, ...){
  if(parallel & (parallel::detectCores()>1)){
    cl <- parallel::makeCluster(max(2, parallel::detectCores() - free_cores))
    parallel::clusterExport(cl, list(deparse(substitute(X_npls)), deparse(substitute(Y_npls))))
    parallel::clusterCall(cl, function() require(sNPLS))
    rep_cv<-parallel::parSapply(cl, 1:times, function(x) cv_snpls(X_npls, Y_npls, ncomp=ncomp, keepJ = keepJ, keepK = keepK, parallel = FALSE, nfold = nfold))
    parallel::stopCluster(cl)
  } else {
    rep_cv<-pbapply::pbreplicate(times, cv_snpls(X_npls, Y_npls, ncomp=ncomp, keepJ = keepJ, keepK = keepK, parallel = FALSE, nfold = nfold))
  }
  resdata<-data.frame(ncomp=sapply(rep_cv[1,], function(x) x[[1]]), keepJ=sapply(rep_cv[1,], function(x) x[[2]]),
                      keepK=sapply(rep_cv[1,], function(x) x[[3]]))
  class(resdata)<-c("repeatcv", "data.frame")
  return(resdata)
}

#' Density plot for repat_cv results
#'
#' @description Plots a grid of slices from the 3-D kernel denity estimates of the repeat_cv function
#' @param x A repeatcv object
#' @param ... Further arguments passed to plot
#' @return A grid of slices from of a 3-D density plot of the results of the repeated cross-validation
#' @importFrom grDevices colorRampPalette
#' @importFrom stats ftable
#' @export
plot.repeatcv <- function(x, ...){
  x.old <- x
  x<-x[,sapply(x, function(x) var(x)>0), drop=FALSE]
  if(ncol(x) < ncol(x.old)) warning(paste("\n", colnames(x.old)[!colnames(x.old) %in% colnames(x)], "is constant at", x.old[1,colnames(x.old)[!colnames(x.old) %in% colnames(x)]]))
  if(ncol(x) == 1){
    ggplot2::ggplot(x, ggplot2::aes_string(x=colnames(x)))+ggplot2::geom_density(color="gray", fill="gray", alpha=0.3)+ggplot2::theme_classic()
  } else{
    H.pi <- ks::Hpi(x)
    fhat <- ks::kde(x, H=H.pi, compute.cont=TRUE, gridsize = rep(151, ncol(x)))
    if(ncol(x) == 3){
      ncomp_values <- sapply(sort(unique(fhat$x[,1])), function(x) which.min(abs(fhat$eval.points[[1]]-x)))
      positions <- as.data.frame(fhat$x)
      positions$Ncomp <- paste("Ncomp =", positions$ncomp)
      df_grid <- expand.grid(keepJ=fhat$eval.points[[2]], keepK=fhat$eval.points[[3]])
      combl <- nrow(df_grid)
      df_grid <- df_grid[rep(1:nrow(df_grid), length(ncomp_values)),]
      df_grid$density <- unlist(lapply(ncomp_values, function(x) as.numeric(matrix(fhat$estimate[x,,], ncol=1))))
      df_grid$Ncomp <- rep(paste("Ncomp =", sort(unique(positions$ncomp))), each=combl)
      ggplot2::ggplot(df_grid, ggplot2::aes_string("keepJ", "keepK", fill="density"))+ggplot2::geom_raster()+
        ggplot2::scale_fill_gradientn(colours =colorRampPalette(c("white", "blue", "red"))(10))+ggplot2::theme_classic()+
        ggplot2::geom_count(inherit.aes = FALSE, ggplot2::aes_string(x="keepJ", y="keepK"), data=positions) +ggplot2::facet_grid(~Ncomp)+
        ggplot2::scale_x_continuous(breaks=if(round(diff(range(df_grid$keepJ)))<=10) round(max(0, min(df_grid$keepJ)):max(df_grid$keepJ)) else round(seq(max(0, min(df_grid$keepJ)), max(df_grid$keepJ), by= ceiling(max(df_grid$keepJ)/20)*2)))+
        ggplot2::scale_y_continuous(breaks=if(round(diff(range(df_grid$keepK)))<=10) round(max(0, min(df_grid$keepK)):max(df_grid$keepK)) else round(seq(max(0, min(df_grid$keepK)), max(df_grid$keepK), by= ceiling(max(df_grid$keepK)/20)*2)))+
        ggplot2::scale_size_area(breaks=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))],
                                 labels=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))])
    } else {
      positions <- as.data.frame(fhat$x)
      df_grid <- expand.grid(V1=fhat$eval.points[[1]], V2=fhat$eval.points[[2]])
      names(df_grid)<-colnames(positions)
      df_grid$density <- as.numeric(matrix(fhat$estimate, ncol=1))
      ggplot2::ggplot(df_grid, ggplot2::aes_string(colnames(df_grid)[1], colnames(df_grid)[2], fill="density"))+ggplot2::geom_raster()+
        ggplot2::scale_fill_gradientn(colours =colorRampPalette(c("white", "blue", "red"))(10))+ggplot2::theme_classic()+
        ggplot2::geom_count(inherit.aes = FALSE, ggplot2::aes_string(x=colnames(df_grid)[1], y=colnames(df_grid)[2]), data=positions)+
        ggplot2::scale_x_continuous(breaks=if(round(diff(range(df_grid[,1])))<=10) round(max(0, min(df_grid[,1])):max(df_grid[,1])) else round(seq(max(0, min(df_grid[,1])), max(df_grid[,1]), by= ceiling(max(df_grid[,1])/20)*2)))+
        ggplot2::scale_y_continuous(breaks=if(round(diff(range(df_grid[,2])))<=10) round(max(0, min(df_grid[,2])):max(df_grid[,2])) else round(seq(max(0, min(df_grid[,2])), max(df_grid[,2]), by= ceiling(max(df_grid[,2])/20)*2)))+
        ggplot2::scale_size_continuous(breaks=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))],
                                 labels=if(length(sort(unique(as.numeric(ftable(positions))))[-1])<=7) sort(unique(as.numeric(ftable(positions))))[-1] else (1:max(as.numeric(ftable(positions))))[round(seq(1, max(as.numeric(ftable(positions))), length.out = 7))])
    }
  }
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
  cat("sNPLS model with", object$ncomp, "components and squared error of", round(object$SqrdE,3), "\n", "\n")
  cat("Coefficients: \n")
  round(coef(object, as.matrix=TRUE),3)
}
