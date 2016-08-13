
list.of.packages <- c("ggplot2", "circular", "gridExtra", "grid", 
  "latex2exp", "shinydashboard", "shiny", "knitr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, 
  "Package"])]
if (length(new.packages)) install.packages(new.packages)


library(knitr)
library("shiny")
library(shinydashboard)
library("foreign")
library(circular)
library(ggplot2)
library("gridExtra")
library("grid")
library(latex2exp)
#' Consensus Model
#' 
#' Function that fit a consensus model for angular variables and its \code{print} method.
#' 
#' The control argument is a list that can supply any of the following components: 
#' \describe{
#'   \item{\code{pginit}}{ The approximate number of points on the grid of possible initial beta values tried when 
#'                         \code{initbeta} is not given. The default is 1000, which runs quickly. A large value of 
#'                         \code{pginit} makes the function slower.  }
#'   \item{\code{maxiter}}{ The maximum number of iterations. The default is 1000.  }
#'   \item{\code{mindiff}}{ The minimum difference between two max cosine to be reached. It defines the convergence criteria:
#'                          if the difference between the max cosine for the updated parameters values and the max
#'                          cosine for the parameters values at the previous iteration is below \code{mindiff}, convergence is
#'                          reached. The default is 0.000001.   }
#' }
#' 
#' @param formula  A formula with the dependent angle on the left of the ~ operator and terms specifying
#'                 the explanatory variables on the right. These terms must be written \code{x:z}, where
#'                 \code{x} is an explanatory angle which relative importance migth depend on the
#'                 positive variable \code{z}. It is not mandatory to specify a \code{z} variable for each
#'                 explanatory angle. For \code{model='simplified'}, the first explanatory angle listed is 
#'                 the reference direction (if a \code{z} variable was specified for this angle, it is ignored).
#' @param data  An optional data frame, list or environment 
#'              containing the variables in the model formula. If not found in data, the variables are taken from 
#'              \code{environment(formula)}, typically the environment from which \code{angular} is called. 
#' @param model  A character string, either \code{'complete'} for the complete model with an intercept (the default) 
#'               or \code{'simplified'} for the simplified model without an intercept.
#' @param initparam  A numerical vector, initial values for the parameters. The default is to use the best initial
#'                   values found among some values tried on a grid of possible values for the parameters.
#' @param control A list of control parameters. See Details.
#' 
#' @return \item{MaxLL}{ the maximum value of the log likeliohood } 
#' @return \item{parameters}{ the parameter estimates and their standard errors (obtained from two definitions) } 
#' @return \item{varcov1}{ the estimated variance covariance matrix for the parameter estimates (obtained from the first definition) }
#' @return \item{varcov2}{ the estimated variance covariance matrix for the parameter estimates (obtained from the second definition) }
#' @return \item{parambeta}{ the beta parameter estimates and their standard errors (obtained by linearization) } 
#' @return \item{varcovbeta1}{ the estimated variance covariance matrix for the beta parameter estimates (obtained by linearization) }
#' @return \item{varcovbeta2}{ the estimated variance covariance matrix for the beta parameter estimates (Sandwhich form) }
#' @return \item{autocorr}{ the autocorrelation of the residuals \eqn{\sin(y_i-\mu_i)}{sin(yi-mui)}  }
#' @return \item{iter.detail}{ the iteration details }
#' @return \item{converge}{ an indicator of convergence }
#' @return \item{call}{ the function call }
#' 
#' @author Sophie Baillargeon, Louis-Paul Rivest and Aurélien Nicosia
consensus <- function(formula, data, model = "simplified", weights = NULL, 
  initbeta = NULL, control = list()) {
  
  call <- mfcall <- match.call()
  model <- model[1]
  
  
  
  ### information from formula and data
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  # useful objects
  nobs <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  nterms <- length(nomterms)
  p <- if ("simplified" == model) 
    nterms - 1 else nterms
  nparam <- if ("simplified" == model) 
    p + 1 else p + 2
  # paramname <- paste0('kappa', 0:p)
  paramname <- nomterms
  if ("complete" == model) 
    paramname <- c(paramname, paste0("beta", p + 1))
  # first column = response variable
  y <- as.vector(mf[, 1])
  # explanatory variables
  noms <- strsplit(nomterms, split = ":")
  noms <- do.call(rbind, noms)
  if ("simplified" == model) {
    x0 <- mf[, noms[1, 1]]
    noms <- noms[-1, , drop = FALSE]
  }
  matx <- as.matrix(mf[, noms[, 1], drop = FALSE])
  if (ncol(noms) == 1) {
    matz <- matrix(1, ncol = ncol(matx), nrow = nrow(matx))  # all z are 1
  } else {
    matz <- as.matrix(mf[, noms[, 2], drop = FALSE])
    matz[, noms[, 2] == noms[, 1]] <- 1  # unspecified z are 1
  }
  weight = rep(1, nobs) * (is.null(weights)) + (!is.null(weights)) * 
    weights
  
  ### log-likelihood function
  LL <- function(param) {
    angleref <- if ("simplified" == model) 
      x0 else rep(param[p + 2], nobs)
    # length of the vector
    sinmu <- param[1] * sin(angleref) + (matz * sin(matx)) %*% 
      param[2:(p + 1)]
    cosmu <- param[1] * cos(angleref) + (matz * cos(matx)) %*% 
      param[2:(p + 1)]
    long <- as.vector(sqrt(sinmu^2 + cosmu^2))
    # predicted value form the model
    mui <- as.vector(atan2(sinmu, cosmu))
    # log likelihood
    term1 <- param[1] * cos(y - angleref) + (matz * cos(y - 
      matx)) %*% param[2:(p + 1)]
    # LL <- sum(weight*term1) - sum(weight*log(besselI(long, 0,
    # expon.scaled = FALSE)))
    LL <- sum(term1) - sum(log(besselI(long, 0, expon.scaled = FALSE)))
    
    
    list(LL = LL, long = long, mui = mui)
  }
  # Function that update parameter of the log-likelihood
  # function
  paramUpdate <- function(paramk, long, mui) {
    angleref <- if ("simplified" == model) 
      x0 else rep(paramk[p + 2], nobs)
    matx0 <- cbind(angleref, matx)
    matz0 <- cbind(rep(1, nobs), matz)
    # score vector
    Along <- as.vector(besselI(long, 1, expon.scaled = FALSE)/besselI(long, 
      0, expon.scaled = FALSE))
    matu <- matz0 * (cos(y - matx0) - cos(matx0 - mui) * 
      Along)
    if ("complete" == model) {
      
      matu <- cbind(matu, paramk[1] * sin(y - angleref) - 
        sin(mui - angleref) * Along)
    }
    vecs <- colSums(matu)
    names(vecs) <- paramname
    # Fisher information matrix
    Xc <- matz0 * cos(matx0 - mui)
    Xs <- matz0 * sin(matx0 - mui)
    if ("complete" == model) {
      
      Xc <- cbind(Xc, paramk[1] * sin(mui - paramk[p + 
        2]))
      Xs <- cbind(Xs, paramk[1] * cos(mui - paramk[p + 
        2]))
    }
    Dc <- diag(1 - Along/long - Along^2, nrow = nobs, ncol = nobs)
    Ds <- diag(Along/long, nrow = nobs, ncol = nobs)
    matI <- t(Xc) %*% Dc %*% Xc + t(Xs) %*% Ds %*% Xs
    colnames(matI) <- rownames(matI) <- paramname
    # update of parameters
    dparam <- as.vector(solve(matI, vecs))
    paramk1 <- paramk + dparam
    list(paramk1 = paramk1, dparam = dparam, matu = matu, 
      matI = matI)
  }
  # initial values of parameters beta we try 10 000 differents
  # possible values
  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 
      1000 else control$pginit
    pg <- round(pginit^(1/nparam))
    possparam <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, 
      pg + 2)]), p + 1)
    if ("complete" == model) 
      possparam[[nparam]] <- seq(0, 2 * pi, length.out = pg + 
        2)[-c(1, pg + 2)]
    possVal <- cbind(expand.grid(possparam), NA)
    colnames(possVal) <- c(paramname, "LL")
    maxLL <- function(param) LL(param = param)$LL
    possVal[, nparam + 1] <- apply(possVal[, 1:nparam], 1, 
      maxLL)
    paramk <- unlist(possVal[which.max(possVal[, nparam + 
      1]), 1:nparam])
  } else {
    if (length(initbeta) != nparam) 
      stop("for the requested model, 'initparam' must be of length ", 
        nparam)
    paramk <- initbeta
  }
  # log likelihood value with respect to initial parameter
  calcul <- LL(param = paramk)
  maxLLk <- calcul$LL
  long <- calcul$long
  mui <- calcul$mui
  # Initialization
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 
    1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 
    1e-06 else control$mindiff
  conv <- FALSE
  # Initialisation of the matrix of informations during the
  # iterations
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 
    3)
  colnames(iter.detail) <- c(paramname, "maxLL", "iter", "nitersh")
  iter.detail[1, ] <- c(paramk, maxLLk, iter, iter.sh)
  # start of the fit
  while (!conv && iter <= maxiter) {
    # update of the parameters
    maj <- paramUpdate(paramk = paramk, long = long, mui = mui)
    paramk1 <- maj$paramk1
    dparam <- maj$dparam
    # computation of the log-likelihood
    calcul <- LL(param = paramk1)
    maxLLk1 <- calcul$LL
    long <- calcul$long
    mui <- calcul$mui
    # if the criteria as decrease then do step halving
    iter.sh <- 0
    while (maxLLk1 < maxLLk) {
      iter.sh <- iter.sh + 1
      paramk1 <- paramk + dparam/(2^iter.sh)
      calcul <- LL(param = paramk1)
      maxLLk1 <- calcul$LL
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) 
        break
    }
    # Does the criteria increase more than mindiff?
    if (maxLLk1 < maxLLk) {
      conv <- FALSE
      warning("the algorithm did no converge, it failed to maximize the log likelihood")
      break
    } else {
      conv <- if (maxLLk1 - maxLLk > mindiff) 
        FALSE else TRUE
      paramk <- paramk1
      maxLLk <- maxLLk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(paramk, maxLLk, iter, 
        iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning("the algorithm did not converge, the maximum number of iterations was reached")
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }
  
  ### Computation of standard errors and Fisher information
  ### matrix for the final parameters.
  if (maxLLk == maxLLk1) {
    maj <- paramUpdate(paramk = paramk, long = long, mui = mui)
  }
  matu <- maj$matu
  matI <- maj$matI
  # parametric estimation of the covariance matrix
  v1 <- solve(matI)
  # non parametric estimation
  mid <- matrix(0, ncol = nparam, nrow = nparam)
  for (i in 1:nobs) {
    mid <- mid + t(matu[i, , drop = FALSE]) %*% matu[i, , 
      drop = FALSE]
  }
  v2 <- v1 %*% mid %*% v1
  
  ### Results for the betas
  paramb <- paramk[2:(p + 1)]/paramk[1]
  matDeriv <- rbind(-paramk[2:(p + 1)]/paramk[1]^2, diag(1/paramk[1], 
    nrow = p, ncol = p))
  vb <- t(matDeriv) %*% v1[1:(p + 1), 1:(p + 1)] %*% matDeriv
  vb2 <- t(matDeriv) %*% v2[1:(p + 1), 1:(p + 1)] %*% matDeriv
  # names(paramb) <- colnames(vb) <- rownames(vb)
  # <-colnames(vb2) <- rownames(vb2)<- paste0('beta', 1:p)
  names(paramb) <- colnames(vb) <- rownames(vb) <- colnames(vb2) <- rownames(vb2) <- paramname[-1]
  
  ### Output
  zvalue <- abs(paramk)/sqrt(diag(v2))
  p <- round(2 * pnorm(abs(paramk)/sqrt(diag(v2)), lower.tail = FALSE), 
    5)
  parameters <- cbind(paramk, sqrt(diag(v2)), zvalue, p)
  # colnames(parameters) <- c('estimate', paste('stderr', 1:2,
  # sep=''))
  colnames(parameters) <- c("estimate", "stderr", "z value", 
    "P(|z|>.)")
  rownames(parameters) <- paramname
  
  
  zvaluebeta <- abs(paramb)/sqrt(diag(vb2))
  pbeta <- round(2 * pnorm(abs(paramb)/sqrt(diag(vb2)), lower.tail = FALSE), 
    5)
  parambeta <- cbind(paramb, sqrt(diag(vb2)), zvaluebeta, pbeta)
  # colnames(parameters) <- c('estimate', paste('stderr', 1:2,
  # sep=''))
  colnames(parambeta) <- c("estimate", "stderr", "z value", 
    "P(|z|>.)")
  rownames(parambeta) <- names(paramb)
  # parambeta <- cbind(paramb, sqrt(diag(vb)), sqrt(diag(vb2)))
  # colnames(parambeta) <- c('estimate', paste('stderr', 1:2,
  # sep = ''))
  res <- sin(y - mui)
  autocorr <- acf(res, plot = FALSE)
  out <- list(MaxLL = maxLLk, parameters = parameters, varcov1 = v1, 
    varcov2 = v2, parambeta = parambeta, varcovbeta1 = vb, 
    varcovbeta2 = vb2, autocorr = autocorr, matx = matx, 
    matz = matz, y = y, long = long, mui = mui, iter.detail = iter.detail, 
    call = call)
  class(out) <- "consensus"
  out
}


#' @param x An object, produced by the \code{\link{consensus}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
print.consensus <- function(x, ...) {
  cat("\nMaximum log-likelihood :", x$MaxLL, "\n")
  cat("\nKappa Parameters:\n")
  print.default(x$parameters, print.gap = 2, quote = FALSE, 
    right = TRUE, ...)
  cat("\n")
  # cat('\nBeta Parameters:\n') print.default(x$parambeta,
  # print.gap = 2, quote = FALSE, right=TRUE, ...) cat('\n')
  invisible(x)
}



gg_qq <- function(fitted, resid, distribution = "norm", ..., 
  line.estimate = NULL, conf = 0.95, labels = names(x)) {
  
  # resid = ScaledResid(object)
  
  # fitted vs predicted
  df = data.frame(predicted.mean = fitted, residual = resid)
  p1 = ggplot(df, aes(x = predicted.mean, y = residual)) + 
    geom_point(colour = "blue") + geom_abline(intercept = 0, 
    slope = 0) + labs(x = "Predicted Mean", y = " Residual")
  
  
  # histogram
  residual = data.frame(resid = resid)
  gg <- ggplot(residual, aes(x = resid))
  range.residual = range(residual$resid)
  gg <- gg + geom_histogram(binwidth = (range.residual[2] - 
    range.residual[1])/10, colour = "black", aes(y = ..density.., 
    fill = ..count..))
  # gg <- gg + scale_fill_gradient('Count')
  p2 <- gg + stat_function(fun = dnorm, color = "red", args = list(mean = mean(residual$resid), 
    sd = sd(residual$resid)))
  
  
  
  
  
  # QQ-plot
  x = resid
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if (is.null(line.estimate)) {
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if (!is.null(labels)) {
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, 
      labels[ord], "")
  }
  
  p3 <- ggplot(df, aes(x = z, y = ord.x)) + geom_point(colour = "Blue") + 
    geom_abline(intercept = coef[1], slope = coef[2]) + geom_ribbon(aes(ymin = lower, 
    ymax = upper), alpha = 0.2) + labs(x = "Theoretical Quantiles", 
    y = " Residuals")
  if (!is.null(labels)) 
    p <- p + geom_text(aes(label = label))
  
  
  # boxplot of residuals
  residual$boxplot = as.factor(rep(1, length(residual$resid)))
  p4 = ggplot(residual, aes(x = boxplot, y = resid)) + geom_boxplot(colour = "blue")
  
  
  
  multiplot(p1, p2, p3, p4, cols = 2)
  
  
  grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, top = "Diagnostic plots for scaled residuals")
  
  
  
  
}
hist.resid <- function(resid) {
  residual = data.frame(resid = resid)
  gg <- ggplot(residual, aes(x = resid))
  gg <- gg + geom_histogram(bins = 12, colour = "black", aes(y = ..density.., 
    fill = ..count..))
  gg <- gg + scale_fill_gradient("Count")
  gg <- gg + stat_function(fun = dnorm, color = "red", args = list(mean = mean(residual$resid), 
    sd = sd(residual$resid)))
  
  return(gg)
}


multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel ncol: Number of columns of plots nrow:
    # Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
      ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
      ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain
      # this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
        layout.pos.col = matchidx$col))
    }
  }
}


# Fonction pour ajuster le mod?le de r?gression angulaire
# acceptant des covariables Date de la derni?re mise ? jour :
# 27 juin 2013 Auteur: Sophie Baillargeon Nicosia Aurelien:
# Ajout de l'estimateur de variance parametrique (v_0)
# (01-05-2014) Variance comme dans l'article (04-02-2015)
# Toujours mod?le simplifi? (04-02-2015)


#' Angular Regression Model
#' 
#' Function that fit a regression model for angular variables and its \code{print} method.
#' 
#' The control argument is a list that can supply any of the following components: 
#' \describe{
#'   \item{\code{pginit}}{ The approximate number of points on the grid of possible initial beta values tried when 
#'                         \code{initbeta} is not given. The default is 1000, which runs quickly. A large value of 
#'                         \code{pginit} makes the function slower.  }
#'   \item{\code{maxiter}}{ The maximum number of iterations. The default is 1000.  }
#'   \item{\code{mindiff}}{ The minimum difference between two max cosine to be reached. It defines the convergence criteria:
#'                          if the difference between the max cosine for the updated parameters values and the max
#'                          cosine for the parameters values at the previous iteration is below \code{mindiff}, convergence is
#'                          reached. The default is 0.000001.   }
#' }
#' 
#' @param formula  A formula with the dependent angle on the left of the ~ operator and terms specifying
#'                 the explanatory variables on the right. These terms must be written \code{x:z}, where
#'                 \code{x} is an explanatory angle which relative importance migth depend on the
#'                 positive variable \code{z}. It is not mandatory to specify a \code{z} variable for each
#'                 explanatory angle. For \code{model='simplified'}, the first explanatory angle listed is 
#'                 the reference direction (if a \code{z} variable was specified for this angle, it is ignored).
#' @param data  An optional data frame, list or environment 
#'              containing the variables in the model formula. If not found in data, the variables are taken from 
#'              \code{environment(formula)}, typically the environment from which \code{angular} is called. 
#' @param model  A character string, either \code{'complete'} for the complete model with an intercept (the default) 
#'               or \code{'simplified'} for the simplified model without an intercept.
#' @param initbeta  A numerical vector, initial values for the beta parameters. The default is to use the best initial
#'                  values found among some values tried on a grid of possible values for beta.
#' @param control A list of control parameters. See Details.
#' 
#' @return \item{MaxCosine}{ the maximum value of the cosine } 
#' @return \item{parameters}{ the parameter estimates and their standard errors (obtained from four definitions) } 
#' @return \item{autocorr}{ the autocorrelation of the residuals \eqn{\sin(y_i-\mu_i)}{sin(yi-mui)}  }
#' @return \item{varcov1}{ the estimated variance covariance matrix for the parameter estimates (obtained from the first definition) }
#' @return \item{varcov2}{ the estimated variance covariance matrix for the parameter estimates (obtained from the second definition) }
#' @return \item{varcov3}{ the estimated variance covariance matrix for the parameter estimates (obtained from the third definition) }
#' @return \item{varcov4}{ the estimated variance covariance matrix for the parameter estimates (obtained from the fourth definition) }
#' @return \item{mui}{ the vector of the predicted mean angles }
#' @return \item{long}{ the vector of the predicted concentrations }
#' @return \item{iter.detail}{ the iteration details }
#' @return \item{converge}{ an indicator of convergence }
#' @return \item{call}{ the function call }
#' 
#' @author Sophie Baillargeon, Louis-Paul Rivest and Aurelien Nicosia 
#' @references L.-P. Rivest, T. Duchesne, A. Nicosia & D. Fortin. A general angular regression model for the analysis of data on animal movement in ecology. Journal of the Royal Statistical Society, series C, to appear
#' @examples 
#' library(circular)
#' data(wind) # example on wind directions
#' n=length(wind)
#' dat.hom=data.frame(wind.t=wind[-c(1,2)],wind.t_1=wind[-c(1,n)],wind.t_2=wind[-c(n-1,n)])
#' an=angular(wind.t~wind.t_1+wind.t_2,data=dat.hom)
#' an
angular <- function(formula, data, model = "simplified", initbeta = NULL, 
  control = list()) {
  
  # library(circular)
  
  call <- mfcall <- match.call()
  model <- model[1]
  
  # browser() # pendant le d?veloppement seulement
  
  ### Afin de tirer les info des arguments formula et data
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  # Quelques objets utiles
  nobs <- nrow(mf)
  nomterms <- attr(attr(mf, "terms"), "term.labels")
  paramname <- nomterms
  nterms <- length(nomterms)
  p <- if ("simplified" == model) 
    nterms - 1 else nterms
  nparam <- if ("simplified" == model) 
    p else p + 1
  # betaname <- paste('beta', 1:nparam, sep='')
  betaname <- paramname[-1]
  # Premi?re colonne = variable r?ponse
  y <- mf[, 1]
  # Objets relatifs aux variables explicatives
  noms <- strsplit(nomterms, split = ":")
  noms <- do.call(rbind, noms)
  if ("simplified" == model) {
    x0 <- mf[, noms[1, 1]]
    noms <- noms[-1, , drop = FALSE]
  }
  matx <- as.matrix(mf[, noms[, 1], drop = FALSE])
  if (ncol(noms) == 1) {
    matz <- matrix(1, ncol = ncol(matx), nrow = nrow(matx))  # tous les z sont des colonnes de uns
  } else {
    matz <- as.matrix(mf[, noms[, 2], drop = FALSE])
    matz[, noms[, 2] == noms[, 1]] <- 1  # les z non sp?cifi?s sont des colonnes de uns
  }
  
  ### Ajustement du mod?le Fonction pour le calcul du crit?re
  ### max-cosine ? maximiser pas en argument car ne change pas :
  ### y, x0, matx, p, model
  maxcos <- function(beta) {
    sinmu <- sin(if ("simplified" == model) x0 else beta[p + 
      1]) + as.vector((matz * sin(matx)) %*% beta[1:p])  # ?quation 1, 1e ?l?ment du vecteur
    cosmu <- cos(if ("simplified" == model) x0 else beta[p + 
      1]) + as.vector((matz * cos(matx)) %*% beta[1:p])  # ?quation 1, 2e ?l?ment du vecteur
    long <- sqrt(sinmu^2 + cosmu^2)  # ?quation 2
    mui <- atan2(sinmu, cosmu)  # ?quation 3
    maxcos <- mean(cos(y - mui))  # indice ? maximiser
    list(maxcos = maxcos, long = long, mui = mui)
  }
  # Fonction pour la mise ? jour des beta pas en argument car
  # ne change pas : y, matx, p, model
  betaUpdate <- function(betak, long, mui) {
    matd <- cbind(matz * sin(matx - mui), if ("simplified" == 
      model) 
      NULL else cos(betak[p + 1] - mui))/long  # Matrice Xi de l'?quation 11
    res <- sin(y - mui)  # sin des r?sidus dans l'?quation 11, not? d
    dbeta <- as.vector(solve(t(matd) %*% matd, t(matd) %*% 
      res))
    # Probleme inversion matrice: on utilise l<inverse generalise
    # dbeta<-as.vector(ginv(t(matd)%*%matd)%*%t(matd)%*%res)
    betak1 <- betak + dbeta
    list(betak1 = betak1, dbeta = dbeta, matd = matd, res = res)
  }
  # Valeurs initiales des param?tres beta Quelques valeurs
  # possibles sur une grille sont essay?es Objectif : essayer
  # environ 10000 vecteurs beta (c'est trop long d'en essayer
  # plus)
  if (is.null(initbeta)) {
    pginit <- if (is.null(control$pginit)) 
      1000 else control$pginit
    pg <- round(pginit^(1/nparam))
    possbeta <- rep(list(seq(-1, 1, length.out = pg + 2)[-c(1, 
      pg + 2)]), p)
    if ("complete" == model) 
      possbeta[[p + 1]] <- seq(0, 2 * pi, length.out = pg + 
        2)[-c(1, pg + 2)]
    possVal <- cbind(expand.grid(possbeta), NA)
    colnames(possVal) <- c(betaname, "maxcosine")
    # for (i in 1:nrow(possVal)) { # la boucle est un peu plus
    # lente que le apply possVal[i,nparam+1] <-
    # maxcos(beta=unlist(possVal[i, 1:nparam]))$maxcos }
    maxcos1 <- function(beta) maxcos(beta = beta)$maxcos
    possVal[, nparam + 1] <- apply(possVal[, 1:nparam, drop = FALSE], 
      1, maxcos1)
    betak <- unlist(possVal[which.max(possVal[, nparam + 
      1]), 1:nparam])
  } else {
    if (length(initbeta) != nparam) 
      stop("for the requested model, 'initparam' must be of length ", 
        nparam)
    betak <- initbeta
  }
  # Calcul du max-cosine pour le beta initial
  calcul <- maxcos(beta = betak)
  maxcosk <- calcul$maxcos
  long <- calcul$long
  mui <- calcul$mui
  # Initialisation de variables pour la boucle
  iter <- iter.sh <- 0
  maxiter <- if (is.null(control$maxiter)) 
    1000 else control$maxiter
  mindiff <- if (is.null(control$mindiff)) 
    1e-06 else control$mindiff
  conv <- FALSE
  # Initialisation de la matrice des d?tails des it?rations
  iter.detail <- matrix(NA, nrow = maxiter + 1, ncol = nparam + 
    3)
  colnames(iter.detail) <- c(betaname, "maxcosine", "iter", 
    "nitersh")
  iter.detail[1, ] <- c(betak, maxcosk, iter, iter.sh)
  # D?but de la boucle d'ajustement du mod?le
  while (!conv && iter <= maxiter) {
    # Mise ? jour de beta
    maj <- betaUpdate(betak = betak, long = long, mui = mui)
    betak1 <- maj$betak1
    dbeta <- maj$dbeta
    # Calcul du crit?re ? optimiser pour beta mis ? jour (betak1)
    calcul <- maxcos(beta = betak1)
    maxcosk1 <- calcul$maxcos
    long <- calcul$long
    mui <- calcul$mui
    # Si le crit?re a diminu?, faire du step halving
    iter.sh <- 0
    while (maxcosk1 < maxcosk) {
      iter.sh <- iter.sh + 1
      betak1 <- betak + dbeta/(2^iter.sh)
      calcul <- maxcos(beta = betak1)
      maxcosk1 <- calcul$maxcos
      long <- calcul$long
      mui <- calcul$mui
      if (iter.sh >= maxiter) 
        break
    }
    # Est-ce que le crit?re a augment? de plus de mindiff?
    if (maxcosk1 < maxcosk) {
      conv <- FALSE
      warning("the algorithm did no converge, it failed to maximize the max-cosine")
      break
    } else {
      conv <- if (maxcosk1 - maxcosk > mindiff) 
        FALSE else TRUE
      betak <- betak1
      maxcosk <- maxcosk1
      iter <- iter + 1
      iter.detail[iter + 1, ] <- c(betak, maxcosk, iter, 
        iter.sh)
    }
  }
  if (iter > maxiter + 1) {
    warning("the algorithm did not converge, the maximum number of iterations was reached")
  } else {
    iter.detail <- iter.detail[1:(iter + 1), , drop = FALSE]
  }
  
  ### Calcul d'erreurs type pour les param?tres Calcul de matd
  ### (matrice Xi) et res (sin(y-mui)) pour les param?tres finaux
  ### Si non convergence parce que max-cosine diminue, on a d?j?
  ### tout.  mais si convergence ou non converge car maxiter
  ### atteint, on n'a pas matd et res pour la toute derni?re mise
  ### ? jour de beta.
  if (maxcosk == maxcosk1) {
    maj <- betaUpdate(betak = betak, long = long, mui = mui)
  }
  matd <- maj$matd
  res <- maj$res
  Akappa <- mean(maxcosk)
  kappahat <- circular::A1inv(Akappa)
  # Calcul des matrices de variance-covariance
  dlong <- cbind(cos(matx - mui), if ("simplified" == model) 
    NULL else -sin(betak[p + 1] - mui))  # chaque ligne de cette matrice est un dlongi transpos?
  DDX <- t(dlong) %*% diag(res/long) %*% matd
  A <- (-1/nobs) * (t(matd) %*% diag(cos(y - mui)) %*% matd + 
    DDX + t(DDX))
  Ainv <- solve(A)
  v0 <- solve(t(matd) %*% matd)/(Akappa * kappahat)
  v1 <- (1/(nobs * (nobs - 1))) * (Ainv %*% (t(matd) %*% diag(res^2) %*% 
    matd) %*% Ainv)
  A1 <- (-1/nobs) * (t(matd) %*% diag(cos(y - mui)) %*% matd)
  A1inv <- solve(A1)
  # A1inv<-ginv(A1)
  v2 <- (1/(nobs * (nobs - 1))) * (A1inv %*% (t(matd) %*% diag(res^2) %*% 
    matd) %*% A1inv)
  Mautocorr <- t(matd[-nobs, ]) %*% diag(res[-nobs] * sin(y[-1] - 
    mui[-nobs])) %*% matd[-1, ]
  v3 <- (1/(nobs * (nobs - 1))) * (Ainv %*% ((t(matd) %*% diag(res^2) %*% 
    matd) + Mautocorr + t(Mautocorr)) %*% Ainv)
  v4 <- solve(t(matd) %*% matd) * mean(res^2)/(maxcosk^2)
  # v4 <- ginv(t(matd) %*% matd)*mean(res^2)/(maxcosk^2)
  
  ### Sortie des r?sultats parameters <-
  ### cbind(betak,sqrt(diag(v0)), sqrt(diag(v1)), sqrt(diag(v2)),
  ### sqrt(diag(v3)), sqrt(diag(v4))) parameters <-
  ### cbind(betak,sqrt(diag(v0)), sqrt(diag(v1)))
  zvalue <- abs(betak)/sqrt(diag(v0))
  p <- round(2 * pnorm(abs(betak)/sqrt(diag(v0)), lower.tail = FALSE), 
    5)
  parameters <- cbind(betak, sqrt(diag(v0)), zvalue, p)
  # colnames(parameters) <- c('estimate', paste('stderr', 1:2,
  # sep=''))
  colnames(parameters) <- c("estimate", "stderr", "z value", 
    "P(|z|>.)")
  rownames(parameters) <- paramname[-1]
  # colnames(parameters) <- c('estimate', paste('stderr', 0:1,
  # sep=''))
  rownames(parameters) <- betaname
  autocorr <- acf(res, plot = FALSE)
  
  # out <- list(MaxCosine=maxcosk,
  # parameters=parameters,varcov0=v0, varcov1=v1, varcov2=v2,
  # varcov3=v3, varcov4=v4, autocorr=autocorr,
  # iter.detail=iter.detail, call=call)
  out <- list(MaxCosine = maxcosk, parameters = parameters, 
    varcov0 = v0, varcov1 = v1, autocorr = autocorr, long = long, 
    mui = mui, y = y, iter.detail = iter.detail, call = call)
  
  class(out) <- "angular"
  out
}

#' @rdname angular
#' @method print angular
#' @param x An object, produced by the \code{\link{angular}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
print.angular <- function(x, ...) {
  cat("\nMaximum cosine:", x$MaxCosine, "\n")
  cat("\nParameters:\n")
  print.default(x$parameters, print.gap = 2, quote = FALSE, 
    right = TRUE, ...)
  cat("\n")
  invisible(x)
}


