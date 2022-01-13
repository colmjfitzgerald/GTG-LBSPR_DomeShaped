#'
#' @title 
#' Parameter estimation/optimisation model
#' 
#' @description 
#' Optimises fishing parameters by minimising NLL of length data
#' 
#' @param StockPars Stock life-history parameters in a list.
#' @param fixedFleetPars The list of fixed fleet parameter, e.g. selectivity-at-length parameters.
#' @param LenDat Observed fish counts in each length bin in an integer vector.
#' @param SizeBins A list containing ToSize and Linc (bin width) for specifying length bins.
#' @param mod Specification of simulation mode: either GTG or original LBSPR.
#' 
#' @details 
#' Calculates the initial estimates of fishing parameters for nonlinear optimisation:
#'  * \code{log(F/M) = log(0.5)}
#'  * \code{log(SL50/Linf) = log(Lmax/Linf)}
#'  * \code{log((SL95-SL50)/Linf = log(0.2*Lmax/Linf)}
#' and runs optimisation routine \code{optim} with these initial estimates. By setting \code{hessian = TRUE} 
#' we can calculate the hessian, its inverse the variance-covariance matrix of the estimates and
#' delta method approximations to standard errors.
#' 
#' Single parameter optimisation log(F/M) using the "Brent" method supports the specification of 
#' lower and upper bound constraints. Multiple parameter estimation through the "BFGS" method 
#' does not support bound constraints. 
#' 
#' Control parameters are specified through a control list and include
#' \describe{
#' \item{\code{maxit}}{Maximum number of iterations in optimisation.} 
#' \item{\code{reltol}}{Relative convergence tolerance.}
#' \item{\code{REPORT}}{Frequency of reports for the BFGS and other methods.}
#' \item{\code{trace}}{Non-negative integer 0 - no tracing, 1-6 higher levels produce more tracing information.}
#' } 
#' 
#' @return A list containing: 
#' \describe{
#' \item{\code{lbPars}}{Optimised fishing parameters.}
#' \item{\code{lbStdErrs}}{Standard errors of optimised parameters.}
#' \item{\code{fixedFleetPars}}{Fixed fishing parameters.}
#' \item{\code{PredLen}}{Predicted proportions-at-length in catch.}
#' \item{\code{NLL}}{Nonlinear log-likelihood from optimisation.}
#' \item{\code{MLE}}{Maximum likelihood estimates (log transformed values).}
#' }
#' 
#' @seealso 
#' * [optfunLBSPRDome()] calculates negative log likelihood of data given per-recruit model prediction.
#' 
#' @author K. Hommik, C. J. Fitzgerald
#' 
#' @export

optLBSPRDome <- function(StockPars, fixedFleetPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  

  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1
    SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
    SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  # control parameters
  control_opt <- list(maxit = 500, reltol = 1e-8, REPORT = 10, trace = 1)
  
    # Starting guesses
  # sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf
  # sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
 
  selectivityCurve <- fixedFleetPars$selectivityCurve
  sFM <- 0.5
  
  if(fixedFleetPars$selectivityCurve=="Logistic" && length(fixedFleetPars) == 1){ # general logistic
    # Starting guesses
    sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf 
    sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
    Start <- log(c(sFM, sSL50, sDel))  #tryFleetPars
    lowerBound <- c(-Inf, log(0.01), log(1e-5) )
    upperBound <- c(log(5), log(1+StockPars$CVLinf*StockPars$MaxSD), 0.0)
    methodOpt <- "L-BFGS-B"
    opt <- stats::optim(par = Start, fn = optfunLBSPRDome, gr = NULL, 
                 fixedFleetPars=fixedFleetPars, LenDat=LenDat, StockPars=StockPars, SizeBins=SizeBins, mod=mod, 
                 method = methodOpt, control= list(maxit=500, factr=1E9),
                 hessian = TRUE)
    mleNames <-  c("log(F/M)", "SL50/Linf", "Sdelta/Linf")
  } else{ # dome-shaped or fixed selectivity logistic
    Start <- log(c(sFM))  #tryFleetPars
    lowerBound <- -20
    upperBound <- 20
    methodOpt <- "Brent"
    opt <- stats::optim(par = Start, fn = optfunLBSPRDome, gr = NULL, 
                 fixedFleetPars=fixedFleetPars, LenDat=LenDat, StockPars=StockPars, SizeBins=SizeBins, mod=mod, 
                 method = methodOpt, lower = lowerBound, upper = upperBound, 
                 control= list(maxit=500, abstol=1E-20),
                 hessian = TRUE)
    mleNames <- c("log(F/M)")
  }
  
  # negative log-likelihood
  newNLL<-opt$value # replaces objective in nlminb

  # variance-covariance matrix for std error calculation
  varcov <- solve(opt$hessian) # inverse of hessian matrix
  
  # maximum likelihood estimators and fishing parameters
  MLE <- data.frame(Parameter = mleNames, "Initial" = Start, "Estimate" = opt$par, "Std. Error" = diag(varcov),
                    check.names = FALSE)
  
  
  # back-transform MLE to obtain fishing parameters
  newFleet <- NULL 
  newFleet$selectivityCurve <- selectivityCurve
  newFleet$FM <- exp(opt$par[1])
  lbPars <- c("F/M"  = exp(opt$par[1]))
  if(fixedFleetPars$selectivityCurve=="Logistic"){
    if(length(fixedFleetPars) == 1){
      newFleet$SL1 <- exp(opt$par[2]) * StockPars$Linf 
      newFleet$SL2 <- newFleet$SL1 + exp(opt$par[3]) * StockPars$Linf
      lbPars <- c(lbPars, 
                  "SL50" = exp(opt$par[2])*StockPars$Linf, 
                  "SL95" = (exp(opt$par[2]) + exp(opt$par[3]))*StockPars$Linf)
    } else{
      newFleet$SL1 <- fixedFleetPars$SL1
      newFleet$SL2 <- fixedFleetPars$SL2
    }
  }else if(selectivityCurve=="Knife"){
    newFleet$MLLKnife <- fixedFleetPars$MLLKnife
  } else if(selectivityCurve %in% c("Normal.sca", "Normal.loc", "logNorm")){ # prescribed values, not optimised
    newFleet$SL1 <- fixedFleetPars$SL1
    newFleet$SL2 <- fixedFleetPars$SL2
    newFleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$SLMin)) newFleet$SLMin <- fixedFleetPars$SLMin
  }
  
  # delta method to approximate standard error, CIs of estimates

  # ML estimators are log(F/M)
  sderr <- c(sqrt(exp(opt$par[1])*varcov[1,1]))
  names(sderr) = "F/M"
  if(fixedFleetPars$selectivityCurve=="Logistic" && length(fixedFleetPars) == 1){
    # log(SL50/Linf), log((SL95-SL50)/Linf) in log-space
    sderrSL50 <- sqrt((StockPars$Linf*exp(opt$par[2]))^2*varcov[2,2])
    sderrSL95 <- sqrt((StockPars$Linf^2)*exp(opt$par[2])^2*varcov[2,2] + 
                     exp(opt$par[3])^2*varcov[3,3] + 2*exp(opt$par[3])*exp(opt$par[2])*varcov[2,3])
    sderr <-  c(sderr, SL50 = sderrSL50, SL95 = sderrSL95)
  } 
  
  if (mod == "GTG") runMod <-  simLBSPRDome(StockPars, newFleet, SizeBins)
  if (mod == "LBSPR") runMod <- simLBSPRDome(StockPars, newFleet, SizeBins)

  lbPars <- c(lbPars, "SPR" = runMod$SPR)
  
  Out <- NULL 
  Out$lbPars <- lbPars      # fishing mortality, selectivity (where applicable), SPR
  Out$lbStdErrs <- sderr    # standard error for fishing mortality,... 
  Out$fixedFleetPars <- fixedFleetPars
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  Out$NLL <- newNLL
  Out$optimOut <- opt
  Out$MLE <- MLE
  return(Out)
}
