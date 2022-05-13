#'
#' @title 
#' Likelihood function for estimation model
#' 
#' @description 
#' Calculates Negative Log Likelihood of length composition data
#' 
#' @param tryFleetPars A list of fleet parameters varied as part of optimisation process.
#' @param fixedFleetPars A list of fixed fleet parameter, e.g. selectivity-at-length parameters.
#' @param LenDat Observed fish counts in each length bin in an integer vector.
#' @param StockPars Stock life-history parameters in a list.
#' @param SizeBins A list containing ToSize and Linc (bin width) for specifying length bins.
#' @param mod Specification of simulation mode: either GTG or original LBSPR.
#' 
#' @details 
#' Calculates: 
#' 
#' * Predicted probability/proportions of fish caught in fishery in each length bin 
#' for the estimated fishing parameters from the simulation routine.
#' * The negative log likelihood (multinomial) of the observed data for those fishing parameters.
#' 
#' using -sum(LenDat * log(predProb/LenProb)). Penalises SL50/Linf values close to 1 using beta distribution.
#' 
#' @return NLL 
#' @seealso 
#' * [optLBSPRDome()] optimises the likelihood (this function) of data by varying fishing parameters to obtain LBSPR fit
#' * [simLBSPRDome()] simulates per-recruit numbers-at-length in catch.
#' 
#' @author K. Hommik, C. J. Fitzgerald

optfunLBSPRDome <- function(tryFleetPars, fixedFleetPars, LenDat, StockPars, SizeBins=NULL, 
                     mod=c("GTG", "LBSPR")) {
  
  Fleet <- NULL
  Fleet$selectivityCurve <- fixedFleetPars$selectivityCurve
  if(Fleet$selectivityCurve=="Logistic"){
    if(length(tryFleetPars) == 3 & length(fixedFleetPars) == 1 & 
       c("selectivityCurve") %in% names(fixedFleetPars)){
      Fleet$SL1 <- exp(tryFleetPars[2]) * StockPars$Linf
      Fleet$SL2 <- Fleet$SL1  + (exp(tryFleetPars[3]) * StockPars$Linf)
    } else {
      Fleet$SL1 <- fixedFleetPars$SL1
      Fleet$SL2 <- fixedFleetPars$SL2
    }
  }else if(Fleet$selectivityCurve=="Knife"){
    Fleet$MLLKnife <- fixedFleetPars$MLLKnife
  }else if(Fleet$selectivityCurve %in% c("Normal.sca", "Normal.loc", "logNorm")){
    Fleet$SL1 <- fixedFleetPars$SL1 
    Fleet$SL2 <- fixedFleetPars$SL2
    Fleet$SLmesh <- fixedFleetPars$SLmesh
    if(!is.null(fixedFleetPars$SLMin)) Fleet$SLMin <- fixedFleetPars$SLMin
  } else if(Fleet$selectivityCurve == "exponentialLogistic"){
    Fleet$SL1 <- fixedFleetPars$SL1 
    Fleet$SL2 <- fixedFleetPars$SL2
    Fleet$SL3 <- fixedFleetPars$SL3
  }
  
  
  Fleet$FM <- exp(tryFleetPars[1]) # changed to 1 from 3, as other parameters are fixed

  if (mod == "GTG") runMod <-  simLBSPRDome(StockPars, Fleet, SizeBins)
  if (mod == "LBSPR") runMod <- simLBSPRDome(StockPars, Fleet, SizeBins)
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))

  #Penalty functions are used when model estimates the selectivity  
  # add penalty for SL50
  if(length(tryFleetPars) == 3 & fixedFleetPars$selectivityCurve == "Logistic"){#is.null(fixedFleetPars)){
    trySL50 <- exp(tryFleetPars[2])
    PenVal <- NLL
    Pen <- stats::dbeta(trySL50, shape1=5, shape2=0.01) * PenVal  #penalty for trySL50 values close to 1/SL50 close to Linf 
    if(!is.finite(NLL)) return(1E9 + stats::runif(1, 1E4, 1E5))
    if (Pen == 0) {Pen <- PenVal * trySL50}
    # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
    NLL <- NLL+Pen
  }

  return(NLL)
}