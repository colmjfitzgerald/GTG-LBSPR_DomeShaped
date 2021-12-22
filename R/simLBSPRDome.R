#' @title
#' LBSPR-dome simulation model
#' 
#' @description 
#' Per recruit simulation of exploited stock
#' 
#' @details 
#' Simulates length composition, SPR, YPR of exploited stock 
#' based on inputs including:
#' \describe{
#' \item{\code{StockPars}}{Life history parameters: growth, maturity, natural mortality.}
#' \item{\code{FM}}{Relative fishing mortality F/M.}
#' \item{\code{SL50, SL95}}{Fishing gear selectivity-at-length parameters.}
#' \item{\code{SizeBins}}{Specifications for a set of length bins.}
#' } 
#' Per recruit theory is used to derive length compositions based 
#' on the specified set of length bins.
#' 
#' Based on original GTG-LBSPR code developed by A. Hordyk 
#' accessible in <https://github.com/AdrianHordyk/GTG_LBSPR>. The 
#' original code was modified to include dome-shaped selectvity.
#' Details of modification in Hommik et al. (2020) accessible at
#' <https://doi.org/10.1016/j.fishres.2020.105574> User can pre-specify 
#' selectivity parameters based on dome-shaped parameterisations from 
#' SELECT method (Millar & Fryer, 1999). Minimum landing limit (MLL) 
#' can be defined. If MLL is defined then model assumes that selection      
#' probability for lengths < MLL is 0. 
#' 
#' @param StockPars Stock life-history parameters in a list.
#' @param FleetPars A list of fishing fleet parameters: fishing mortality and selectivity-at-length.
#' @param SizeBins List specifying how to create size bins including Linc length increment.
#' 
#' @seealso 
#' * [optLBSPRDome()] optimises the likelihood of data by varying fishing parameters to obtain LBSPR fit.
#' * [optfunLBSPRDome()] calculates negative log likelihood of data given per-recruit model prediction.
#' 
#' @return Output containing SPR, YPR, catch-length, population-length compositions for fished and unfished stocks
#' @author K. Hommik, C. J. Fitzgerald
#' 

simLBSPRDome <- function(StockPars, FleetPars, SizeBins=NULL)  {

  sink(stdout(), type="message")  
  # Assign Variables 
  NGTG <- StockPars$NGTG 
  GTGLinfBy <- StockPars$GTGLinfBy 
  if (!exists("GTGLinfBy")) GTGLinfBy <- NA
  if (is.null(GTGLinfBy)) GTGLinfBy <- NA
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  MaxSD <- StockPars$MaxSD 
  MK <- StockPars$MK 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Walpha <- StockPars$Walpha 
  Wbeta <- StockPars$Wbeta 
  FecB <- StockPars$FecB 
  Steepness <- StockPars$Steepness 
  Mpow <- StockPars$Mpow
  R0 <- StockPars$R0 
  
  MLLKnife <- FleetPars$MLLKnife
  selectivityCurve <- FleetPars$selectivityCurve
  
  # extract selectivity-at-length parameters
  if(selectivityCurve == "Logistic"){
    SL50 <- FleetPars$SL1
    SL95 <- FleetPars$SL2
  } else if(selectivityCurve == "Normal.loc"){             # normal selectivity with fixed spread
    SLk <- FleetPars$SL1
    SLsigma <- FleetPars$SL2
    SLmesh <- FleetPars$SLmesh                # mesh sizes
    SLMin <- FleetPars$SLMin          # minimum landing limit, if necessary
    if (is.null(SLMin)) SLMin <- NA
  } else if(selectivityCurve == "Normal.sca"){             # normal selectivity with proportional spread
    SLk1 <- FleetPars$SL1
    SLk2 <- FleetPars$SL2
    SLmesh <- FleetPars$SLmesh 
    SLMin <- FleetPars$SLMin
    if (is.null(SLMin)) SLMin <- NA
  } else if(selectivityCurve == "logNorm"){                # lognormal selectivity
    SLk <- FleetPars$SL1
    SLsigma <- FleetPars$SL2
    SLmesh <- FleetPars$SLmesh 
    SLMin <- FleetPars$SLMin
    if (is.null(SLMin)) SLMin <- NA
  } else if(selectivityCurve=="Knife"){                    # Knife-edge selectivity
    MLLKnife <- FleetPars$MLLKnife
  }  
  FM <- FleetPars$FM 
  
  
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age # Assumed constant CV here
  
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 1
    SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  # Error Catches #
  if (!(exists("NGTG") | exists("GTGLinfBy"))) stop("NGTG or GTGLinfBy must be specified")
  if (!exists("R0")) R0 <- 1E6
  if (is.null(R0)) R0 <- 1E6
  
  # Set up Linfs for the different GTGs
  if (exists("NGTG") & !exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
    GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
  } else  if (!exists("NGTG") & exists("GTGLinfBy")) {
    DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
    NGTG <- length(DiffLinfs)
  } else if (exists("NGTG") & exists("GTGLinfBy")) {
    if (!is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, by=GTGLinfBy)
      NGTG <- length(DiffLinfs)
    } 
    if (is.na(GTGLinfBy)) {
      DiffLinfs <- seq(from=Linf-MaxSD*SDLinf, to=Linf+MaxSD*SDLinf, length=NGTG)
      GTGLinfBy <- DiffLinfs[2]-DiffLinfs[1]
    }  
  } 
  # Distribute Recruits across GTGS 
  RecProbs <- stats::dnorm(DiffLinfs, Linf, sd=SDLinf) / 
    sum(stats::dnorm(DiffLinfs, Linf, sd=SDLinf)) 
  
  # Length Bins 
  if (is.null(ToSize)) ToSize <- max(DiffLinfs, Linf + MaxSD * SDLinf)
  LenBins <- seq(from=0, by=Linc, to=ToSize)
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=(length(LenBins)-1))
  
  Weight <- Walpha * LenMids^Wbeta
  
  # Maturity and Fecundity for each GTG 
  L50GTG <- L50/Linf * DiffLinfs # Maturity at same relative size
  L95GTG <- L95/Linf * DiffLinfs # Assumes maturity age-dependant 
  DeltaGTG <- L95GTG - L50GTG
  MatLenGTG <- sapply(seq_along(DiffLinfs), function (X) 
    1.0/(1+exp(-log(19)*(LenMids-L50GTG[X])/DeltaGTG[X])))
  FecLenGTG <- MatLenGTG * LenMids^FecB # Fecundity across GTGs 
  
  if(selectivityCurve=="Logistic"){
    VulLen <- 1.0/(1+exp(-log(19)*((LenBins+0.5*Linc)-SL50)/((SL95)-(SL50)))) # Selectivity-at-Length
    
  }else if(selectivityCurve=="Normal.sca"){
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*(((LenBins+0.5*Linc)-((SLk1)*SLmesh[j]))/((SLk2)^0.5*SLmesh[j]))^2)
    }
    if(!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    VulLen <- VulLen/max(VulLen)
    
  }else if(selectivityCurve=="Normal.loc"){
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*(((LenBins+0.5*Linc)-((SLk)*SLmesh[j]))/((SLsigma)))^2)
    }
    if(!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    VulLen <- VulLen/max(VulLen)
    
  }else if(selectivityCurve=="logNorm"){
    VulLen <- 0
    for (j in seq_along(SLmesh)){
      VulLen <- VulLen + exp(-0.5*((log(LenBins+0.5*Linc)-log((SLk)*SLmesh[j]))/(SLsigma))^2)
    }
    if(!is.na(SLMin)) VulLen[LenBins < SLMin] <- 0
    VulLen <- VulLen/max(VulLen)
    
  }else if(selectivityCurve=="Knife"){    # knife-edge selectivity
    VulLen <- 0
    VulLen[(LenBins+0.5*Linc) < MLLKnife] <- 0
    VulLen[(LenBins+0.5*Linc) > MLLKnife] <- 1
    SL95 <- SL50 <- NA 
  }
  

  # Add F-mortality below MLL
  SelLen <- VulLen # Selectivity is equal to vulnerability currently
  
  # Life-History Ratios 
  MKL <- MK * (Linf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
  # Matrix of MK for each GTG
  # MKMat <- sapply(seq_along(DiffLinfs), function(X) 
  # MKL + Mslope*(DiffLinfs[X] - CentLinf))
  MKMat <- matrix(rep(MKL, NGTG), nrow=length(MKL), byrow=FALSE)
  
  FK <- FM * MK # F/K ratio 
  FKL <- FK * SelLen # F/K ratio for each length class   
  # FkL[Legal == 0] <- FkL[Legal == 0] * DiscardMortFrac 
  ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
  
  # Set Up Empty Matrices 
  # number-per-recruit at length
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) 
  
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- 
    NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), 
                                                ncol=NGTG) # number per GTG in each length class 
  # Distribute Recruits into first length class
  NPRFished[1, ] <- NPRUnfished[1, ] <- RecProbs * R0 
  for (L in 2:length(LenBins)) { # Calc number at each size class
    NPRUnfished[L, ] <- NPRUnfished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^MKMat[L-1, ]
    NPRFished[L, ] <- NPRFished[L-1, ] * ((DiffLinfs-LenBins[L])/(DiffLinfs-LenBins[L-1]))^ZKLMat[L-1, ]
    ind <- DiffLinfs  < LenBins[L]
    NPRFished[L, ind] <- 0
    NPRUnfished[L, ind] <- 0
  } 
  NPRUnfished[is.nan(NPRUnfished)] <- 0
  NPRFished[is.nan(NPRFished)] <- 0
  NPRUnfished[NPRUnfished < 0] <- 0
  NPRFished[NPRFished < 0] <- 0
  
  for (L in 1:length(LenMids)) { # integrate over time in each size class
    NatLUnFishedPop[L, ] <- (NPRUnfished[L,] - NPRUnfished[L+1,])/MKMat[L, ]
    NatLFishedPop[L, ] <- (NPRFished[L,] - NPRFished[L+1,])/ZKLMat[L, ]  
    FecGTGUnfished[L, ] <- NatLUnFishedPop[L, ] * FecLenGTG[L, ]
  }
  
  if(selectivityCurve=="Logistic"){
    VulLen2 <- 1.0/(1+exp(-log(19)*(LenMids-(SL50))/((SL95)-(SL50))))# Selectivity-at-Length
    
  }else if(selectivityCurve=="Normal.sca"){   # normal selectivity with proportional spread
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((LenMids-((SLk1)*SLmesh[j]))/((SLk2)^0.5*SLmesh[j]))^2)
    }
    if(!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(selectivityCurve=="Normal.loc"){   # normal selectivity with fixed spread
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((LenMids-((SLk)*SLmesh[j]))/(SLsigma))^2)
    }
    if(!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(selectivityCurve=="logNorm"){   # lognormal selectivity
    VulLen2 <- 0
    for (j in seq_along(SLmesh)){
      VulLen2 <- VulLen2 + exp(-0.5*((log(LenMids)-log((SLk)*SLmesh[j]))/(SLsigma))^2)
    }
    if(!is.na(SLMin)) VulLen2[LenMids < SLMin] <- 0
    VulLen2 <- VulLen2/max(VulLen2)
    
  }else if(selectivityCurve=="Knife"){   # knife-edge selectivity
    VulLen2 <- 0
    VulLen2[LenMids < MLLKnife] <- 0
    VulLen2[LenMids > MLLKnife] <- 1
    SL95 <- SL50 <- NA
  }
  
  
  
  # print(LenMids)
  # print(c(SL50, SL95))
  # plot(LenMids, VulLen2)
  
  
  # print(cbind(LenMids, VulLen2))
  NatLUnFishedCatch <- NatLUnFishedPop * VulLen2 # Unfished Vul Pop
  NatLFishedCatch <- NatLFishedPop * VulLen2 # Catch Vul Pop
  
  # plot(LenMids, apply(NatLFishedCatch, 1, sum), type="p")
  # matplot(LenMids, (NatLFishedCatch), type="l")
  
  # Expected Length Structure - standardised 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  # Calc SPR
  EPR0 <- sum(NatLUnFishedPop * FecLenGTG) # Eggs-per-recruit Unfished
  EPRf <- sum(NatLFishedPop * FecLenGTG) # Eggs-per-recruit Fished
  SPR <- EPRf/EPR0 
  
  # Equilibrium Relative Recruitment
  recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
  # RelRec/R0 - relative recruitment 
  YPR <- sum(NatLFishedPop  * Weight * VulLen2) * FM 
  Yield <- YPR * RelRec
  
  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  FitPR <- FitPR/stats::median(FitPR)
  ## Debugging
  # plot(FitPR, ylim=c(0,2)) # Should be relatively flat for equal fitness across GTG
  
  # Mslope ignored in this version 
  ObjFun <- sum((FitPR - stats::median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting Mslope 
  Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E12 # Penalty for optimising Mslope
  
  ObjFun <- ObjFun + Pen
  # print(cbind(Mslope, ObjFun, Pen))
  
  # Calculate spawning-per-recruit at each size class
  SPRatsize <- cumsum(rowSums(NatLUnFishedPop * FecLenGTG))
  SPRatsize <- SPRatsize/max(SPRatsize)
  
  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$YPR <- YPR
  Output$LCatchFished <- ExpectedLenCatchFished
  Output$LPopFished <- ExpectedLenPopFished
  Output$LCatchUnfished <- ExpectedLenCatchUnfished
  Output$LPopUnfished <- ExpectedLenPopUnfished
  Output$NatLPopFished <- NatLFishedPop
  Output$NatLPopUnFish <- NatLUnFishedPop
  Output$NatLCatchUnFish <- NatLUnFishedCatch
  Output$NatLCatchFish <- NatLFishedCatch
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- DiffLinfs[2] - DiffLinfs[1]
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$Winf <- Walpha * Linf^Wbeta
  Output$FecLen <- FecLenGTG 
  Output$MatLen <- MatLenGTG 
  Output$SelLen <- SelLen
  Output$MKL <- MKL
  Output$MKMat <- MKMat 
  Output$FKL <- FKL 
  Output$ZKLMat <- ZKLMat 
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
  Output$L50GTG <- L50GTG 
  Output$L95GTG <- L95GTG
  Output$SPRatsize <- SPRatsize
  Output$RelRec <- RelRec
  sink(type = "message")
  return(Output)
}