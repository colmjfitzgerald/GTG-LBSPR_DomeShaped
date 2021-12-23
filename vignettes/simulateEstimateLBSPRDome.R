## ----packages_preliminaries, include = FALSE----------------------------------
library(plyr) 
library(dplyr)
library(ggplot2)
library(kableExtra)

# load GTG LB-SPR routines
#source("../R/optfunLBSPRDome.R")
#source("../R/optLBSPRDome.R")
#source("../R/simLBSPRDome.R")
set.seed(seed = 9999)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----StockPars, include = TRUE------------------------------------------------
# simulated stock parameters (similar to some Irish brown trout populations) ====
StockPars <- NULL 
StockPars$MK <- 2.0
StockPars$NGTG <- 17
StockPars$Linf <- 45.0
StockPars$CVLinf <- 0.1
StockPars$MaxSD <- 2
StockPars$L50 <- 25.0  
StockPars$L95 <- 30.0 
StockPars$FecB <- 3 
StockPars$Walpha <- 0.0084
StockPars$Wbeta <- 3.1115
StockPars$Steepness <- 0.8 
StockPars$Mpow <- 0

## ----SizeBins, include = TRUE-------------------------------------------------
SizeBins <- NULL
SizeBins$Linc <- 1
SizeBins$ToSize <- StockPars$Linf*(1+StockPars$MaxSD*StockPars$CVLinf)

## ----FleetPars----------------------------------------------------------------
FleetPars <- list(selectivityCurve = "logNorm", SLmesh = 1.0)  # lognorm, single mesh size
FleetPars$SL1 <- 23.0     # k_mode from SELECT (log-normal) - mode of lognormal distribution
FleetPars$SL2 <- 0.4      # standard deviation at log scale
FleetPars$SLMin <- 0   # minimum landing limit (MLL)

## ----calc_selectivity_at_length-----------------------------------------------

# calculate selectivity at length *IF* appropriate parameters are specified
lengthFish <- seq(0, StockPars$Linf*(1 + StockPars$CVLinf*StockPars$MaxSD), length.out = 101)
if(!is.null(FleetPars$SLmesh)) meshSize <- FleetPars$SLmesh
  
if(FleetPars$selectivityCurve=="logNorm"){ 
    SLk <- FleetPars$SL1
    SLsigma <- FleetPars$SL2
    SLMin <- FleetPars$SLMin
    
    gearSelLen <- 0
    for (j in seq_along(meshSize)){
      gearSelLen <- gearSelLen + exp(-0.5*((log(lengthFish)-log((SLk)*meshSize[j]))/(SLsigma))^2)
    }
    if(!is.na(SLMin)) gearSelLen[lengthFish < SLMin] <- 0 
    gearSelLen <- gearSelLen/max(gearSelLen)
    
}

## ----plot_selectivity---------------------------------------------------------
ggplot2::ggplot() + 
  ggplot2::geom_line(data = data.frame(length = lengthFish, selectivity = gearSelLen),
            aes(x = length, y = selectivity), colour = "black", linetype = 1, size = 1.25) + 
  ggplot2::scale_x_continuous(limits = c(0, 60), breaks = seq(0,60,10)) +
  ggplot2::theme_bw()

