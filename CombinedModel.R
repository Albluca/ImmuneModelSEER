library(cowplot)

library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(lattice)

# Define Support Functions ------------------------------------------------




Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}

GetTick <- function(Min, Max) {
  LowTick <- 10 ^ floor(log10(Min))
  LowTick <- floor(Min / LowTick) * LowTick

  HighTick <- 10 ^ floor(log10(Max))
  HighTick <- ceiling(Max / HighTick) * HighTick

  Pwr <- floor(log10(LowTick))
  Idx <- LowTick / (10 ^ Pwr)

  Ticks <- NULL
  TickLab <- NULL

  while (Idx * 10 ^ Pwr <= HighTick) {
    Ticks <- c(Ticks, Idx * 10 ^ Pwr)

    if (Idx == 1) {
      TickLab <- c(TickLab, bquote("10" ^  ~ .(Pwr)))
    } else {
      TickLab <- c(TickLab, "")
    }

    if (Idx < 9) {
      Idx <- Idx + 1
    } else {
      Idx <- 1
      Pwr <- Pwr + 1
    }

  }

  # TickLab[1] <- LowTick
  # TickLab[length(TickLab)] <- HighTick

  AllInfo <- list()
  AllInfo$Ticks <- Ticks
  AllInfo$TickLab <- as.vector(TickLab)

  return(AllInfo)

}





k_mid <- 0.044


# Other support functions -------------------------------------------------


Fit2ParModel <- function(FittedData, a=NULL, b=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ log(a / (exp(
        b * exp(-k_mid * x)
      ) - 1)),
      data = FittedData,
      start = StartList,
      algorithm = "port",
      lower = LowVect,
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out2_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out2_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  if(Rounds > 0){
    for (i in 1:Rounds) {
      
      out2_m <-
        nls(
          log(y) ~ log(a / (exp(
            b * exp(-k_mid * x)
          ) - 1)),
          data = FittedData,
          start = StartList,
          lower = LowVect,
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      
      PEst <- summary(out2_m)$parameters[, "Estimate"]
      names(PEst) <- rownames(summary(out2_m)$parameters)
      StartList[names(PEst)] <- PEst
      
    }
  }
  
  R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA, NA)
  names(PVals) <- c("a", "b", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}



Fit3ParModel <- function(FittedData, a=NULL, b=NULL, c=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1),c = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(c), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ c*log(x)+log(a / (exp(
        b * exp(-k_mid * x)
      ) - 1)),
      data = FittedData,
      start = StartList,
      algorithm = "port",
      lower = LowVect,
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out2_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out2_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  if(Rounds > 0){
    for (i in 1:Rounds) {
      
      out2_m <-
        nls(
          log(y) ~ c*log(x)+log(a / (exp(
            b * exp(-k_mid * x)
          ) - 1)),
          data = FittedData,
          start = StartList,
          lower = LowVect,
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      
      PEst <- summary(out2_m)$parameters[, "Estimate"]
      names(PEst) <- rownames(summary(out2_m)$parameters)
      StartList[names(PEst)] <- PEst
      
    }
  }
  
  R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA, NA, NA)
  names(PVals) <- c("a", "b","c", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, c, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "c", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}



Fit3ParNegModel <- function(FittedData, a=NULL, b=NULL, c=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1),c = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(c), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out2_m <-
    nls(
      log(y) ~ c*log(x)+log(-a / (exp(
        -b * exp(-k_mid * x)
      ) - 1)),
      data = FittedData,
      start = StartList,
      algorithm = "port",
      lower = LowVect,
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out2_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out2_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  if(Rounds > 0){
    for (i in 1:Rounds) {
      
      out2_m <-
        nls(
          log(y) ~ c*log(x)+log(-a / (exp(
            -b * exp(-k_mid * x)
          ) - 1)),
          data = FittedData,
          start = StartList,
          lower = LowVect,
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      
      PEst <- summary(out2_m)$parameters[, "Estimate"]
      names(PEst) <- rownames(summary(out2_m)$parameters)
      StartList[names(PEst)] <- PEst
      
    }
  }
  
  R2Mb <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA, NA, NA)
  names(PVals) <- c("a", "b","c", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, b, c, k_mid)
  
  ReturnVect <- c(PVals, AIC(out2_m), R2Mb, deviance(out2_m), sum(summary(out2_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "c", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}








Fit1ParModel <- function(FittedData, a=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(k_mid))]
  LowVect <- rep(1e-10, times=length(StartList))
  
  out1_m <-
    nls(
      log(y) ~ log(a * exp(k_mid * x)),
      data = FittedData,
      start = StartList,
      lower = LowVect,
      algorithm = "port",
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  
  PEst <- summary(out1_m)$parameters[, "Estimate"]
  names(PEst) <- rownames(summary(out1_m)$parameters)
  StartList[names(PEst)] <- PEst
  
  for (i in 1:Rounds) {
    
    
    out1_m <-
      nls(
        log(y) ~ log(a * exp(k_mid * x)),
        data = FittedData,
        start = StartList,
        lower = LowVect,
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    
    PEst <- summary(out1_m)$parameters[, "Estimate"]
    names(PEst) <- rownames(summary(out1_m)$parameters)
    StartList[names(PEst)] <- PEst
    
  }
  
  
  R1Mb <- 1 - sum((summary(out1_m)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  
  PVals <- c(NA, NA)
  names(PVals) <- c("a", "k_mid")
  PVals[names(PEst)] <- PEst
  PVals[is.na(PVals)] <- c(a, k_mid)
  
  ReturnVect <- c(PVals, AIC(out1_m), R1Mb, deviance(out1_m), sum(summary(out1_m)$residuals^2))
  
  names(ReturnVect) <- c("a", "k", "AIC", "R^2", "Dev", "Resid")
  
  return(ReturnVect)
  
}

TryFitting <- function(Model, FittedData, a = NULL, b = NULL, c = NULL, k_mid = NULL,
                       t_0 = NULL, t_0_low = -Inf, t_0_high = Inf,
                       Rounds = 20, Retires = 10){
  
  AllFits <- NULL
  
  if(Model == "A"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit1ParModel(FittedData = FittedData, a = a, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "B"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit2ParModel(FittedData = FittedData, a = a, b = b, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  if(Model == "PLIM"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit3ParModel(FittedData = FittedData, a = a, b = b,c = c, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  if(Model == "PLIMNeg"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- Fit3ParNegModel(FittedData = FittedData, a = a, b = b,c = c, k_mid = k_mid, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "C"){
    for(i in 1:Retires){
      TR <- try(ModelPars <- FitPLModel(FittedData = FittedData, Rounds = Rounds))
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  
  if(Model == "TM"){
    for(i in 1:Retires){
      if(is.null(k_mid)){
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = NULL,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      } else{
        TR <- try(ModelPars <- FitTimeModel(FittedData = FittedData, a = a, t_0 = t_0, k = -k_mid,
                                            t_0_low = t_0_low, t_0_high = t_0_high, Rounds = Rounds))
      }
      if(length(TR)!=1){
        AllFits <- rbind(AllFits, ModelPars)
      }
    }
  }
  
  return(AllFits[which.min(AllFits[,"Resid"]),])
  
}



Predict.Model <- function(Model, Pars, NewX){
  
  if(Model == "A"){
    Model.Pred <- Pars["a"] * exp(Pars["k"] * NewX)
  }
  
  if(Model == "B"){
    Model.Pred <- Pars["a"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)
  }
  
  if(Model == "PLIM"){
    Model.Pred <- Pars["a"]*NewX^Pars["c"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)
  }
  
  if(Model == "PLIMNeg"){
    Model.Pred <- Pars["a"]*NewX^Pars["c"] / (exp(Pars["b"] * exp(-Pars["k"] * NewX)) - 1)
  }
  
  if(Model == "C"){
    Model.Pred <- exp(Pars["b"])*NewX^Pars["a"]
  }
  
  if(Model == "TM"){
    Model.Pred <- Pars["a"]/(exp(exp(Pars["k"]*(NewX-Pars["t_0"])))-1)
  }
  
  return(Model.Pred)
  
}


GetTick <- function(Min, Max) {
  LowTick <- 10 ^ floor(log10(Min))
  LowTick <- floor(Min / LowTick) * LowTick
  
  HighTick <- 10 ^ floor(log10(Max))
  HighTick <- ceiling(Max / HighTick) * HighTick
  
  Pwr <- floor(log10(LowTick))
  Idx <- LowTick / (10 ^ Pwr)
  
  Ticks <- NULL
  TickLab <- NULL
  
  while (Idx * 10 ^ Pwr <= HighTick) {
    Ticks <- c(Ticks, Idx * 10 ^ Pwr)
    
    if (Idx == 1) {
      TickLab <- c(TickLab, bquote("10" ^  ~ .(Pwr)))
    } else {
      TickLab <- c(TickLab, "")
    }
    
    if (Idx < 9) {
      Idx <- Idx + 1
    } else {
      Idx <- 1
      Pwr <- Pwr + 1
    }
    
  }
  
  # TickLab[1] <- LowTick
  # TickLab[length(TickLab)] <- HighTick
  
  AllInfo <- list()
  AllInfo$Ticks <- Ticks
  AllInfo$TickLab <- as.vector(TickLab)
  
  return(AllInfo)
  
}


Model.To.Use <- "PLIM"






# Process Infection Data --------------------------------------------------




setwd("~/Downloads")



# Process Cancer Data -----------------------------------------------------

library("SEER2R")

AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic", TXTfileName = "Incidence-Full-1Y-18.txt", UseVarLabelsInData = TRUE)

# unique(AllData$`Site_recode_ICDO3/WHO_2008`)

AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]


AllData$`Site_recode_ICDO3/WHO_2008`=trimws(AllData$`Site_recode_ICDO3/WHO_2008`)

SumData <- NULL
PlotList <- list()

PlotList <- sapply(trimws(unique(AllData$`Site_recode_ICDO3/WHO_2008`)),function(x) NULL)


k_mid <- 0.044


R2=matrix(nrow = 10,ncol = 8)


# Heatmap -----------------------------------------------------------------


for(Cancer in  "Brain"){
  
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  

    
#     
#     Tks <- GetTick(Min = max(c(0.1, min(DataSubset$Lower_Confidence_Interval, na.rm = TRUE))),
#                    Max = max(DataSubset$Upper_Confidence_Interval, na.rm = TRUE))
    
    Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
    Converted.Age[Converted.Age==84] <- NA
    
    DataSubset <- cbind(DataSubset, Converted.Age)
    
    XFil <- DataSubset$Converted.Age
    YFil <- DataSubset$AgeAdjusted_Rate
    
    DataSubset <- DataSubset[-which(is.na(XFil)),]
    
    ToFil <- YFil == 0 | XFil < 18 | is.na(XFil) | is.na(YFil)
    
    XFil <- XFil[!ToFil]
    YFil <- YFil[!ToFil]
    for(d in 1:8){
    for(c in 1:4){
      
      out2_m <-
        nls(
          log(YFil) ~ log(a)+(d-1)*log(XFil)-log(exp(
            10^(c-1) * exp(-k_mid * XFil)
          ) - 1),
          start = list(a = 1),
          algorithm = "port",
          lower = c(1e-10),
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      
      for (i in 1:10) {
        out2_m <-
          nls(
            log(YFil) ~ log(a)+(d-1)*log(XFil)-log(exp(
              10^(c-1) * exp(-k_mid * XFil)
            ) - 1),
            start = list(a = a2_m),
            lower = c(1e-10),
            algorithm = "port",
            control = nls.control(
              maxiter = 1000,
              tol = 1e-8/i,
              warnOnly = TRUE
            )
          )
        a2_m <- summary(out2_m)$parameters["a", "Estimate"]
        
      }
      #     
      #     ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
      #     ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
      #     DataSubset <- cbind(DataSubset, ModelB.Pred)
      # R2[c+6,1]=c
      R2[c+6,d] <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    }
    
    for(c in 1:5){
      
      out2_m <-
        nls(
          log(YFil) ~ (d-1)*log(XFil) +log(-a/(exp(
            -10^(c-1) * exp(-k_mid * XFil)
          ) - 1)),
          start = list(a = 1),
          algorithm = "port",
          lower = c(1e-10),
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      
      for (i in 1:10) {
        out2_m <-
          nls(
            log(YFil) ~ (d-1)*log(XFil) +log(-a/(exp(
              -10^(c-1) * exp(-k_mid * XFil)
            ) - 1)),
            start = list(a = a2_m),
            lower = c(1e-10),
            algorithm = "port",
            control = nls.control(
              maxiter = 1000,
              tol = 1e-8/i,
              warnOnly = TRUE
            )
          )
        a2_m <- summary(out2_m)$parameters["a", "Estimate"]
        
      }
      #     
      #     ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
      #     ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
      #     DataSubset <- cbind(DataSubset, ModelB.Pred)
      # R2[6-c,1]=-c
      R2[6-c,d] <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    }
    
    out2_m <-
      nls(
        log(YFil) ~ (d-1)*log(XFil) +log(a*exp(k_mid * XFil)),
        start = list(a = 1),
        algorithm = "port",
        lower = c(1e-10),
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8,
          warnOnly = TRUE
        )
      )
    a2_m <- summary(out2_m)$parameters["a", "Estimate"]
    
    for (i in 1:10) {
      out2_m <-
        nls(
          log(YFil) ~ (d-1)*log(XFil) +log(a*exp(k_mid * XFil)),
          
          start = list(a = a2_m),
          lower = c(1e-10),
          algorithm = "port",
          control = nls.control(
            maxiter = 1000,
            tol = 1e-8/i,
            warnOnly = TRUE
          )
        )
      a2_m <- summary(out2_m)$parameters["a", "Estimate"]
      
    }
    #     
    #     ModelB.Pred <- c(a2_m / (exp(b2_m * exp(-k_mid * DataSubset$Converted.Age)) - 1))
    #     ModelB.Pred[DataSubset$Converted.Age<18 & DataSubset$Converted.Age>83] <- NA
    #     DataSubset <- cbind(DataSubset, ModelB.Pred)
    # R2[6-c,1]=-c
    R2[6,d] <- 1 - sum((summary(out2_m)$residuals)^2)/sum((log(YFil)-mean(log(YFil)))^2)
    }
  }

R2p <- ifelse(R2<0,0,R2)

levelplot(R2p,col.regions=terrain.colors(100))



# Sweet spot --------------------------------------------------------------

Model.To.Use="PLIM"

  ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  
  ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
                        Rounds = 4, Retires = 100)
  
  
  Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:85)
  
  # Model.Pred[XFil<12 | XFil>83] <- NA
  
  p <-
    ggplot(
      data = DataSubset,
      mapping = aes(
        x = Converted.Age,
        y = AgeAdjusted_Rate
      )
    )+geom_line()+scale_y_log10()
  p <- p +  
    geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred),
              mapping = aes(x = X, y = Model.Pred2),
              color = adjustcolor("green", alpha.f = 0.8),
              size = 2)
  p <- p +  theme(
    plot.title=element_text(face="bold", size=25),
    axis.text=element_text(size=11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(),
    plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
  ) + labs(x = "Age \n", y = "Incidence")+
    ggtitle(Cancer) +
    geom_point(color="blue",size=0.2)
  
  p
  
  
  Model.To.Use="PLIMNeg"
  
  ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  
  ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
                        Rounds = 4, Retires = 100)
  
  
  Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:85)
  
  # Model.Pred[XFil<12 | XFil>83] <- NA
  
  p <-
    ggplot(
      data = DataSubset,
      mapping = aes(
        x = Converted.Age,
        y = AgeAdjusted_Rate
      )
    )+geom_line()+scale_y_log10()
  p <- p +  
    geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred),
              mapping = aes(x = X, y = Model.Pred2),
              color = adjustcolor("green", alpha.f = 0.8),
              size = 2)
  p <- p +  theme(
    plot.title=element_text(face="bold", size=25),
    axis.text=element_text(size=11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(),
    plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
  ) + labs(x = "Age \n", y = "Incidence")+
    ggtitle(Cancer) +
    geom_point(color="blue",size=0.2)
  
  p
  
  ModData
