library(dplyr)
library(ggplot2)
library(cowplot)
library(xlsx)
library(grid)

setwd("~/Downloads")


k_mid=0.044

Fit2ParModel <- function(FittedData, a=NULL, b=NULL, k_mid=NULL, Rounds=10){
  
  StartList <- list(a = runif(1), b = runif(1), k_mid = 0.05)
  
  StartList <- StartList[c(is.null(a), is.null(b), is.null(k_mid))]
  LowVect <- rep(-Inf, times=length(StartList))
  
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

TryFitting <- function(Model, FittedData, a = NULL, b = NULL, k_mid = NULL,
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


Model.To.Use <- "B"




Dis1 <- read.xlsx("Diseases2.xlsx", sheetIndex = 1,stringsAsFactors=FALSE)

Dis1$Age=as.numeric(as.character(Dis1$Age))
Dis1$Inc=as.numeric(as.character(Dis1$Inc))
Dis1$LowerCI=as.numeric(as.character(Dis1$LowerCI))
Dis1$UpperCI=as.numeric(as.character(Dis1$UpperCI))
Dis1$Pop=as.numeric(as.character(Dis1$Pop))


PlotList <- list()
orderlist = list()
PlotList <- sapply(trimws(unique(Dis1$Name)[!unique(Dis1$Name)%in%c(NA,"West Nile Virus Disease Male","West Nile Virus Disease Female")]
),function(x) NULL)


for(Disease in na.omit(unique(Dis1$Name))){

  if(Disease %in% c("West Nile Virus Disease Male","West Nile Virus Disease Female")){
    next()
  }
  
  print(Disease)
  
  
  XFil <- filter(Dis1,Name==Disease)$Age
  YFil <- filter(Dis1,Name==Disease)$Inc
  
  Tks <-
    GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
            Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))

  ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  
  ModData <- TryFitting(Model = Model.To.Use, FittedData = FittedData, k_mid = k_mid,
                        Rounds = 4, Retires = 50)
  
  ModData2 <- TryFitting(Model = "A", FittedData = FittedData, k_mid = k_mid,
                                                                        Rounds = 4, Retires = 50)
  
  Model.Pred <- Predict.Model(Model = Model.To.Use, Pars = ModData, NewX = 18:92)
  
  Model.Pred2 <- Predict.Model(Model = "A", Pars = ModData2, NewX = 18:92)
  # Model.Pred[XFil<12 | XFil>83] <- NA
  
  orderlist=rbind(orderlist,data_frame(ModData["R^2"]))
  
  df=mutate(filter(Dis1,Name==Disease),LowerCI=replace(LowerCI,LowerCI<Tks$Ticks[1],Tks$Ticks[1]))

  
  p <- ggplot(data = df, aes(x = Age, y = Inc))
  p <- p +geom_line(color="blue")+ theme(
    plot.title=element_text(face="bold", size=25),
    axis.text=element_text(size=25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = if(Disease%in% c("West Nile Virus Disease","Influenza A")){
    element_rect(fill = adjustcolor("green", alpha.f = 0.2))
  }  else{element_rect(fill = adjustcolor("yellow", alpha.f = 0.2))}
  ,
    plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm"),
    plot.title = element_text(face = "bold")
  )+
    scale_y_log10(
      breaks = Tks$Ticks,
      labels = Tks$TickLab,
      limits = c(Tks$Ticks[1], Tks$Ticks[27])
    ) +
   
    geom_line(data=data.frame(X=18:92,Model.Pred=Model.Pred),
              mapping = aes(x = X, y = Model.Pred),
              color = adjustcolor("red", alpha.f = 1),
              size = 2)+
    geom_line(data=data.frame(X=18:92,Model.Pred2=Model.Pred2),
              mapping = aes(x = X, y = Model.Pred2),
              color = adjustcolor("orange", alpha.f = 1),
              size = 2)+
    scale_x_continuous(breaks = c(0, 45, 90), limits = c(0, 92)) +
    labs(x = NULL, y = NULL) +
    geom_linerange(color="blue",
                   aes(
                     ymin = LowerCI,
                     ymax = UpperCI
                   )
                   , size = 1.5)+
    geom_point(na.rm = TRUE, size = 2, color="blue") +
    ggtitle(Disease)
  
  p
  
  # PlotList[[length(PlotList) + 1]] <- p1 + p2 + p3
  
  
    PlotList[[trimws(Disease)]] <- p

  
}


# quartz()

ml=plot_grid(plotlist=PlotList[order(orderlist$`ModData["R^2"]`,decreasing = TRUE)], ncol = 3)

# for(i in 1:4){
# print(
#   prop.test(filter(Dis1,Name==Disease)$Inc[i]*filter(Dis1,Name==Disease)$Pop[i]*1e-5,filter(Dis1,Name==Disease)$Pop[i],conf.level=0.9999)
#       $conf.int[2])
# }

ggsave(filename = "Infection.V9.ggplot.pdf", plot = ml, width = 8, height = 8, scale = 2)

