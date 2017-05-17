library(dplyr)
library(ggplot2)
library(cowplot)
library(xlsx)
library(grid)
library(SEER2R)
setwd("~/Downloads")


k_mid=0.044


# Define support functions

Convert1Year <- function(Vect) {
  as.integer(trimws(strsplit(
    x = Vect, split = "years", fixed = TRUE
  )))
}



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



FitPLModel <- function(FittedData, Rounds=10){
  
  out_pl <-
    nls(
      log(y) ~ a*log(x) + b,
      start = list(a = 1, b = 1),
      data = FittedData,
      algorithm = "port",
      control = nls.control(
        maxiter = 1000,
        tol = 1e-8,
        warnOnly = TRUE
      )
    )
  a_pl <- summary(out_pl)$parameters["a", "Estimate"]
  b_pl <- summary(out_pl)$parameters["b", "Estimate"]
  
  
  for (i in 1:Rounds) {
    out_pl <-
      nls(
        log(y) ~ a*log(x) + b,
        data = FittedData,
        start = list(a = a_pl, b = b_pl),
        algorithm = "port",
        control = nls.control(
          maxiter = 1000,
          tol = 1e-8/i,
          warnOnly = TRUE
        )
      )
    a_pl <- summary(out_pl)$parameters["a", "Estimate"]
    b_pl <- summary(out_pl)$parameters["b", "Estimate"]
    
  }
  
  RPLMb <- 1 - sum((summary(out_pl)$residuals)^2)/sum((log(FittedData$y)-mean(log(FittedData$y)))^2)
  ReturnVect <- c(a_pl, b_pl, AIC(out_pl), RPLMb, deviance(out_pl), sum(summary(out_pl)$residuals^2))
  
  names(ReturnVect) <- c("a", "b", "AIC", "R^2", "Dev", "Resid")
  
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


AllData <- read.SeerStat(DICfileName = "Incidence-Full-1Y-18.dic",
                         TXTfileName = "Incidence-Full-1Y-18.txt",
                         UseVarLabelsInData = TRUE)
AllData$AgeAdjusted_Rate <- as.numeric(levels(AllData$AgeAdjusted_Rate))[AllData$AgeAdjusted_Rate]
AllData$Lower_Confidence_Interval <- as.numeric(levels(AllData$Lower_Confidence_Interval))[AllData$Lower_Confidence_Interval]
AllData$Upper_Confidence_Interval <- as.numeric(levels(AllData$Upper_Confidence_Interval))[AllData$Upper_Confidence_Interval]



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


SumData <- NULL


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
  
  
    ModData.A <- TryFitting(Model = "A", FittedData = FittedData,k_mid = k_mid,
                            Rounds = 4, Retires = 50)
    
    ModData.B <- TryFitting(Model = "B", FittedData = FittedData,k_mid = k_mid,
                            Rounds = 9, Retires = 100)
    
    
    
    names(ModData.A) <- paste("IM-I",names(ModData.A), sep = ":")
    names(ModData.B) <- paste("IM-II",names(ModData.B), sep = ":")
    SumData <- rbind(SumData, c(Disease, ModData.A, ModData.B))
    
    
}


SumData[,1] <- trimws(SumData[,1])
colnames(SumData)[1] <- "Disease"

write.csv(x = SumData, file = "SuppT1Diseases.csv", row.names = FALSE)


# MalevsFemale Cancers ----------------------------------------------------

Model.To.Use <- "A"
SumData <- NULL

for(Cancer in unique(AllData$`Site_recode_ICDO3/WHO_2008`)){
  
  if(Cancer == "All Sites"){
    next()
  }
  
   print(Cancer)
   
   DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                           AllData$Race_recode_W_B_AI_API == "White" &
                           AllData$Sex == "  Male",]
   
   if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
     
     SumData <- rbind(SumData, c(Cancer, NA, NA))
     next
   }
   
   Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
   Converted.Age[Converted.Age==84] <- NA
   
   DataSubset <- cbind(DataSubset, Converted.Age)
   
   XFil <- DataSubset$Converted.Age
   YFil <- DataSubset$AgeAdjusted_Rate
   
   #    Tks <-
   #      GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
   #              Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
   
   ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
   
   XFil <- XFil[!ToFil]
   YFil <- YFil[!ToFil]
   
   FittedData <- data.frame(cbind(XFil, YFil))
   colnames(FittedData) <- c("x", "y")
   
   ModData.M <- TryFitting(Model = Model.To.Use, FittedData = FittedData,
                           Rounds = 4, Retires = 50)
   
   
   DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                           AllData$Race_recode_W_B_AI_API == "White" &
                           AllData$Sex == "  Female",]
   
   if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
     
     SumData <- rbind(SumData, c(Cancer, NA, NA))
     next
   }
   
   Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
   Converted.Age[Converted.Age==84] <- NA
   
   DataSubset <- cbind(DataSubset, Converted.Age)
   
   XFil <- DataSubset$Converted.Age
   YFil <- DataSubset$AgeAdjusted_Rate
   
   #    Tks <-
   #      GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
   #              Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
   
   ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
   
   XFil <- XFil[!ToFil]
   YFil <- YFil[!ToFil]
   
   FittedData <- data.frame(cbind(XFil, YFil))
   colnames(FittedData) <- c("x", "y")
   
   ModData.F <- TryFitting(Model = Model.To.Use, FittedData = FittedData,
                           Rounds = 4, Retires = 50)
  
   MaleAlpha=ModData.M[2]
   FemaleAlpha=ModData.F[2]
   
#    
#    names(ModData.F) <- paste(names(ModData.F), "F", sep = "_")
#    names(ModData.M) <- paste(names(ModData.M), "M", sep = "_")
   SumData <- rbind(SumData, c(Cancer, FemaleAlpha, MaleAlpha))
   
   
}


SumData[,1] <- trimws(SumData[,1])
colnames(SumData) <- c("Cancer/Site","Female alpha","Male alpha")

write.csv(x = SumData, file = "MaleVSFemale3.csv", row.names = FALSE)

barplot(table(SumData[,2] < SumData[,3])[c("TRUE", "FALSE")], names.arg = c("Male larger k", "Female larger k"),
        main = "Cancer with R^2 > .85")


# Make Table ----------------------------------------------------

Model.To.Use <- "A"
SumData <- NULL

for(Cancer in unique(AllData$`Site_recode_ICDO3/WHO_2008`)){
  
  if(Cancer == "All Sites"){
    next()
  }
  
  print(Cancer)
  
  DataSubset <- AllData[AllData$`Site_recode_ICDO3/WHO_2008` == Cancer &
                          AllData$Race_recode_W_B_AI_API == "White" &
                          AllData$Sex == "Male and female",]
  
  if(sum(DataSubset$AgeAdjusted_Rate, na.rm = TRUE)==0){
    next
  }
  
  Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
  Converted.Age[Converted.Age==84] <- NA
  
  DataSubset <- cbind(DataSubset, Converted.Age)
  
  XFil <- DataSubset$Converted.Age
  YFil <- DataSubset$AgeAdjusted_Rate
  
  #    Tks <-
  #      GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
  #              Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))
  
  ToFil <- YFil == 0 | XFil < 16 | is.na(XFil) | is.na(YFil)
  
  XFil <- XFil[!ToFil]
  YFil <- YFil[!ToFil]
  
  FittedData <- data.frame(cbind(XFil, YFil))
  colnames(FittedData) <- c("x", "y")
  
  ModData.A <- TryFitting(Model = "A", FittedData = FittedData,k_mid = k_mid,
                          Rounds = 4, Retires = 50)
  
  ModData.B <- TryFitting(Model = "B", FittedData = FittedData,k_mid = k_mid,
                          Rounds = 9, Retires = 100)
  ModData.C <- TryFitting(Model = "C", FittedData = FittedData,
                          Rounds = 9, Retires = 100)
  
  
  
  names(ModData.A) <- paste("IM-I",names(ModData.A), sep = ":")
  names(ModData.B) <- paste("IM-II",names(ModData.B), sep = ":")
  names(ModData.C) <- paste("PLM",names(ModData.C), sep = ":")
  SumData <- rbind(SumData, c(Cancer, ModData.A, ModData.B, ModData.C))
  
  
}


SumData[,1] <- trimws(SumData[,1])
colnames(SumData)[1] <- "Cancer/Site"

write.csv(x = SumData, file = "SuppT1temp.csv", row.names = FALSE)


# T-lymphoblastic leukemia ------------------------------------------------

tlym=read.csv("lymphoblastic.csv", header=TRUE, stringsAsFactors=FALSE)


DataSubset <- tlym[ tlym$Sex == "Male and female",]



Converted.Age <- Convert1Year(DataSubset$Age_recode_with_single_ages_and_85)
Converted.Age[Converted.Age==84] <- NA

DataSubset <- cbind(DataSubset, Converted.Age)  
head(DataSubset)

p <-
  ggplot(
    data = DataSubset,
    mapping = aes(
      x = Converted.Age,
      y = Count/Population*100000
    )
  )
p <- p +  
  geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred.C),
            mapping = aes(x = X, y = Model.Pred2),
            color = adjustcolor("green", alpha.f = 0.8),
            size = 2)+
  geom_line(data=data.frame(X=18:85,Model.Pred=Model.Pred.B),
            mapping = aes(x = X, y = Model.Pred),
            color = adjustcolor("red", alpha.f = 1),
            size = 2)+
  geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred.A),
            mapping = aes(x = X, y = Model.Pred2),
            color = adjustcolor("red", alpha.f = 0.3),
            size = 2)
p <- p +  theme(
  plot.title=element_text(face="bold", size=25),
  axis.text=element_text(size=11),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
) + labs(x = "Age \n", y = "Incidence")+
  ggtitle("T-Lymphoblastic Leukemia") +
  geom_point(color="blue",size=0.2)

p <- p + 
#   scale_y_log10(
#     #     breaks = Tks$Ticks,
#     #     labels = Tks$TickLab,
#     #     limits = range(Tks$Ticks)
#   ) + 
  scale_x_continuous(
    breaks = c(0, 45, 90),
    limits = c(0, 90)
  )

p <- p  + geom_pointrange(color="blue",
                          aes(
                            ymin = Lower_Confidence_Interval,
                            ymax = Upper_Confidence_Interval
                          ),size=0.2
)

p



# TB Prevalence -----------------------------------------------------------


Dis1 <- read.xlsx("Diseases.xlsx", sheetIndex = 4,stringsAsFactors=FALSE)


XFil <- Dis1$Age
YFil <- Dis1$TB.Prevalence

   Tks <-
     GetTick(Min = (10 ^ -1.2) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)),
             Max = (10 ^ 1.9) * sqrt(min(YFil, na.rm = TRUE)*max(YFil, na.rm = TRUE)))



FittedData <- data.frame(cbind(XFil, YFil))
colnames(FittedData) <- c("x", "y")
  
ModData.A <- TryFitting(Model = "A", FittedData = FittedData,k_mid = k_mid,
                        Rounds = 4, Retires = 50)

Model.Pred.A <- Predict.Model(Model = "A", Pars = ModData.A, NewX = 18:70)  

ModData.B <- TryFitting(Model = "B", FittedData = FittedData,k_mid = k_mid,
                                                                                                  Rounds = 4, Retires = 50)

Model.Pred.B <- Predict.Model(Model = "B", Pars = ModData.B, NewX = 18:70)



p <-
  ggplot(
    data = Dis1,
    mapping = aes(
      x = Age,
      y = TB.Prevalence
    )
  )
p <- p +  
#   geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred.C),
#             mapping = aes(x = X, y = Model.Pred2),
#             color = adjustcolor("green", alpha.f = 0.8),
#             size = 2)+
  geom_line(data=data.frame(X=18:70,Model.Pred.B=Model.Pred.B),
            mapping = aes(x = X, y = Model.Pred.B,color="Immunological Model II"),
            size = 2,show.legend = TRUE)+
  geom_line(data=data.frame(X=18:70,Model.Pred.A=Model.Pred.A),
            mapping = aes(x = X, y = Model.Pred.A,color="Immunological Model I"),
            size = 2,show.legend = TRUE) + scale_color_manual(values = c(adjustcolor("red", alpha.f = 0.3), adjustcolor("red", alpha.f = 1)))

p <- p +  theme(
  plot.title=element_text(face="bold", size=25),
  axis.text=element_text(size=11),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
) + labs(x = "Age \n", y = "Prevalence")+
  ggtitle("Tuberculosis in Cambodia") +
  geom_point(color="blue",size=2)

p <- p + 
    scale_y_log10(
          breaks = Tks$Ticks,
          labels = Tks$TickLab,
          limits = c(100,5000)
    ) + 
  scale_x_continuous(
    breaks = c(0, 45, 90),
    limits = c(0, 90)
  )


p


# TREC plot ---------------------------------------------------------------


Dis1 <- read.xlsx("Healthies.xls", sheetIndex = 1,stringsAsFactors=FALSE)

head(Dis1)
p <-
  ggplot(
    data = Dis1,
    mapping = aes(
      x = Age,
      y = CD4.sjTREC
    )
  )
# p <- p +  
#   #   geom_line(data=data.frame(X=18:85,Model.Pred2=Model.Pred.C),
#   #             mapping = aes(x = X, y = Model.Pred2),
#   #             color = adjustcolor("green", alpha.f = 0.8),
#   #             size = 2)+
#   geom_line(data=data.frame(X=18:70,Model.Pred=Model.Pred.B),
#             mapping = aes(x = X, y = Model.Pred),
#             color = adjustcolor("red", alpha.f = 1),
#             size = 2)+
#   geom_line(data=data.frame(X=18:70,Model.Pred=Model.Pred.A),
#             mapping = aes(x = X, y = Model.Pred),
#             color = adjustcolor("red", alpha.f = 0.3),
#             size = 2)
p <- p +   theme_bw() +theme(
  plot.title=element_text(face="bold", size=25),
  axis.text=element_text(size=11),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm"),
  panel.border = element_rect(colour = "black", fill=NA, size=1)
) + labs(x = "Age \n", y = "TRECs per μg of DNA")+
  ggtitle("T-cell production byproducts") +
  geom_point(color="black",size=3,shape=1)+geom_smooth(method = "lm",se=FALSE,linetype=2,color="black")+
  geom_point(data = Dis1,
             mapping = aes(
               x = Age,
               y = CD8.sjTREC
             ),color="black",size=3,shape=2)

p <- p + 
  scale_y_log10(
#     breaks = Tks$Ticks,
#     labels = Tks$TickLab,
#     limits = c(100,5000)
  ) + 
  scale_x_continuous(
    breaks = c(0, 45, 90),
    limits = c(0, 90)
  )


p




# Survivability Pivot Age inverse correlation -----------------------------


Dis1 <- read.xlsx("SuppT1.xlsx", sheetIndex = 1,stringsAsFactors=FALSE)

Dis2=Dis1[!Dis1$Pivot.Age==0,]

Dis3=Dis1[Dis1$IM.II.b>-400,]



p <-
  ggplot(
    data = Dis2,
    mapping = aes(
      x = Pivot.Age,
      y = X5.yr.survival
    )
  )

p <- p +  theme(
  plot.title=element_text(face="bold", size=25),
  axis.text=element_text(size=11),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
) + labs(x = "Pivot Age \n", y = "5 year survival (%)")+
  geom_point(color="blue",size=2)




p


p <-
  ggplot(
    data = Dis3,
    mapping = aes(
      x = IM.II.b,
      y = X5.yr.survival
    )
  )

p <- p +  theme(
  plot.title=element_text(face="bold", size=25),
  axis.text=element_text(size=11),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(),
  plot.margin = unit(c(0.2, 0.1, 0.1, 0.1), "cm")
) + labs(x = "Immune model II parameter B \n", y = "5 year survival (%)")+
  geom_point(color="blue",size=2)




p

cor.test(Dis2$Pivot.Age,Dis2$X5.yr.survival)
cor.test(Dis3$IM.II.b,Dis3$X5.yr.survival)


