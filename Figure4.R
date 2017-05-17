library(ggplot2)
library(cowplot)
library(scales)
library(dplyr)

setwd("~/Downloads")

DelmeTomasetti <- read.csv("DelmeTomasetti3.csv", header=FALSE, stringsAsFactors=FALSE)

D.T.df <- as.data.frame(DelmeTomasetti)
D.T.df <- D.T.df[,-1]

DelmeTomasettiOv <- read.csv("DelmeTomasetti4.csv", header=FALSE, stringsAsFactors=FALSE)

D.T.df.Ov <- as.data.frame(DelmeTomasettiOv)
D.T.df.Ov <- D.T.df.Ov[,-1]


# fig4a -------------------------------------------------------------------


p <-
  ggplot(
    data = D.T.df.Ov,
    mapping = aes(
      x = V8,
      y = V3
    )
  )

p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 2, linetype = "solid")
)  + geom_smooth(method = "lm", mapping = aes(),colour="black",se=FALSE) +
  labs(x = "\nTotal stem cell divisions\n", y = "\nLifetime risk\n") +
  geom_point(mapping = aes(), na.rm = TRUE)

p <- p + scale_y_log10() + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))

p1 <- p

p


x <- runif(100, -3, 3) 
y <- 2 + x + rnorm(100) 

lm(log(D.T.df.Ov$V3) ~ log(D.T.df.Ov$V8)) 
lm(log(D.T.df.Ov$V3) ~ offset(0.4352*log(D.T.df.Ov$V8))) 
lm(log(D.T.df.Ov$V3[which(D.T.df.Ov$V9=='red')]) ~ offset(0.4352*log(D.T.df.Ov$V8[which(D.T.df.Ov$V9=='red')]))) 
lm(log(D.T.df.Ov$V3[which(D.T.df.Ov$V9=='blue')]) ~ offset(0.4352*log(D.T.df.Ov$V8[which(D.T.df.Ov$V9=='blue')]))) 






PolsPosRed <- cbind(c(0, 0.7e9, 6e9, 6e9),
                    c(0, 700, 700, 0))
colnames(PolsPosRed) <- c("x", "y")
PolsPosRed <- as.data.frame(PolsPosRed)


PolsPosBlue <- cbind(c(0, 0.7e9, 0.7e9, 0),
                    c(0, 700, 6e3, 6e3))
colnames(PolsPosBlue) <- c("x", "y")
PolsPosBlue <- as.data.frame(PolsPosBlue)


# Panel groups ------------------------------------------------------------


p <-
  ggplot(
    data = D.T.df,
    mapping = aes(
      x = V5,
      y = V7,
      color=V9
    )
  )

p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="none",
  axis.line = element_line(colour = "black", size = 2, linetype = "solid")
)  + labs(x = "\nTotal stem cell number (s)\n", y = "\nLifetime divisions per stem cell (d)\n")
  # geom_point(position = position_jitter(width = 1e8, height = 100))


p <- p + scale_color_manual(values = c('blue', 'red'))

p <- p + geom_polygon(data = PolsPosRed, alpha=0.4, color=NA, fill=I("blue"), mapping = aes(x = x, y = y))

p <- p + geom_polygon(data = PolsPosBlue, alpha=0.4, color=NA, fill=I("red"), mapping = aes(x = x, y = y))

p <- p + geom_point()

p <- p + stat_function(fun = function(x) {4e12/x}, colour="grey") + 
  stat_function(fun = function(x) {8e12/x}, colour="grey",size=1) + 
  stat_function(fun = function(x) {12e12/x}, colour="grey",size=1.5) + 
  stat_function(fun = function(x) {16e12/x}, colour="grey",size=2.5) + 
  
  
  xlim(0,6e9)+ylim(0,6e3)+scale_x_continuous(labels = c(0,expression(paste("2x10"^"9")),expression(paste("4x10"^"9")),expression(paste("6x10"^"9"))))

p3 <- p

p












# Panel groups risk -------------------------------------------------------






p <-
  ggplot(
    data = D.T.df.Ov,
    mapping = aes(
      x = V8,
      y = V3,
      color=V9
    )
  )

p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="none",
  axis.line = element_line(colour = "black", size = 2, linetype = "solid")
)  +
  # geom_smooth(method = "lm", mapping = aes(color=I("black"))) +
  labs(x = "\nTotal stem cell divisions\n", y = "\nLifetime risk\n") +
  geom_point()

p <- p + scale_y_log10() + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
p <- p + scale_color_manual(values = c('blue', 'red'))+ stat_smooth(data=filter(D.T.df.Ov,V9=='red'),method = "lm", mapping = aes(), se = FALSE,fullrange = TRUE)+ stat_smooth(data=filter(D.T.df.Ov,V9=='blue'),method = "lm", mapping = aes(), se = FALSE,fullrange = TRUE)

p4 <- p

p

D.T.dfL=cbind(D.T.df,1,D.T.df$V5/D.T.df$V7)

pg=plot_grid(p1,p3, p4, ncol = 3)

ggsave(filename = "Fig4.pdf", plot = pg, width = 7, height = 2, scale = 2)


# line --------------------------------------------------------------------


p <-
  ggplot(
    data = D.T.df,
    mapping = aes(
      x = D.T.df$V5/D.T.df$V7,
      y = 1,
      color=V9,label=V2
    )
  )
p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position="none",
  axis.line = element_line(colour = "black", size = 2, linetype = "solid")
)  +
  # geom_smooth(method = "lm", mapping = aes(color=I("black"))) +
  labs(x = "\nTotal stem cell divisions/Lifetime divisions per stem cell\n", y = "") +
  geom_point(
    # colour="White"
    )

p <- p + scale_y_log10() + scale_x_log10(limits=c(1e1,3e8))
p <- p + scale_color_manual(values = c('blue', 'red')) 
p <- p +  geom_text(
                                                                    x = 2.4+c(1:16,20:25)*0.23,
                                                                   y = (1:22)*0-0.03,angle=45,hjust = 1,size=3.3
                                                                  )+
  geom_segment(xend = 2.4+c(1:16,20:25)*0.23,
               yend = (1:22)*0-0.03,size=0.2)
p5 <- p

p



# Power law type model ----------------------------------------------------



D.T.df.2 <- cbind(D.T.df[,c("V8", "V3")], rep(NA, nrow(D.T.df)), rep("Data", nrow(D.T.df)))
D.T.df.3 <- cbind(D.T.df[,c("V8", "V10", "V3")], rep("Model", nrow(D.T.df)))
colnames(D.T.df.2) <- c("x", "y","lim", "Type")
colnames(D.T.df.3) <- c("x", "y","lim", "Type")

D.T.df.A <- rbind(D.T.df.2, D.T.df.3)




p <-
  ggplot(
    data = D.T.df.A,
    mapping = aes(
      x = x,
      y = y,
      color=Type
    )
  )

p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 2, linetype = "solid"),
  legend.position="none"
)  + labs(x = "\nTotal stem cell divisions\n", y = "\nLifetime risk\n")

p <- p + scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
p <- p + 
  geom_pointrange(color="green",
                  aes(
                    ymin = y,
                    ymax = lim
                  )
                  , size = 0.5, linetype = "dashed")+ 
  geom_point(size = 2.5)+
  scale_color_manual(values = c("black", adjustcolor("green", alpha.f = 0.5)))


p6 <- p

p

D.T.df$V10=1.37e-13*(D.T.df$V5)*(D.T.df$V7)^0.91

D.T.df.2 <- cbind(D.T.df[,c("V8", "V3")], rep(NA, nrow(D.T.df)), rep("Data", nrow(D.T.df)))
D.T.df.3 <- cbind(D.T.df[,c("V8", "V10", "V3")], rep("Model", nrow(D.T.df)))
colnames(D.T.df.2) <- c("x", "y","lim", "Type")
colnames(D.T.df.3) <- c("x", "y","lim", "Type")

D.T.df.A <- rbind(D.T.df.2, D.T.df.3)


p <-
  ggplot(
    data = D.T.df.A,
    mapping = aes(
      x = x,
      y = y,
      color=Type
    )
  )

p <- p + theme(
  plot.title = element_text(face="bold", size=15),
  axis.text = element_text(size=10),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size = 2, linetype = "solid"),
  legend.position="none"
)  + labs(x = "\nTotal stem cell divisions\n", y = "\nLifetime risk\n")

p <- p + scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
p <- p + 
  geom_pointrange(color="green",
                  aes(
                    ymin = y,
                    ymax = lim
                  )
                  , size = 0.5, linetype = "dashed")+ 
  geom_point(size = 2.5)+
  scale_color_manual(values = c("black", adjustcolor("green", alpha.f = 0.5)))


p7 <- p

p

