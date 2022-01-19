library(vegan)
library(adespatial)
library(tidyverse)
library(reshape)
library(ade4)
library(here)

##get quickMEM function (author-provided function to replace old quickPCNM function)
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:sr.value?do=export_code&codeblock=1')
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/scalog.R')

#Bird data
BD <- read.csv(here("Data/BirdData.csv"), header=TRUE)
BDTime <- read.csv(here("Data/BirdData_TimeVector.csv"))

#Phytoplankton data
PD <- read.csv(here("Data/PhytoData.csv"), header=T)
PDTime <- read.csv(here("Data/PhytoData_TimeVector.csv"), header=T)

#Invertebrate data
ID <- read.csv(here("Data/InvertData.csv"), header=T)
IDSpace <- read.csv(here("Data/InvertData_SpaceVector.csv"), header=T)

par(mar=c(1,1,1,1))
#quickMEM model for birds data
modelBD <- quickMEM(BD, BDTime, detrend = TRUE)

#quickMEM model for phyto data
modelPD <- quickMEM(PD, PDTime, detrend = TRUE)

#quickMEM model for invert data
modelID <- quickMEM(ID, IDSpace, detrend = TRUE)

#Extract lc scores for RDA 1 and 2 (significant axes) for bird data
fitted.scores.BD <- data.frame(scores(modelBD$RDA,display="lc",choices=1:2))

#Plot these RDA axes for bird data
png(here("Plots/BirdsRDAs.png"), units="in", width=12, height=7, pointsize=20, res=500)
par(mfrow=c(1,2))
plot(fitted.scores.BD$RDA1~BDTime$Year, type="l", ylab="Fitted score", xlab="Year", main="RDA 1 birds")
plot(fitted.scores.BD$RDA2~BDTime$Year, type="l", ylab="Fitted score", xlab="Year", main="RDA 2 birds")
dev.off()

#Extract lc scores for RDA 1,2,3,4 (significant axes) for phyto data
fitted.scores.PD <- data.frame(scores(modelPD$RDA,display="lc",choices=1:4))

#Plot these RDA axes for phyto data
png(here("Plots/PhytoRDAs.png"), units="in", width=12, height=12, pointsize=20, res=500)
par(mfrow=c(2,2))
plot(fitted.scores.PD$RDA1~PDTime$Year, type="l", ylab="Fitted scores", xlab="Year", main="RDA 1 phyto")
plot(fitted.scores.PD$RDA2~PDTime$Year, type="l", ylab="Fitted scores", xlab="Year", main="RDA 2 phyto")
plot(fitted.scores.PD$RDA3~PDTime$Year, type="l", ylab="Fitted scores", xlab="Year", main="RDA 3 phyto")
plot(fitted.scores.PD$RDA4~PDTime$Year, type="l", ylab="Fitted scores", xlab="Year", main="RDA 4 phyto")
dev.off()

#Extract lc scores for RDA 1 and 2 (significant axes) for invert data
fitted.scores.ID <- data.frame(scores(modelID$RDA,display="lc",choices=1:2))

#Plot these RDA axes for invert data
png(here("Plots/InvertsRDAs.png"), units="in", width=20, height=22, pointsize=20, res=500)
par(mfrow=c(1,2))
sr.value(IDSpace[,c(2,1)],fitted.scores.ID$RDA1, ylim=c(5850000,7720000), xlim=c(384000,597000), addaxes=F, possub = "topright", grid=F, csize=1.75)
sr.value(IDSpace[,c(2,1)],fitted.scores.ID$RDA2, ylim=c(5850000,7720000), xlim=c(384000,597000), addaxes=F, possub = "topright", grid=F, csize=1.75)
dev.off()

##Plot everything together with different panels
png(here("Plots/CombinedPlot.png"), units="in", width=24, height=15, pointsize=20, res=500)
par(mfrow=c(3,5))
plot.new()
text(x = 0.5, y = 0.5, paste("Birds"), 
     cex = 3, col = "black")
par(mar=c(5,5,4,2))
plot(fitted.scores.BD$RDA1~BDTime$Year, type="l", ylab="lc score RDA 1", xlab="Year", cex.axis=1.4, lwd=3, cex.lab=1.5, ylim=c(-2,2))
plot(fitted.scores.BD$RDA2~BDTime$Year, type="l", ylab="lc score RDA 2", xlab="Year", cex.axis=1.4, lwd=3, cex.lab=1.5, ylim=c(-2,2))
plot.new()

plot.new()
par(mar=c(0,0,0,0))
plot.new()
text(x = 0.5, y = 0.5, paste("Phytoplankton"), 
     cex = 3, col = "black")

par(mar=c(5,5,4,2))
plot(fitted.scores.PD$RDA1~PDTime$Year, type="l", ylab="lc score RDA 1", xlab="Year", cex.axis=1.5, lwd=3, cex.lab=1.6, ylim=c(-2,3))
plot(fitted.scores.PD$RDA2~PDTime$Year, type="l", ylab="lc score RDA 2", xlab="Year", cex.axis=1.5, lwd=3, cex.lab=1.6, ylim=c(-2,3))
plot(fitted.scores.PD$RDA3~PDTime$Year, type="l", ylab="lc score RDA 3", xlab="Year", cex.axis=1.5, lwd=3, cex.lab=1.6, ylim=c(-2,3))
plot(fitted.scores.PD$RDA4~PDTime$Year, type="l", ylab="lc score RDA 4", xlab="Year", cex.axis=1.5, lwd=3, cex.lab=1.6, ylim=c(-2,3))

par(mar=c(0,0,0,0))
plot.new()
text(x = 0.5, y = 0.5, paste("Invertebrates"), 
     cex = 3, col = "black")
par(mar=c(5,5,4,2))
sr.value(IDSpace[,c(2,1)],fitted.scores.ID$RDA1, ylim=c(5450000,7720000), xlim=c(404000,577000), addaxes=F, possub = "topright", grid=F, csize=1.75, clegend=2)
sr.value(IDSpace[,c(2,1)],fitted.scores.ID$RDA2, ylim=c(5450000,7720000), xlim=c(404000,577000), addaxes=F, possub = "topright", grid=F, csize=1.75, clegend=2)

dev.off()

#Spearman rank correlation analyses for birds
n <- rep(NA,ncol(BD))
lc.scores <- scores(modelBD$RDA, choices=1:2, scaling=3, display="lc")

d.tmpRDA1<-data.frame(species=n, RhoRDA1=n)

for (i in 1:ncol(BD)){
  test <- cor.test(lc.scores[,1], BD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA1$species[i]=names(BD[i])
  d.tmpRDA1$RhoRDA1[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA1 <- d.tmpRDA1

d.tmpRDA2<-data.frame(species=n, RhoRDA2=n)

for (i in 1:ncol(BD)){
  test <- cor.test(lc.scores[,2], BD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA2$species[i]=names(BD[i])
  d.tmpRDA2$RhoRDA2[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
  
}

StochDetRDA2 <- d.tmpRDA2
StochDetRDAsBIRDS <- cbind(StochDetRDA1,StochDetRDA2[,2])

#Spearman rank correlation analyses for phytoplankton
n <- rep(NA,ncol(PD))
lc.scores <- scores(modelPD$RDA, choices=1:4, scaling=3, display="lc")
d.tmpRDA1<-data.frame(species=n, RhoRDA1=n)

for (i in 1:ncol(PD)){
  test <- cor.test(lc.scores[,1], PD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA1$species[i]=names(PD[i])
  d.tmpRDA1$RhoRDA1[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA1 <- d.tmpRDA1

d.tmpRDA2<-data.frame(species=n, RhoRDA2=n)

for (i in 1:ncol(PD)){
  test <- cor.test(lc.scores[,2], PD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA2$species[i]=names(PD[i])
  d.tmpRDA2$RhoRDA2[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA2 <- d.tmpRDA2

d.tmpRDA3<-data.frame(species=n, RhoRDA3=n)

for (i in 1:ncol(PD)){
  test <- cor.test(lc.scores[,3], PD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA3$species[i]=names(PD[i])
  d.tmpRDA3$RhoRDA3[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA3 <- d.tmpRDA3

d.tmpRDA4<-data.frame(species=n, RhoRDA4=n)

for (i in 1:ncol(PD)){
  test <- cor.test(lc.scores[,4], PD[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA4$species[i]=names(PD[i])
  d.tmpRDA4$RhoRDA4[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA4 <- d.tmpRDA4

StochDetRDAsPHYTO <- cbind(StochDetRDA1,StochDetRDA2[,2], StochDetRDA3[,2], StochDetRDA4[,2])


#Spearman rank correlation analyses for Invertebrates
n <- rep(NA,ncol(ID))
lc.scores <- scores(modelID$RDA, choices=1:2, scaling=3, display="lc")
d.tmpRDA1<-data.frame(species=n, RhoRDA1=n)

for (i in 1:ncol(ID)){
  test <- cor.test(lc.scores[,1], ID[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA1$species[i]=names(ID[i])
  d.tmpRDA1$RhoRDA1[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA1 <- d.tmpRDA1

d.tmpRDA2<-data.frame(species=n, RhoRDA2=n)

for (i in 1:ncol(ID)){
  test <- cor.test(lc.scores[,2], ID[,i],alternative = "two.sided",
                   method = "spearman",
                   ALPHA = 0.05)
  d.tmpRDA2$species[i]=names(ID[i])
  d.tmpRDA2$RhoRDA2[i] <- if (test$p.value >= 0.05 | is.na(test$p.value)){
    print('')
  } else {
    paste(round(test$estimate,3), "*")
  }
  
}

StochDetRDA2 <- d.tmpRDA2

StochDetRDAsINVERTS <- cbind(StochDetRDA1,StochDetRDA2[,2])


