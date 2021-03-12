#load libraries

library(vegan)
library(BiodiversityR)
library(MASS)
library(dplyr)
library(vegan3d)
library(ggplot2)
library(here)
library(FD)
library(nlme)
library(lmerTest)
library(car)
library(vegan)
library(lattice)
library(plyr)
library(gridExtra)
library(RColorBrewer)
library(Rmisc)
library(ggplot2)
library(mgcv)
library(MuMIn)
library(lsmeans)

#load csv file (with estimations) with chemico-physical and macroinvertebrate data
#NAs in original measurements estimated as mean values from data collected before and after
#Water Temperature NAs within date for 1 or more ponds estimated from WT in other ponds for that date

Data=read.table("ENV_MACRO_02_03_21.csv",h=T, sep=",")  # sites, abiotic, biotic data
names(Data)
str(Data)
colnames(Data)[which(names(Data) == "Algae")] <- "SubVeg"

str(Data)

###adjust date format
Data$newdate <- strptime(as.character(Data$Date), "%d/%m/%Y")
Data$newdate<-format(Data$newdate, format="%Y-%m-%d")
Data$newdate = as.Date(Data$newdate,"%Y-%m-%d")

#order by date
Data<-Data[ order(Data$newdate),]

#convert in factor (from character)
Data$sampler_ID<-as.factor(Data$sampler_ID)
Data$sampler_ID2<-as.factor(Data$sampler_ID2)
Data$SeasonUnique<-as.factor(Data$SeasonUnique)
Data$Year<-as.factor(Data$Year)
Data$Date2<-as.factor(Data$Date2)
Data$Date3<-as.factor(Data$Date3)
Data$Season<-as.factor(Data$Season)
Data$Pond<-as.factor(Data$Pond)
Data$Treat<-as.factor(Data$Treat)

#re-order ponds from 1 to 12
levels(Data$Pond)
Data$Pond <- factor(Data$Pond, levels = c("P1", "P2", "P3", "P4",
                                          "P5", "P6", "P7", "P8", 
                                          "P9", "P10", "P11", "P12"))
levels(Data$Pond)

#re-order seasons over the years
Data$Season_Year <- ordered(Data$Season_Year, levels = c("W_17", "SP_18", "SU_18", "A_18", "W_18", "SP_19", "SU_19"))
Data$Season_Year2 <- ordered(Data$Season_Year2, levels = c("Win1", "Sp1", "Su1", "Au1", "Win2", "Sp2", "Su2"))
levels(Data$Season_Year)
levels(Data$Season_Year2)

###removed variables with too many NAs (not representative) and variables not related to water-environment
Datared<-Data[, !(colnames(Data) %in% c("AbuZoo","TOC","LeafDec", "LeafDec_mic","LeafDec_shr",
                                        "GreenTea","RedTea","Anura_Bufonidae","Anura_Ranidae","Urodela_Salamandridae",
                                        "NI", "Coleoptera_Curculionidae_adult", "Lepidoptera_Crambidae", "Diptera_Chironomidae_pupae"))]

####AIM 1: analysis env diversity (Time points ca. on a Month scale)
names(Datared)
SubFactorsEnv<- Datared[,c(23:40)]
row.names(Datared)



DataAggEnv=as.data.frame(aggregate(SubFactorsEnv, by=list(Datared$Month, 
                                                          Datared$Year, 
                                                          Datared$Pond,
                                                          Datared$Distance_lake,
                                                          Datared$Position,
                                                          Datared$MonthFromBeg,
                                                          Datared$Season,
                                                          Datared$Season_Year,
                                                          Datared$Season_Year2,
                                                          Datared$Treat), FUN=mean,na.rm = T))
head(DataAggEnv)
library(data.table)
head(DataAggEnv)
setnames(DataAggEnv, old = c(1:10), new = c('Month','Year','Pond','Dist_Lake','Position',
                                            'MonthFromBeg',"Season",'Season_Year','Season_Year2','Treat'))


#check variables
op<-par(mfrow=c(4,5), mar=c(3,3,3,0.5))
boxplot(DataAggEnv$Cond~as.factor(DataAggEnv$MonthFromBeg), main="Cond")
boxplot(DataAggEnv$pH~as.factor(DataAggEnv$MonthFromBeg), main="pH")
boxplot(DataAggEnv$O2~as.factor(DataAggEnv$MonthFromBeg), main="O2")
boxplot(DataAggEnv$WT~as.factor(DataAggEnv$MonthFromBeg), main="WT")
boxplot(DataAggEnv$DOC~as.factor(DataAggEnv$MonthFromBeg), main="DOC")
boxplot(DataAggEnv$FL~as.factor(DataAggEnv$MonthFromBeg), main="FL")
boxplot(DataAggEnv$CL~as.factor(DataAggEnv$MonthFromBeg), main="Cl")
boxplot(DataAggEnv$NO3~as.factor(DataAggEnv$MonthFromBeg), main="NO3")
boxplot(DataAggEnv$PO4~as.factor(DataAggEnv$MonthFromBeg), main="PO4")
boxplot(DataAggEnv$SO4~as.factor(DataAggEnv$MonthFromBeg), main="SO4")
boxplot(DataAggEnv$Na~as.factor(DataAggEnv$MonthFromBeg), main="Na")
boxplot(DataAggEnv$NH4~as.factor(DataAggEnv$MonthFromBeg), main="NH4")
boxplot(DataAggEnv$K~as.factor(DataAggEnv$MonthFromBeg), main="K")
boxplot(DataAggEnv$Mg~as.factor(DataAggEnv$MonthFromBeg), main="Mg")
boxplot(DataAggEnv$Ca~as.factor(DataAggEnv$MonthFromBeg), main="Ca")
boxplot(DataAggEnv$WatLev~as.factor(DataAggEnv$MonthFromBeg), main="WatLev")
boxplot(DataAggEnv$SubVeg~as.factor(DataAggEnv$MonthFromBeg), main="SubVeg")


dev.off()

###check outliers
names(DataAggEnv)

names(EnvAgg)

op<-par(mfrow=c(4,5), mar=c(3,3,3,0.5))
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$Cond, main="Cond")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$pH, main="pH")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$O2, main="O2")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$WT, main="WT")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$DOC, main="DOC")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$FL, main="FL")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$CL, main="Cl")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$NO3, main="NO3")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$PO4, main="PO4")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$SO4, main="SO4")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$Na, main="Na")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$NH4, main="NH4")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$K, main="K")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$Mg, main="Mg")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$Ca, main="Ca")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$WatLev, main="WatLev")
plot(y=DataAggEnv$MonthFromBeg,x= DataAggEnv$SubVeg, main="SubVeg")


dev.off()

###Plots over time
names(DataAggEnv)
EnvAggTime<- DataAggEnv


library(tidyr)

names(EnvAggTime)
long_data<-EnvAggTime%>%
  pivot_longer(c(11:26, 28), names_to = "variables", values_to = "value")

long_data$variables<-as.factor(long_data$variables)
levels(long_data$variables)       
long_data$variables<-factor(long_data$variables, levels = c("Cond","pH","O2","WT","DOC","FL", "Cl","NO3","PO4",     
                                                            "SO4","Na","NH4","K","Mg","Ca","WatLev","SubVeg"))

tiff("All_Treat.tif", res=600, compression = "lzw", height=13, width=20, units="in")

ggplot(long_data,aes(x=MonthFromBeg, y = value)) +
  geom_point() +
  geom_smooth(aes())+
  facet_wrap(~variables, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 10,15, 20))+
  labs(x = "Months from Beg", y = "response")

dev.off()

tiff("All_Ponds.tif", res=600, compression = "lzw", height=13, width=20, units="in")

ggplot(long_data,aes(x=MonthFromBeg, y = value)) +
  geom_point(alpha=0.2) +
  geom_smooth(aes(col=Pond), size=0.5,se=F)+
  #geom_smooth(fill="red", col="red", linetype = "solid", size=1)+
  facet_wrap(.~variables, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 10,15, 20))+
  labs(x = "Months from Beg", y = "response")+
  theme(strip.text.x = element_text(size = 11, 
                                    color = "black", face = "bold.italic"))
dev.off()



##select environmental parameters  to assess correlations and VIF
names(DataAggEnv)
EnvAgg<- DataAggEnv[,c(11:26, 28)]

names(EnvAgg)
str(EnvAgg)
library(psych)

#pairs.panels(ENVAgg, 
#             method = "spearman", # correlation method
#             hist.col = "#00AFBB",
#             density = TRUE,  # show density plots
 #            ellipses = F,
#             cex.labels=0.8# show correlation ellipses
#)


library(gclus)
##remove rows for which 1 parameter contain NAs! These will be the dataset I will run the analysis later
DataAggEnv_NA<-DataAggEnv[complete.cases(DataAggEnv), ]
EnvAgg_NA<-EnvAgg[complete.cases(EnvAgg), ]
##log bse 10 of env dataset

LogEnvAgg_NA<-log10(EnvAgg_NA+1)

cor.res <- cor(LogEnvAgg_NA,method = "spearman")
round(cor.res, 2)
env.order<-order.single(cor.res)
op<-par(mfrow=c(1,1),oma=c(0,0,0,0))

panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt), cex=1.2)
  text(.5, .75, Signif, cex=1.2)
}

panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


pairs(LogEnvAgg_NA[, env.order], lower.panel=panel.smooth, 
      upper.panel=panel.cor, diag.panel=panel.hist, 
      main="Spearman", cex.labels=1.2, cex=1.2)

par(op)


library(psycho)

cor.res <- cor(LogEnvAgg_NA,method = "spearman") #rank correlation
library("corrplot")
corrplot.mixed(cor.res, tl.cex = 0.8, tl.col = "black", addCoefasPercent = TRUE)

names(LogEnvAgg_NA)
#ENVAgg_red <- subset(ENVAgg_NA, select = -c(3,5,13))###removed O2, DOC, K


##VIF


#VIF FUNCTION (Zuur et al.)
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  #cat("Correlations of the variables\n\n")
  #tmp_cor <- cor(dataz,use="complete.obs")
  #print(tmp_cor)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. 
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}


corvif(LogEnvAgg_NA)

########################MULTIVARIATE ANALYSIS


##perform DCA to see whether to use linear or uni-modal ordination methods
Env.dca <- decorana(LogEnvAgg_NA)
Env.dca
#Axis length DCA1<3 -> use of Linear approaches PCA-RDA

############################

#PCA

library(ggfortify)
library(factoextra)
library(FactoMineR)

###estimate missing data based on relation between individuals and variable correlation

#library(missMDA)
#nb<-estim_ncpPCA(LogENVAgg, scale=TRUE)
#comp<-imputePCA(LogENVAgg, ncp=2, scale=TRUE)
#res.pca<-PCA(comp$completeObs)
#mi<-MIPCA(LogENVAgg, scale=TRUE, ncp=2)
#op<-par(mfrow=c(2,2), mar=c(3,3,3,1))
#plot(mi)

######run pca and find correlations Var with axis
ENVAgg_NA_cent <- scale(LogEnvAgg_NA,scale=T,center=T)
pca_res <- PCA(ENVAgg_NA_cent,  graph=T)
#plot(pca_res,choix="ind",habillage="Na")

##find correlation and p-value for varaible-axisPCA
dimdesc(pca_res, axes=1:2)

###pca with vegan - broken stick approach for selection of PCA axes to consider in the analysis
env.pca<-rda(LogEnvAgg_NA,scale=T, center=T)
biplot(env.pca)
summary(env.pca)
ev<-env.pca$CA$eig
ev[ev>mean(ev)]
n<-length(ev)
bsm<-data.frame(j=seq(1:n), p=0)
bsm$p[1]<- 1/n

for (i in 2:n){
  bsm$p[i] = bsm$p[i-1]+ (1/(n+1-i))  
}

bsm$p<-100*bsm$p/n
bsm
barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, las=2, ylim=c(0,35))


####Other exploration with prcomp package
pca_res2 <- prcomp(LogEnvAgg_NA, scale.=TRUE, center=TRUE)
names(pca_res2)

######pca score coordinates and Variance PCA scree plot

summary(pca_res2)
var <- get_pca_var(pca_res2)
fviz_eig(pca_res2, addlabels = TRUE, ylim = c(0, 40))+theme_minimal()

#DataAggTot_NA$Treat_Year<-do.call(paste,c(DataAggTot_NA[c("Treat","Year")],sep="-"))
#DataAggTot_NA$Cat_Y<-ifelse(DataAggTot_NA$newdate<"2018-11-16","1stYear","2ndYear")
#DataAggTot_NA$Treat_Cat_Y<-do.call(paste,c(DataAggTot_NA[c("Treat","Cat_Y")],sep="-"))
#DataAggTot_NA$Pond_seas<-do.call(paste,c(DataAggTot_NA[c("Pond","Season_Year")],sep="-"))

#autoplot(pca_res2, data = DataAggTot_NA, colour = 'Cat_Y' ,frame = T)
#autoplot(pca_res2, data = DataAggTot_NA, colour = 'Pond' ,frame = T)
#autoplot(pca_res2, data = DataAggTot_NA, colour = 'Season_Year2' ,frame = T)
#autoplot(pca_res2, data = DataAggTot_NA, colour = 'Pond' ,frame = T)
#autoplot(pca_res2, data = DataAggTot_NA, colour = 'Treat_Cat_Y' ,frame = T)

library(corrplot)
###correlation and contribution of variable to PCA axis
corrplot(var$cos2, is.corr=F)
corrplot(var$contrib, is.corr=F)  

fviz_contrib(pca_res2, choice = "var", axes = 1, top = 10)+
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))
fviz_contrib(pca_res2, choice = "var", axes = 2, top = 10)+
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

fviz_contrib(pca_res, choice = "var", axes = 1:2, top = 10)


######pca variables plot of correlation factor and contribution
#fviz_pca_var(pca_res,col.var="cos2",axes = c(1, 2),
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
#        theme(text = element_text(size = 18),
#        axis.title = element_text(size = 18),
#        axis.text = element_text(size = 18))

fviz_pca_var(pca_res,col.var="contrib",axes = c(1, 2),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme(text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

##From PCA:
#Cond, DOC, K, pH, algae corr with axis 1
#Na, SO4, No3, PO4, Algae axis 2
#Cl, Wt, Ca, PO4, NH4

#add additional combinations to ENV dataset

DataAggEnv_NA$Season_Year_Month<-do.call(paste,c(DataAggEnv_NA[c("Season_Year","Month")],sep="-"))
DataAggEnv_NA$Season_Year_Month <- ordered(DataAggEnv_NA$Season_Year_Month, levels = c("W_17-12", "W_17-2", "SP_18-5", "SU_18-6", "SU_18-8", "A_18-10", "A_18-11",
"W_18-12", "W_18-1","SP_19-3", "SP_19-5", "SU_19-6"))



names(DataAggEnv_NA)


fviz_pca_biplot(pca_res2,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),#,"text"),
                pointshape=21,
                #habillage=DataAggEnv_NA$Season_Year_Month,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA$Season_Year_Month,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")

fviz_pca_biplot(pca_res2,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),#,"text"),
                pointshape=21,
                #habillage=DataAggTot_NA$Season_Year,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA$Season_Year,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")
#gradient.cols = c("red","yellow","springgreen","royalblue"))



####################season divided by first and second year
####first year
###homogenised month_season for the first year and second year to have also comparable months between the 2 years

head(DataAggEnv_NA)
DataAggEnv_NA_1stYear<-subset(DataAggEnv_NA,Season_Year=="W_17"|Season_Year== "SP_18"|Season_Year== "SU_18")
DataAggEnv_NA_1stYear$Season_Year_Month<-do.call(paste,c(DataAggEnv_NA_1stYear[c("Season_Year","Month")],sep="-"))
DataAggEnv_NA_1stYear<-subset(DataAggEnv_NA_1stYear,Season_Year_Month!="SU_18-8")#removed because there is sampling in August in 2019


names(DataAggEnv_NA_1stYear)
Env_NA_1stYear<- DataAggEnv_NA_1stYear[,c(11:26,28)]
LogEnv_NA_1stYear<-log10(Env_NA_1stYear+1)
names(LogEnv_NA_1stYear)
pca_res2_1stY <- prcomp(LogEnv_NA_1stYear,  scale.=TRUE,center=TRUE)


First<-fviz_pca_biplot(pca_res2_1stY,axes = c(1, 2),
                       #select.var = list(contrib = 7),
                       geom=c("point"),
                       pointshape=21,
                       #habillage=DataAggTot_NA$Season_Year,
                       #habillage = 12,
                       pointsize=4,
                       fill.ind=as.factor(DataAggEnv_NA_1stYear$Season_Year_Month),
                       col.ind="black",
                       #col.var="contrib",
                       arrowsize = 0.8,
                       labelsize=6,
                       #palette = "jco",
                       repel = FALSE,
                       addEllipses = TRUE,
                       #palette = "aaas",
                       ellipse.level=0.95,
                       ellipse.type="convex",
                       col.var = "red")



################################2ndyear


head(DataAggEnv_NA)
DataAggEnv_NA_2ndYear<-subset(DataAggEnv_NA,Season_Year=="W_18"|Season_Year== "SP_19"|Season_Year== "SU_19")
DataAggEnv_NA_2ndYear$Season_Year_Month<-do.call(paste,c(DataAggEnv_NA_2ndYear[c("Season_Year","Month")],sep="-"))
DataAggEnv_NA_2ndYear<-subset(DataAggEnv_NA_2ndYear,Season_Year_Month!="SP_19-3")##removed as no sampling in march 2018 (comparability)
names(DataAggEnv_NA_2ndYear)

Env_NA_2ndYear<- DataAggEnv_NA_2ndYear[,c(11:17, 19:26,28)] ###removed from env variables NO3 because all values costant to minimum detection
LogEnv_NA_2ndYear<-log10(Env_NA_2ndYear+1)
names(LogEnv_NA_2ndYear)
pca_res2_2ndY <- prcomp(LogEnv_NA_2ndYear,  scale.=TRUE,center=TRUE)

names(pca_res2_2ndY)

SecY<-fviz_pca_biplot(pca_res2_2ndY,axes = c(1, 2),
                      #select.var = list(contrib = 7),
                      geom=c("point"),
                      pointshape=21,
                      #habillage=DataAggTot_NA$Season_Year,
                      #habillage = 12,
                      pointsize=4,
                      fill.ind=as.factor(DataAggEnv_NA_2ndYear$Season_Year_Month),
                      col.ind="black",
                      #col.var="contrib",
                      arrowsize = 0.8,
                      labelsize=6,
                      #palette = "jco",
                      repel = FALSE,
                      addEllipses = TRUE,
                      #palette = "aaas",
                      ellipse.level=0.95,
                      ellipse.type="convex",
                      col.var = "red")
library(ggpubr)

ggarrange(First, SecY, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

###ponds

fviz_pca_biplot(pca_res2,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),
                pointshape=21,
                #habillage=DataAggTot_NA$Season_Year,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA$Pond,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")


###subsets Pond within years
#First year


First<-fviz_pca_biplot(pca_res2_1stY,axes = c(1, 2),
                     #select.var = list(contrib = 7),
                     geom=c("point"),
                     pointshape=21,
                     #habillage=DataAggTot_NA$Season_Year,
                     #habillage = 12,
                     pointsize=4,
                     fill.ind=DataAggEnv_NA_1stYear$Pond,
                     col.ind="black",
                     #col.var="contrib",
                     arrowsize = 0.8,
                     labelsize=6,
                     #palette = "jco",
                     repel = FALSE,
                     addEllipses = TRUE,
                     #palette = "aaas",
                     ellipse.level=0.95,
                     ellipse.type="convex",
                     col.var = "red")


###subsets Pond within seasons


SecY<-fviz_pca_biplot(pca_res2_2ndY,axes = c(1, 2),
                        #select.var = list(contrib = 7),
                        geom=c("point"),
                        pointshape=21,
                        #habillage=DataAggTot_NA$Season_Year,
                        #habillage = 12,
                        pointsize=4,
                        fill.ind=DataAggEnv_NA_2ndYear$Pond,
                        col.ind="black",
                        #col.var="contrib",
                        arrowsize = 0.8,
                        labelsize=6,
                        #palette = "jco",
                        repel = FALSE,
                        addEllipses = TRUE,
                        #palette = "aaas",
                        ellipse.level=0.95,
                        ellipse.type="convex",
                        col.var = "red")

library(ggpubr)

ggarrange(First, SecY, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

####Position in the mesocosm

fviz_pca_biplot(pca_res2,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),
                pointshape=21,
                #habillage=DataAggTot_NA$Season_Year,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA$Position,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")


###position 1st year

firstP<-fviz_pca_biplot(pca_res2_1stY,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),
                pointshape=21,
                #habillage=DataAggTot_NA$Season_Year,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA_1stYear$Position,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")


###position 2nd year

secondP<-fviz_pca_biplot(pca_res2_2ndY,axes = c(1, 2),
                #select.var = list(contrib = 7),
                geom=c("point"),
                pointshape=21,
                #habillage=DataAggTot_NA$Season_Year,
                #habillage = 12,
                pointsize=4,
                fill.ind=DataAggEnv_NA_2ndYear$Position,
                col.ind="black",
                #col.var="contrib",
                arrowsize = 0.8,
                labelsize=6,
                #palette = "jco",
                repel = FALSE,
                addEllipses = TRUE,
                #palette = "aaas",
                ellipse.level=0.95,
                ellipse.type="convex",
                col.var = "red")


ggarrange(firstP, secondP, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)


######RDA total ###############################
####done on reduced dataset
####removed Fall 2018 and not comparable months from the analysis to have comparable seasons among first and second year
names(DataAggEnv_NA)
DataAggEnv_NA$Season_Year_Month<-do.call(paste,c(DataAggEnv_NA[c("Season_Year","Month")],sep="-"))
DataAggEnv_NA_red<-subset(DataAggEnv_NA,Season_Year!="A_18")
DataAggEnv_NA_red<-subset(DataAggEnv_NA_red,Season_Year_Month!="SP_19-3")
DataAggEnv_NA_red<-subset(DataAggEnv_NA_red,Season_Year_Month!="SU_18-8")

DataAggEnv_NA_red$Season_Year<-as.factor(DataAggEnv_NA_red$Season_Year)
DataAggEnv_NA_red$Season_Year<-droplevels(DataAggEnv_NA_red$Season_Year)
DataAggEnv_NA_red$Season_Year_Month<-as.factor(DataAggEnv_NA_red$Season_Year_Month)
DataAggEnv_NA_red$Season_Year_Month<-droplevels(DataAggEnv_NA_red$Season_Year_Month)


DataAggEnv_NA_red$MonthFromBeg
DataAggEnv_NA_red$Cat_Y<-ifelse(DataAggEnv_NA_red$MonthFromBeg<8,"1st","2nd")

names(DataAggEnv_NA_red)

EnvAgg_NA_red<- DataAggEnv_NA_red[,c(11:26, 28)]                           
LogEnvAgg_NA_red<-log10(EnvAgg_NA_red+1)

corvif(LogEnvAgg_NA_red)
##although some collinearity among variables, as here they are response variable, I leave all IN

#####analysis with RDA
names(LogEnvAgg_NA_re)
row.names(LogEnvAgg_NA_red)
#Z_LogENVAgg_NA<-scale(LogEnvAgg_NA3, center=TRUE, scale=TRUE)


taxa.cap1<- rda(LogEnvAgg_NA_red~ Pond + Cat_Y*MonthFromBeg + Season,
                               scale.=TRUE, center=TRUE, data=DataAggEnv_NA_red)

vif.cca(taxa.cap1)

taxa.cap0 <- rda(LogEnvAgg_NA_red ~ 1, data=DataAggEnv_NA_red)
#step.env<-ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)
step.env <- ordistep(taxa.cap0, scope=formula(taxa.cap1),direction = 'forward', permutations = 999)

step.env
screeplot(step.env)
R2adj <- RsquareAdj(step.env)$adj.r.squared
R2adj
anova(step.env)
test_terms<-anova(step.env, by="margin", permu=999)
test_terms.adj<-test_terms
test_terms.adj$`Pr(>F)` <- p.adjust (test_terms$`Pr(>F)`, method = 'holm')
test_terms.adj

test_axis<-anova(step.env, by="axis", permu=999)
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

####plot

with(DataAggEnv_NA_red, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")


with(DataAggEnv_NA_red, levels(Season_Year))
colseas<- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)

with(DataAggEnv_NA_red, points(step.env, display = "sites", col=colvec[Pond],
                     pch = 21, bg = colvec[Pond]))
with(DataAggEnv_NA_red, colvec[Pond])
spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
with(DataAggEnv_NA_red, legend("topright", legend = levels(Pond), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,DataAggEnv_NA_red$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)

dev.off()

#tiff("rda_trout_bl2.tif", res=600, compression = "lzw", height=7, width=9, units="in")


####variance partitioning on significant (RDA selected) predictors
names(DataAggEnv_NA)

varp<-varpart (LogEnvAgg_NA_red, ~Pond,~MonthFromBeg*Cat_Y,~Season, data = DataAggEnv_NA_red,transfo="stand")
plot (varp, digits = 2, Xnames = c("Pond variability","Time","Seasonality"), bg = c('navy',"tomato",  "yellow"), cex=1.5)




###RDA first year (Win, Sp, Su 2018)
######RDA###############################
names(LogEnv_NA_1stYear)
DataAggEnv_NA_1stYear$Season_Year

corvif(LogEnv_NA_1stYear)

taxa.cap1<- rda(LogEnv_NA_1stYear~ Pond + MonthFromBeg + Season,
                scale=TRUE, center=TRUE, DataAggEnv_NA_1stYear)

vif.cca(taxa.cap1)

taxa.cap0 <- rda(LogEnv_NA_1stYear ~ 1, data=DataAggEnv_NA_1stYear)
#step.env<-ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)
step.env <- ordistep(taxa.cap0, scope=formula(taxa.cap1),direction = 'forward', permutations = 999)
step.env
screeplot(step.env)
R2adj <- RsquareAdj(step.env)$adj.r.squared
R2adj


anova(step.env)
test_terms<-anova(step.env, by="margin", permu=999)
test_terms.adj<-test_terms
test_terms.adj$`Pr(>F)` <- p.adjust (test_terms$`Pr(>F)`, method = 'holm')
test_terms.adj

test_axis<-anova(step.env, by="axis", permu=999)
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

with(DataAggEnv_NA_1stYear, levels(Pond))

colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")

DataAggEnv_NA_1stYear$Season_Year<-droplevels(DataAggEnv_NA_1stYear$Season_Year)

with(DataAggEnv_NA_1stYear, levels(Season_Year))
colseas<- c("red2", "chocolate","firebrick")

###plot rda

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5,xlim=c(-1, 1), ylim=c(-1, 1))
with(DataAggEnv_NA_1stYear, points(step.env, display = "sites", col=colvec[Pond],
                           pch = 21, bg = colvec[Pond]))
with(DataAggEnv_NA_1stYear, colvec[Pond])
spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
with(DataAggEnv_NA_1stYear, legend("topright", legend = levels(Pond), bty = "n",
                           col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,DataAggEnv_NA_1stYear$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.7, alpha=0.1)



#tiff("rda_trout_bl2.tif", res=600, compression = "lzw", height=7, width=9, units="in")

names(DataAggEnv_NA_1stYear)
varp<-varpart (LogEnv_NA_1stYear, ~Pond,~MonthFromBeg,~Season, data = DataAggEnv_NA_1stYear,transfo="stand")
plot (varp, digits = 2, Xnames = c("Pond variability","Time","Seasonality"), bg = c('navy',"tomato",  "yellow"))


###RDA second year (Win, SP, SUM 2019)

######RDA###############################
names(LogEnv_NA_2ndYear)
names(DataAggEnv_NA_2ndYear)

corvif(LogEnv_NA_2ndYear)


taxa.cap1<- rda(LogEnv_NA_2ndYear~ Pond + MonthFromBeg + Season ,
                scale.=TRUE, center=TRUE, data=DataAggEnv_NA_2ndYear)

vif.cca(taxa.cap1)

taxa.cap0 <- rda(LogEnv_NA_2ndYear ~ 1, data=DataAggEnv_NA_2ndYear)
#step.env<-ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)
step.env <- ordistep(taxa.cap0, scope=formula(taxa.cap1),direction = 'forward', permutations = 999)
step.env
screeplot(step.env)
R2adj <- RsquareAdj(step.env)$adj.r.squared
R2adj

anova(step.env)
test_terms<-anova(step.env, by="margin", permu=999)
test_terms.adj<-test_terms
test_terms.adj$`Pr(>F)` <- p.adjust (test_terms$`Pr(>F)`, method = 'holm')
test_terms.adj

test_axis<-anova(step.env, by="axis", permu=999)
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj


with(DataAggEnv_NA_2ndYear, levels(Pond))

colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")

with(DataAggEnv_NA_2ndYear, levels(Season_Year))
colseas<- c("red2", "chocolate","firebrick" )

###plot rda

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5, xlim=c(-1, 1), ylim=c(-1, 1))
with(DataAggEnv_NA_2ndYear, points(step.env, display = "sites", col=colvec[Pond],
                                   pch = 21, bg = colvec[Pond]))
with(DataAggEnv_NA_2ndYear, colvec[Pond])
spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
with(DataAggEnv_NA_2ndYear, legend("topright", legend = levels(Pond), bty = "n",
                                   col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,DataAggEnv_NA_2ndYear$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)


#tiff("rda_trout_bl2.tif", res=600, compression = "lzw", height=7, width=9, units="in")

names(DataAggEnv_NA)
varp<-varpart (LogEnv_NA_2ndYear, ~Pond,~MonthFromBeg,~Season, data = DataAggEnv_NA_2ndYear,transfo="stand")
plot (varp, digits = 2, Xnames = c("Pond variability","Time","Seasonality"), bg = c('navy',"tomato",  "yellow"))


