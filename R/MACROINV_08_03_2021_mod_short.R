### ------------------------ ###
### --- SHORT MZB SCRIPT --- ### 
### ------------------------ ###

# --------------- #
# date:  11.03.21
# files in 
# files out
# EREES POND
# Prepare mzb data  
# --------------- #

# 1. data og data set 
# 2. drop columns -> datared
# 3. Aggregate -> dataAggEnv

# read data ---------------------------------------------------------------

data = read.table("data/ENV_MACRO_08_03_21.csv", h = T, sep = ",")  # sites, abiotic, biotic data
data %<>% dplyr::rename(SubVeg = Algae)

# adjust date  ------------------------------------------------------
data %<>%
  mutate(
    newdate = strptime(as.character(Date), "%d/%m/%Y") %>%
      format(format = "%Y-%m-%d") %>%
      as.Date("%Y-%m-%d")
  )

#order by date
data%<>% arrange(newdate)

# Convert to factor  ------------------------------------------------------
data %<>%
  mutate (
    sampler_ID %<>% as.factor(),
    sampler_ID2 %<>% as.factor(),
    SeasonUnique %<>% as.factor(),
    Year %<>% as.factor(),
    Date2 %<>% as.factor(),
    Date3 %<>% as.factor(),
    Season %<>% as.factor(),
    Pond %<>% as.factor(),
    Treat %<>% as.factor()
    )

#re-order ponds from 1 to 12
data$Pond %<>% factor(levels = paste0("P", 1:12))

#re-order seasons over the years
data$Season_Year %<>%
  ordered(levels = c("W_17", "SP_18", "SU_18", "A_18", "W_18", "SP_19", "SU_19"))
data$Season_Year2 %<>% 
  ordered(levels = c("Win1", "Sp1", "Su1", "Au1", "Win2", "Sp2", "Su2"))

# Remove Variables with to many missing entries ---------------------------
# datared used to be created here 
data %<>% dplyr::select(
  !c(
    "AbuZoo",
    "TOC",
    "LeafDec",
    "LeafDec_mic",
    "LeafDec_shr",
    "GreenTea",
    "RedTea",
    "Anura_Bufonidae",
    "Anura_Ranidae",
    "Urodela_Salamandridae",
    "NI",
    "Coleoptera_Curculionidae_adult",
    "Lepidoptera_Crambidae",
    "Diptera_Chironomidae_pupae"
  )
)

# Aggregate MZB Subsamples ------------------------------------------------
fac_env  <- data[, c(23:40)]
fac_tax <- data[, c(41:89)]

agg_env = as.data.frame(aggregate(
  fac_env,
  by = list(
    data$Month,
    data$Year,
    data$Pond,
    data$Distance_lake,
    data$Position,
    data$newdate,
    data$DaysFromBeg,
    data$SeasonUnique,
    data$Season,
    data$Season_Year,
    data$Season_Year2,
    data$Treat
  ),
  FUN = mean,
  na.rm = T
))

agg_tax = as.data.frame(aggregate(
  fac_tax,
  by = list(
    data$Month,
    data$Year,
    data$Pond,
    data$Distance_lake,
    data$Position,
    data$newdate,
    data$DaysFromBeg,
    data$SeasonUnique,
    data$Season,
    data$Season_Year,
    data$Season_Year2,
    data$Treat
  ),
  FUN = sum,
  na.rm = T
))

# create pond date variable 
agg_env %<>% 
  mutate(Pond_dateEnv = paste0(Group.3, "-", Group.6))
agg_tax %<>% 
  mutate(Pond_dateTaxa = paste0(Group.3, "-", Group.6))

setnames(
  agg_env,
  old = c(1:12),
  new = c(
    'Month',
    'Year',
    'Pond',
    'Dist_Lake',
    'Position',
    'newdate',
    'DaysFromBeg',
    'SeasonUnique',
    "Season",
    'Season_Year',
    'Season_Year2',
    'Treat'
  )
)

setnames(
  agg_tax,
  old = c(1:12),
  new = c(
    'Month1',
    'Year1',
    'Pond1',
    'Dist_Lake1',
    'Position1',
    'newdate1',
    'DaysFromBeg1',
    'SeasonUnique1',
    "Season1",
    'Season_Year1',
    'Season_Year21',
    'Treat1'
  )
)

# create dataAggTot
agg_all <- merge(agg_env, agg_tax, by = "row.names")
agg_all %<>% arrange(newdate) 

# remove columns 
agg_all <- subset(agg_all, select = -c(1, 33:44, 94))
agg_all %<>% dplyr::rename(Pond_date = Pond_dateEnv)

# Remove Missing data from total dataset (env+mzb)
agg_all <- agg_all[complete.cases(agg_all),]
# drop A_18 
agg_all %<>% filter(Season_Year != "A_18")
agg_all %<>% mutate(Season_Year_Month = paste0(c("Season_Year","-", "Month")))

agg_all %<>%
  filter(!Season_Year_Month %in% c("SP_19-3", "SU_18-8")) %>%
  mutate(Season_Year = droplevels(as.factor(Season_Year)), 
         Season_Year_Month  = droplevels(as.factor(Season_Year_Month)),
         Cat_Y = ifelse(DaysFromBeg < 239, "1st", "2nd")
  )

# extract taxa from dataAggTot
agg_tax = agg_all[,c(32:80)]                           
# drop rows and column with 0 sums 
agg_tax %<>% dplyr::select(where(function(x) sum(x)>0))
agg_all %<>% dplyr::select(where(function(x) ifelse(is.numeric(x),sum(x)>0,TRUE)))
backup = agg_tax
agg_tax = agg_tax[rowSums(agg_tax) != 0, ]
agg_all = agg_all[rowSums(backup) != 0, ]

# extract environmental data from tot
agg_env <- agg_all[, c(13:28, 30)]

agg_env %<>% dplyr::select(!c("K", "pH", "DOC", "Ca"))
data_exp <- agg_all[,c(3, 7, 9, 80)]
saveRDS(data_exp, "data/taxa_exp.rds")
saveRDS(agg_env, "data/taxa_env.rds")
saveRDS(agg_tax, "data/taxa_bio.rds")

# FIT MODEL 1: MZB + ENV --------------------------------------------------
# taxa.cap1<- rda(Taxadec ~., EnvAgg_NA_red2,  center=T, scale.=TRUE)

# FIT MODEL 2: MZB ~ EXP --------------------------------------------------

###With Season and POND




taxa.cap1<- rda(Taxadec ~Pond+Season+DaysFromBeg*Cat_Y, EnvAgg_NA_red_pond,  center=T, scale.=TRUE)
vif.cca(taxa.cap1)

taxa.cap1<- rda(Taxadec ~., EnvAgg_NA_red_pond,  center=T, scale.=TRUE)
taxa.cap0 <- rda(Taxadec ~ 1, data=EnvAgg_NA_red_pond, scale=TRUE,center=TRUE)
#step.env <- ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)
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

dataAggTot_NA_red$Pond_Year<-do.call(paste,c(dataAggTot_NA_red[c("Pond","Cat_Y")],sep="-"))

with(dataAggTot_NA_red, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")


with(dataAggTot_NA_red, levels(Season_year))
colseas<- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)
with(dataAggTot_NA_red, points(step.env, display = "sites", col=colvec[Pond],
                               pch = 21, bg = colvec[Pond]))
with(dataAggTot_NA_red, colvec[Pond])
#spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
#arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
#text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
#with(dataAggTot_NA_red, legend("topright", legend = levels(Pond), bty = "n",
#                               col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,dataAggTot_NA_red$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)



varp<-varpart (Taxadec, ~Pond, ~Cat_Y, ~Season,  data = dataAggTot_NA_red)
plot (varp, digits = 1, Xnames = c('Pond diversity',"Time", "Season"), bg = c('navy', 'tomato', "yellow"))




summary(step.env)

#Dominant taxa
plot(step.env,type="n",xlim=c(-6,6), ylim=c(-2, 2),cex.axis=1, cex.lab=1, main="Dominant")
text(step.env, display = "species", col="blue", cex=1, select=c(3, 4,  6, 12:16, 18, 19 ,26:28, 30,31, 32),scaling=1, pos=3)
with(step.env, points(step.env, display = "species", select=c(3, 4,  6, 12:16, 18, 19 ,26:28, 30,31, 32),scaling=1, pch=21, bg="blue"))
text(step.env, display = "bp",col="red", cex=1.5, scaling=1)


####################RDA separately for 1st and 2nd year


###1st year

names(dataAggTot_NA_red)
dataAggTot_NA_red_1st<-subset(dataAggTot_NA_red,Cat_Y=="1st")
dataAggTot_NA_red_1st$Season_Year_Month<-droplevels(dataAggTot_NA_red_1st$Season_Year_Month)
names(dataAggTot_NA_red_1st)
TaxaAgg_NA_red_1st<- dataAggTot_NA_red_1st[,c(32:80)]                           

###remove rows and column = 0 in both the total dataset and the taxa dataset
TaxaAgg_NA_red_1st=TaxaAgg_NA_red_1st[,colSums(TaxaAgg_NA_red_1st)!=0]
dataAggTot_NA_red_1st=dataAggTot_NA_red_1st[rowSums(TaxaAgg_NA_red_1st)!=0,]
TaxaAgg_NA_red_1st=TaxaAgg_NA_red_1st[rowSums(TaxaAgg_NA_red_1st)!=0,]



row.names(TaxaAgg_NA_red_1st)
#dataAggTot_NA=dataAggTot_NA[-c(2, 3, 5:7, 11 ,14),]
names(TaxaAgg_NA_red_1st)
Taxadec1st=decostand(TaxaAgg_NA_red_1st, method="hell")
Taxa.dca <- decorana(Taxadec1st)


####fist Axis DCA <3 with hellinger transformed taxa dataset-> use RDA

###RDA TAXA and ENVIRONMENTAL variables

names(dataAggTot_NA_red_1st)
names(EnvAgg_NA_RDA_1st)
EnvAgg_NA_RDA_1st<-dataAggTot_NA_red_1st[,c(13:28, 30)]
corvif(EnvAgg_NA_RDA_1st)

EnvAgg_NA_RDA_1st2<-EnvAgg_NA_RDA_1st[, !(colnames(EnvAgg_NA_RDA_1st) %in% c("WatLev","K", "O2", "Ca", "pH", "Na"))]
corvif(EnvAgg_NA_RDA_1st2)

taxa.cap1<- rda(Taxadec1st ~., EnvAgg_NA_RDA_1st2,  center=T, scale.=TRUE)
vif.cca(taxa.cap1)

taxa.cap1<- rda(Taxadec1st ~., EnvAgg_NA_RDA_1st2,  center=T, scale.=TRUE)
taxa.cap0 <- rda(Taxadec1st ~ 1, data=EnvAgg_NA_RDA_1st2, scale.=TRUE,center=TRUE)
#step.env <- ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)
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

dataAggTot_NA_red_1st$Pond_Year<-do.call(paste,c(dataAggTot_NA_red_1st[c("Pond","Cat_Y")],sep="-"))

with(dataAggTot_NA_red_1st, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")


with(dataAggTot_NA_red_1st, levels(Season_Year))
colseas<- c("red2", "chocolate","firebrick")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)
with(dataAggTot_NA_red_1st, points(step.env, display = "sites", col=colvec[Pond],
                               pch = 21, bg = colvec[Pond]))
with(dataAggTot_NA_red_1st, colvec[Pond])
#spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
#arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
#text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
#with(dataAggTot_NA_red, legend("topright", legend = levels(Pond), bty = "n",
#                               col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,dataAggTot_NA_red_1st$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)


#Dominant taxa
summary(step.env)
plot(step.env,type="n",xlim=c(-6,6), ylim=c(-2, 2),cex.axis=1, cex.lab=1, main="Dominant")
text(step.env, display = "species", col="blue", cex=1, select=c(3, 4,  5, 10:13, 19:23),scaling=1, pos=3)
with(step.env, points(step.env, display = "species", select=c(3, 4,  5, 10:13, 19:23),scaling=1, pch=21, bg="blue"))
text(step.env, display = "bp",col="red", cex=1.5, scaling=1)


###variance partitioning
varp<-varpart (step.env, ~, ~,  data = dataAggTot_NA_red_1st)
plot (varp, digits = 1, Xnames = c('Abiotic',"Biotic"), bg = c('navy', 'tomato', "yellow"), cex=1.5)


###With Season and POND - 1st year

names(dataAggTot_NA_red_1st)
EnvAgg_NA_pond_1st<-dataAggTot_NA_red_1st[,c( 3, 7, 9)]
corvif(EnvAgg_NA_pond_1st)

taxa.cap1<- rda(Taxadec1st ~., EnvAgg_NA_pond_1st,  center=T, scale.=TRUE)
vif.cca(taxa.cap1)

taxa.cap1<- rda(Taxadec1st ~Pond + DaysFromBeg+Season, EnvAgg_NA_pond_1st,  center=T, scale.=TRUE)
taxa.cap0 <- rda(Taxadec1st ~ 1, data=EnvAgg_NA_pond_1st, scale=TRUE,center=TRUE)
#step.env <- ordiR2step(taxa.cap0, taxa.cap1, trace=FALSE)

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

dataAggTot_NA_red_1st$Pond_Year<-do.call(paste,c(dataAggTot_NA_red_1st[c("Pond","Cat_Y")],sep="-"))

with(dataAggTot_NA_red_1st, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")
with(dataAggTot_NA_red_1st, levels(Season_Year))
colseas<- c("coral2", "gold")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)
with(dataAggTot_NA_red_1st, points(step.env, display = "sites", col=colvec[Pond],
                               pch = 21, bg = colvec[Pond]))
with(dataAggTot_NA_red_1st, colvec[Pond])
#spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
#arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
#text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
#with(dataAggTot_NA_red, legend("topright", legend = levels(Pond), bty = "n",
#                               col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,dataAggTot_NA_red_1st$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)



varp<-varpart (Taxadec1st, ~DaysFromBeg, ~Season,  data = dataAggTot_NA_red_1st)
plot (varp, digits = 1, Xnames = c("Time", "Season"), bg = c('navy', 'tomato', "yellow"), cex=1.5)


summary(step.env)

#Dominant taxa
plot(step.env,type="n",xlim=c(-6,6), ylim=c(-2, 2),cex.axis=1, cex.lab=1, main="Dominant")
text(step.env, display = "species", col="blue", cex=1, select=c(4, 6, 12:17, 25:30),scaling=1, pos=3)
with(step.env, points(step.env, display = "species", select=c(4, 6, 12:17, 25:30),scaling=1, pch=21, bg="blue"))
text(step.env, display = "bp",col="red", cex=1.5, scaling=1)



####################RDA 2nd year

names(dataAggTot_NA_red)
dataAggTot_NA_red_2nd<-subset(dataAggTot_NA_red,Cat_Y=="2nd")
names(dataAggTot_NA_red_2nd)
TaxaAgg_NA_red_2nd<- dataAggTot_NA_red_2nd[,c(32:80)]                           

###remove rows and column = 0 in both the total dataset and the taxa dataset
TaxaAgg_NA_red_2nd=TaxaAgg_NA_red_2nd[,colSums(TaxaAgg_NA_red_2nd)!=0]
dataAggTot_NA_red_2nd=dataAggTot_NA_red_2nd[rowSums(TaxaAgg_NA_red_2nd)!=0,]
TaxaAgg_NA_red_2nd=TaxaAgg_NA_red_2nd[rowSums(TaxaAgg_NA_red_2nd)!=0,]



row.names(TaxaAgg_NA_red_2nd)
#dataAggTot_NA=dataAggTot_NA[-c(2, 3, 5:7, 11 ,14),]
names(TaxaAgg_NA_red_2nd)
Taxadec2nd=decostand(TaxaAgg_NA_red_2nd, method="hell")
Taxa.dca <- decorana(Taxadec2nd)


####fist Axis DCA <3 with hellinger transformed taxa dataset-> use RDA

###RDA TAXA and ENVIRONMENTAL variables

names(dataAggTot_NA_red_2nd)
EnvAgg_NA_RDA_2nd<-dataAggTot_NA_red_2nd[,c(13:27, 30)]
corvif(EnvAgg_NA_RDA_2nd)

EnvAgg_NA_RDA_2nd2<-EnvAgg_NA_RDA_2nd[, !(colnames(EnvAgg_NA_RDA_2nd) %in% c("K", "pH", "Ca", "DOC"))]###remove colinear variable
corvif(EnvAgg_NA_RDA_2nd2)

taxa.cap1<- rda(Taxadec2nd ~., EnvAgg_NA_RDA_2nd2,  center=T, scale.=TRUE)
vif.cca(taxa.cap1)

taxa.cap1<- rda(Taxadec2nd ~., EnvAgg_NA_RDA_2nd2,  center=T, scale.=TRUE)
taxa.cap0 <- rda(Taxadec2nd ~ 1, data=EnvAgg_NA_RDA_2nd2, scale=TRUE,center=TRUE)
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

dataAggTot_NA_red_2nd$Pond_Year<-do.call(paste,c(dataAggTot_NA_red_2nd[c("Pond","Cat_Y")],sep="-"))

with(dataAggTot_NA_red_2nd, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")
with(dataAggTot_NA_red_2nd, levels(Cat_Y))
colseas<- c("coral2", "gold")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)
with(dataAggTot_NA_red_2nd, points(step.env, display = "sites", col=colvec[Pond],
                                   pch = 21, bg = colvec[Pond]))
with(dataAggTot_NA_red_2nd, colvec[Pond])
#spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
#arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
#text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
#with(dataAggTot_NA_red, legend("topright", legend = levels(Pond), bty = "n",
#                               col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,dataAggTot_NA_red_2nd$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)


#Dominant taxa

summary(step.env)
plot(step.env,type="n",xlim=c(-6,6), ylim=c(-2, 2),cex.axis=1, cex.lab=1, main="Dominant")
text(step.env, display = "species", col="blue", cex=1, select=c(2, 4, 6:10, 12:13,16:18, 20:22),scaling=1, pos=3)
with(step.env, points(step.env, display = "species", select=c(2, 4, 6:10, 12:13,16:18, 20:22),scaling=1, pch=21, bg="blue"))
text(step.env, display = "bp",col="red", cex=1.5, scaling=1)


###variance partitioning
varp<-varpart (Taxadec2nd, ~FL+NH4+SO4, ~SubVeg,~WT,  data = dataAggTot_NA_red_2nd)
plot (varp, digits = 1, Xnames = c('Abiotic',"Biotic", "Water Temperature"), bg = c('navy', 'tomato', "yellow"), cex=1.5)

###With Season and POND - 2nd year

names(dataAggTot_NA_red_2nd)
EnvAgg_NA_pond_2nd<-dataAggTot_NA_red_2nd[,c(3:5, 7, 9)]
corvif(EnvAgg_NA_pond_2nd)

taxa.cap1<- rda(Taxadec2nd ~., EnvAgg_NA_pond_2nd,  center=T, scale.=TRUE)
vif.cca(taxa.cap1)

taxa.cap1<- rda(Taxadec2nd ~., EnvAgg_NA_pond_2nd,  center=T, scale.=TRUE)
taxa.cap0 <- rda(Taxadec2nd ~ 1, data=EnvAgg_NA_pond_2nd, scale=TRUE,center=TRUE)
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

dataAggTot_NA_red_2nd$Pond_Year<-do.call(paste,c(dataAggTot_NA_red_2nd[c("Pond","Cat_Y")],sep="-"))

with(dataAggTot_NA_red_2nd, levels(Pond))
colvec <- c("red2", "chocolate","firebrick", "orchid1",
            "green4","springgreen4","chartreuse1","seagreen1",
            "mediumblue", "cadetblue", "blue", "steelblue")
with(dataAggTot_NA_red_2nd, levels(Cat_Y))
colseas<- c("coral2", "gold")

tiff("rda.tif", res=600, compression = "lzw", height=7, width=9, units="in")

plot(step.env,type="n",cex.axis=1.5, cex.lab=1.5)
with(dataAggTot_NA_red_2nd, points(step.env, display = "sites", col=colvec[Pond],
                                   pch = 21, bg = colvec[Pond]))
with(dataAggTot_NA_red_2nd, colvec[Pond])
#spe.sc<-scores(step.env, choices=1:2, scaling=2, display="sp")
#arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='darkcyan')
#text(step.env, display = "species", scaling = 2, cex = 1, col = "darkcyan")
#with(dataAggTot_NA_red, legend("topright", legend = levels(Pond), bty = "n",
#                               col = colvec, pch = 21, pt.bg = colvec))
text(step.env, display = "cn",col="blue", cex=1, scaling=2)
ordiellipse(step.env,dataAggTot_NA_red_2nd$Season_Year, draw="polygon",conf=0.95, 
            col=colseas, kind = "ehull", lwd=2, label=T, cex=0.5, alpha=0.1)



varp<-varpart (Taxadec2nd,~Season,~Pond , data = dataAggTot_NA_red_2nd)
plot (varp, digits = 1, Xnames = c( "Season", "Pond variability"), bg = c('navy', 'tomato', "yellow"), cex=1.5)


summary(step.env)

#Dominant taxa
plot(step.env,type="n",xlim=c(-6,6), ylim=c(-2, 2),cex.axis=1, cex.lab=1, main="Dominant")
text(step.env, display = "species", col="blue", cex=1, select=c(2, 4, 6:10, 12:13,16:18, 20:22),scaling=1, pos=3)
with(step.env, points(step.env, display = "species", select=c(2, 4, 6:10, 12:13,16:18, 20:22),scaling=1, pch=21, bg="blue"))
text(step.env, display = "bp",col="red", cex=1.5, scaling=1)


#####indicator species analysis

library(indicspecies)
indval = multipatt(Taxadec, dataAggTot_NA_red$Cat_Y,control = how(nperm=999), duleg=TRUE)
summary(indval, alpha=0.05)

# Beta-disper
####################b-diversity distance from centromer########
##also use to check homoeneity of variance among groups betadisper p<0.05

diss<-vegdist(Taxadec)
bd1<-betadisper(diss,dataAggTot_NA_red$Season_Year, type = "median")
scores(bd1)


tiff("distance_Centromer_lateral.tif", res=600, compression = "lzw", height=5, width=5, units="in")

boxplot(bd1$distances~dataAggTot_NA_red$Season_Year, xlab="", ylab="Distance centromer", main="B-disper", prob=TRUE, 
        cex.lab=1, cex.axis=1, cex.main=1, cex.sub=0.5,las=1)

###test for significance
permutest(bd1)
anova(bd1)
plot(bd1)
lm1<-permutest(bd1, pairwise=TRUE, p.adjust="BY")
lm1

plot(bd1$distances~dataAggTot_NA_red$DaysFromBeg)
abline(lm(bd1$distances~dataAggTot_NA_red$DaysFromBeg), col="blue")
lm1<-lm(bd1$distances~dataAggTot_NA_red$DaysFromBeg)
summary(lm1)



########################UNTILL HERE#######################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









######################nMDS-perMANOVA
mds.diss=metaMDS(Taxadec,dist="bray")

###NMDS
tiff("nMDS_Bti.tif", res=900, compression = "lzw", height=7, width=11, units="in")

#x11(); par(mfrow=c(2,2));

plot(mds.diss, type = "n",las=1, ylab="", axes=F,cex.axis=1.5, cex.main=1.5, cex.lab= 1)
mds.diss
points(mds.diss, pch=16,col=as.numeric(dataAggTot_NA_red$Season_Year,cex=1.5))

identify(mds.diss$points)
taxadec2<-Taxadec[-c(5),] 
dataAggTot_NA_red2<-dataAggTot_NA_red[-c(5),] 
mds.diss=metaMDS(taxadec2,dist="bray")


plot(mds.diss, type = "n",las=1, ylab="", axes=F,cex.axis=1.5, cex.main=1.5, cex.lab= 1)
mds.diss
points(mds.diss, pch=16,col=as.numeric(dataAggTot_NA_red2$Season_Year),cex=1.5)
axis(1,cex.axis=1)
axis(2,cex.axis=1)

#legend("bottomleft",legend=c("C","T"),col=c(1,2),pch=c(16,16),cex=1.2, bty="n" )

#groupz <- sort(unique(dataAgg2$Treat))
#for(i in seq(groupz)) { 
#  ordiellipse(mds.diss, dataAgg2$Treat, kind="se", conf=0.95, label=T, 
#              font=2, cex=1.5, col="black", show.groups=groupz[i]) 
#}

plot.new()
plot.window(xlim=c(-2,2), ylim=c(-1.5,1.5))
#plot(mds.diss,las=1, ylab="NDMS2", xlab="NMDS1", type="n",axes=F,cex.axis=1.5, cex.main=2, cex.lab= 1)
#points(mds.diss, pch = 21, cex=1, col="gray", bg="gray", lwd=1)
#points(mds.diss, display="species", pch = 21, cex=1, col="black", bg="white", lwd=2)
#axis(1,cex.axis=1)
#axis(2,cex.axis=1)
points(mds.diss, pch=16,col=as.numeric(dataAggTot_NA_red2$Season_Year),cex=1.5)
axis(1,cex.axis=1)
axis(2,cex.axis=1)
#ord.fit <- envfit(mds.diss ~ SubVeg+WT+Cond+pH, data=dataAggTot_NA_red, perm=999)
#plot(ord.fit, cex=1.5, size=3)

ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T,
            font=1, cex=0.6, draw="polygon", col="Blue", show.groups="W_17", alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="green", show.groups="SP_18",alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="red", show.groups="SU_18",alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="orange", show.groups="A_18",alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="darkgreen", show.groups="W_18",alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="darkblue", show.groups="SP_19",alpha = 0.2)
ordiellipse(mds.diss, dataAggTot_NA_red2$Season_Year, kind="se", conf=0.95, label=T, 
            font=1, cex=0.6, draw="polygon", col="darkred", show.groups="SU_19",alpha = 0.2)

#orditorp(mds.diss,display="species",col="red",cex=1)

######PONDS

plot.new()
plot.window(xlim=c(-2,2), ylim=c(-1.5,1.5))
#plot(mds.diss,las=1, ylab="NDMS2", xlab="NMDS1", type="n",axes=F,cex.axis=1.5, cex.main=2, cex.lab= 1)
#points(mds.diss, pch = 21, cex=1, col="gray", bg="gray", lwd=1)
#points(mds.diss, display="species", pch = 21, cex=1, col="black", bg="white", lwd=2)
#axis(1,cex.axis=1)
#axis(2,cex.axis=1)
points(mds.diss, pch=16,col=as.numeric(dataAggTot_NA_red2$Season_Year),cex=1.5)
axis(1,cex.axis=1)
axis(2,cex.axis=1)
#ord.fit <- envfit(mds.diss ~ SubVeg+WT+Cond+pH, data=dataAggTot_NA_red, perm=999)
#plot(ord.fit, cex=1.5, size=3)

ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T,
            font=1, cex=0.8, draw="polygon", col="Blue", show.groups="P1")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="green", show.groups="P2")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="red", show.groups="P3")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="orange", show.groups="P4")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="darkgreen", show.groups="P5")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="darkblue", show.groups="P6")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="darkred", show.groups="P7")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="yellow", show.groups="P8")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="brown", show.groups="P9")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="gray", show.groups="P10")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="gold", show.groups="P11")
ordiellipse(mds.diss, dataAggTot_NA_red2$Pond, kind="se", conf=0.95, label=T, 
            font=1, cex=0.8, draw="polygon", col="red", show.groups="P12")

#ordisurf(mds.diss, dataAggTot_NA_red2$DaysFromBeg, add=TRUE)


dev.off()


dev.off()


#############ADONIS - perMANOVA analysis of the variance

names(dataAgg2)


diss<-vegdist(taxadec2, method = "bray")
adonis(diss~Pond*Cat_Y+DaysFromBeg+Season, data=dataAggTot_NA_red2)

###similarity % analysis (Simper)

sim<-simper(taxadec2, dataAggTot_NA_red2$Cat_Y)
summary(sim)
sink()





#####over time of dominant taxa

names(dataAggTot)
long_data_sel<-dataAggTot%>%
  pivot_longer(c(34, 35, 38, 45:49, 51:52, 62:64, 66:68, 73, 81:83), names_to = "variables", values_to = "value")

long_data_sel$variables<-as.factor(long_data_sel$variables)
levels(long_data_sel$variables)       
long_data_sel$variables<-factor(long_data_sel$variables, levels = c("Abund","Rich", "Sh",
                                                                    "Basommatophora_Lymnaeidae", "Basommatophora_Planorbidae",
                                                                    "Coleoptera_Dytiscidae","Diptera_Chaoboridae","Diptera_Chironomidae","Diptera_Culicidae","Diptera_Tabanidae",
                                                                    "Ephemeroptera", "Ephemeroptera_Baetidae","Ephemeroptera_Caenidae", 
                                                                    "Odonata","Odonata_Aeshnidae","Odonata_Coenagrionidae","Odonata_Libellulidae","Odonata_Zygoptera",
                                                                    "Oligochaeta"))

tiff("All_Treat.tif", res=600, compression = "lzw", height=13, width=20, units="in")

ggplot(long_data_sel,aes(x=newdate, y = value)) +
  geom_point() +
  geom_smooth(aes())+
  facet_wrap(~variables, scales = "free", ncol = 4, ) +
  scale_x_date(date_labels = "%b-%y",date_breaks = "3 months")+
  labs(x = "Days from Beg", y = "response")+
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold.italic"))

dev.off()

tiff("All_Ponds.tif", res=600, compression = "lzw", height=13, width=20, units="in")

ggplot(long_data_sel,aes(x=newdate, y = value)) +
  geom_point(alpha=0.2) +
  geom_smooth(aes(col=Pond), size=0.5,se=F)+
  geom_smooth(fill="red", col="red", linetype = "dashed", size=1.4, se=F)+
  facet_wrap(.~variables, scales = "free", ncol = 4) +
  scale_x_date(date_labels = "%b-%y",date_breaks = "3 months")+
  labs(x = "Days from Beg", y = "response")+
  theme(strip.text.x = element_text(size = 12, 
                                    color = "black", face = "bold.italic"))
dev.off()


#measure variation of sd (among ponds) over time

sum_long <- summarySE(long_data_sel, measurevar="value", groupvars=c("newdate","SeasonUnique" , "DaysFromBeg", "variables" ))
str(sum_long)

tiff("All_sd.tif", res=600, compression = "lzw", height=13, width=20, units="in")

ggplot(sum_long,aes(x=DaysFromBeg, y = sd)) +
  geom_point() +
  geom_smooth()+
  geom_smooth(  method = lm, col="red", se=F)+
  facet_wrap(.~variables, scales = "free", ncol = 4) +
  #scale_x_date(date_labels = "%b-%y",date_breaks = "3 months")+
  labs(x = "Days from Beg", y = "response")+
  theme(strip.text.x = element_text(size = 11, 
                                    color = "black", face = "bold.italic"))
dev.off()



################Modeling each variable VS Treatment (Control-Bti)
names(dataAgg)
names(long_data)

#library(lattice)
#xyplot(Algae~DaysFromBeg|Pond,
#data=dataAgg,
#panel=function(x,y){
#panel.points(x,y, col=1)
#panel.loess(x,y, span=0.5,col=1, lwd=2)
#})

algae1<-ggplot(dataAgg, aes(DaysFromBeg, Algae))+
  geom_point(aes())+
  geom_smooth(aes( fill="red", alpha=0.7),method = "lm", lty=1)+
  geom_smooth(aes( fill="blue", col="red",alpha=0.7),method="loess", span=0.7,formula=y~x, se=T,  lty=1)

algae1+facet_wrap(.~Treat, ncol=2)




dataAgg$DaysFromBeg
names(data)

lm <- lme(Algae ~ Treat*DaysFromBeg+
            Treat*Cat_Y,
          #Pond*DaysFromBeg+
          #Pond*Cat_Y,
          random=~1|Pond, 
          method="REML",
          data=dataAgg,
          na.action = na.omit)

#lm1 <- update(lm, weights=varIdent(form=~1|Pond))
#lm2 <- update(lm, weights=varFixed(~DaysFromBeg))
#lm3 <- update(lm, weights=varComb(varIdent(form=~1|Pond),
#                                 varFixed(~DaysFromBeg)))
#lm4 <- update(lm, weights=varComb(varIdent(form=~1|Pond),
#                                  varPower(form=~DaysFromBeg)))
#lm4 <- update(lm, weights=varComb(varIdent(form=~1|Pond),
#varPower(form=~DaysFromBeg)))
#AIC(lm, lm1, lm2, lm3, lm4)

Resid<-resid(lm,type="normalized")
par(mfrow=c(2,2), mar=c(3,3,3,3))
qqnorm(Resid)
qqline(Resid, col="red")
plot(fitted(lm), Resid)
#hist(Resid)
boxplot(Resid~Cat_Y, data=dataAgg)
#boxplot(Resid~Pond, data=dataAgg)
plot(DaysFromBeg~Resid, data=dataAgg)
acf(residuals(lm), lag.max=10,main = "")
pacf(residuals(lm4), lag.max=10,main = "")


lm1 <- lme(Algae ~ Treat*DaysFromBeg+
             Treat*Cat_Y,
           random=~1|Pond, 
           method="REML",
           data=dataAgg,
           na.action = na.omit)

###testing for serial autocorrelation (corARMA)

M2 <- update(lm1, cor=corARMA(p=1, form=~1|Pond))
M3 <- update(lm1, cor=corARMA(p=2, form=~1|Pond))
M4 <- update(lm1, cor=corARMA(p=3, form=~1|Pond))
M5  <- update(lm1, cor=corARMA(p=4, form=~1|Pond))
M6  <- update(lm1, cor=corARMA(p=5, form=~1|Pond))

M7  <- update(lm1, cor=corARMA(q=1, form=~1|Pond))
M8  <- update(lm1, cor=corARMA(q=2, form=~1|Pond))
M9  <- update(lm1, cor=corARMA(q=3, form=~1|Pond))
M10  <- update(lm1, cor=corARMA(q=4, form=~1|Pond))
M11 <- update(lm1, cor=corARMA(q=5, form=~1|Pond))

M12  <- update(lm1, cor=corARMA(p=1,q=1, form=~1|Pond))
M13  <- update(lm1, cor=corARMA(p=2,q=1, form=~1|Pond))
M14  <- update(lm1, cor=corARMA(p=3,q=1, form=~1|Pond))
M15  <- update(lm1, cor=corARMA(p=4,q=1, form=~1|Pond))
#M16  <- update(lm1, cor=corARMA(p=5, q=1,form=~1|Pond))
M17  <- update(lm1, cor=corARMA(p=1, q=2,form=~1|Pond))
M18  <- update(lm1, cor=corARMA(p=2,q=2, form=~1|Pond))
M19  <- update(lm1, cor=corARMA(p=3,q=2, form=~1|Pond))
#M20 <- update(lm1, cor=corARMA(p=4,q=2, form=~1|Pond))
#M21 <- update(lm1, cor=corARMA(p=5,q=2, form=~1|Pond))
M22 <- update(lm1, cor=corARMA(p=1,q=3, form=~1|Pond))
#M23 <- update(lm1, cor=corARMA(p=2,q=3, form=~1|Pond))
M24 <- update(lm1, cor=corARMA(p=3,q=3, form=~1|Pond))
#M25 <- update(lm1, cor=corARMA(p=4,q=3, form=~1|Pond))
#M26 <- update(lm1, cor=corARMA(p=5,q=3, form=~1|Pond))
M27 <- update(lm1, cor=corARMA(p=1,q=4, form=~1|Pond))
M28 <- update(lm1, cor=corARMA(p=2,q=4, form=~1|Pond))
#M29 <- update(lm1, cor=corARMA(p=3,q=4, form=~1|Pond))
#M30 <- update(lm1, cor=corARMA(p=4,q=4, form=~1|Pond))
#M31 <- update(lm1, cor=corARMA(p=5,q=4, form=~1|Pond))
M32 <- update(lm1, cor=corARMA(p=1,q=5, form=~1|Pond))
#M33 <- update(lm1, cor=corARMA(p=2,q=5, form=~1|Pond))
#M34 <- update(lm1, cor=corARMA(p=3,q=5, form=~1|Pond))
#M35 <- update(lm1, cor=corARMA(p=4,q=5, form=~1|Pond))
#M36 <- update(lm1, cor=corARMA(p=5,q=5, form=~1|Pond))


AIC(lm1,M2,M3, M4, M5, M6, M7,M8,M9,M10,
    M11,M12, M13, M14, M15, M17,M18, M19, 
    M22, M24, M27, M28, M32)

M2 <- lme(Algae ~ Treat*DaysFromBeg+
            Treat*Season_Year,
          random=~1|Pond,
          method="REML",
          cor=corARMA(p=1,q=0, form=~1|Pond),
          data=dataAgg,
          na.action = na.omit)

summary(M2)
Resid<-resid(M2, type="normalized")
par(mfrow=c(3,2))
qqnorm(Resid)
qqline(Resid, col="red")
plot(fitted(M2),Resid)
#hist(Resid)
boxplot(Resid~Cat_Y, data=dataAgg)
boxplot(Resid~as.factor(Pond), data=dataAgg)
#boxplot(Resid~Pond, data=dataAgg)
plot(Resid~DaysFromBeg, data=dataAgg)
acf(residuals(M2), lag.max=10,main = "")
pacf(residuals(M2), lag.max=10,main = "")

Eall<-vector(length=length(dataAgg$Algae))
Eall[]<-NA
I1<-!is.na(dataAgg$Algae)
Eall[I1]<-Resid
xyplot(Eall~dataAgg$DaysFromBeg|dataAgg$Pond, col=1)


E1<-Eall[data$Pond=="AP1"]
E2<-Eall[data$Pond=="AP2"]
E3<-Eall[data$Pond=="AP3"]
E4<-Eall[data$Pond=="AP4"]
E5<-Eall[data$Pond=="AP5"]
E6<-Eall[data$Pond=="AP6"]
E7<-Eall[data$Pond=="AP7"]
E8<-Eall[data$Pond=="AP8"]
E9<-Eall[data$Pond=="AP9"]
E10<-Eall[data$Pond=="AP10"]
E11<-Eall[data$Pond=="AP11"]
E12<-Eall[data$Pond=="AP12"]


par(mfrow=c(4,3), mar=c(3,3,3,1))
acf(E1, na.action=na.pass)
acf(E2, na.action=na.pass)
acf(E3, na.action=na.pass)
acf(E4, na.action=na.pass)
acf(E5, na.action=na.pass)
acf(E6, na.action=na.pass)
acf(E7, na.action=na.pass)
acf(E8, na.action=na.pass)
acf(E9, na.action=na.pass)
acf(E10, na.action=na.pass)
acf(E11, na.action=na.pass)
acf(E12, na.action=na.pass)


D<-cbind(E1,E2, E3, E4, E5, E6,
         E7, E8, E9, E10, E11, E12)

cor<-cor(D, use="pairwise.complete.obs")


#select factors
M2 <- lme(Algae ~ Treat*DaysFromBeg+
            Treat*Cat_Y,
          random=~1|Pond, 
          method="ML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)


m1<-update(M2,.~.-DaysFromBeg:Treat)
m2<-update(M2,.~.-Treat:Cat_Y)

anova(M2,m1)
anova(M2,m2)


M2 <- lme(Algae ~ DaysFromBeg+Treat+Cat_Y,
          random=~1|Pond, 
          method="ML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)


m3<-update(M2,.~.- DaysFromBeg)
m4<-update(M2,.~.-Treat)
m5<-update(M2,.~.-Cat_Y)
anova(M2,m3)
anova(M2,m4)
anova(M2,m5)


boxplot(Algae~Cat_Y, data=dataAgg)

M2 <- lme(Algae ~ DaysFromBeg,
          random=~1|Pond, 
          method="REML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)



summary(M2)
anova(M2)
r.squaredGLMM(M2)

Fit<- fitted(M2, level=0)
names(Fit)

library(ggeffects)
PDpr <- ggpredict(M2, "DaysFromBeg")
plot(PDpr)
PDpr<-as.data.frame(PDpr)

PDpr<-rename(PDpr, c("x"="DaysFromBeg", "predicted"="FitAlgae","std.error"="se_fit","conf.low"="LOW_fit","conf.high"="High_fit","group"="group_fit"  ))

fit1<-ggplot(dataAgg, aes(DaysFromBeg, Fit))+
  geom_point(aes(DaysFromBeg, Algae))+
  #geom_smooth(aes(col="gray"),col="red",method = "lm", lty=1)
  geom_smooth( aes(y=Algae,col="gray"),col="red",se=F, method = "lm", lty=1)+
  geom_smooth(aes(y=Fit, col="gray"),method="lm", se=F,col="blue", lty=1)
fit1+facet_wrap(~Pond, ncol=2, nrow=6)

Resid<-resid(M2, type="normalized")
par(mfrow=c(3,2))
qqnorm(Resid)
qqline(Resid, col="red")
plot(fitted(M2),Resid)
#hist(Resid)
boxplot(Resid~Treat*Year, data=dataAgg)
boxplot(Resid~as.factor(DaysFromBeg), data=dataAgg)
#boxplot(Resid~Pond, data=dataAgg)
plot(Resid~DaysFromBeg, data=dataAgg)
acf(residuals(M2), lag.max=10,main = "")
pacf(residuals(M2), lag.max=10,main = "")



################model Ponds
names(dataAgg)

#library(lattice)
#xyplot(Algae~DaysFromBeg|Pond,
#data=dataAgg,
#panel=function(x,y){
#panel.points(x,y, col=1)
#panel.loess(x,y, span=0.5,col=1, lwd=2)
#})

algae1<-ggplot(dataAgg, aes(DaysFromBeg, Algae))+
  geom_point(aes())+
  geom_smooth(aes( fill="red", alpha=0.7),method = "lm", lty=1)
#geom_smooth(aes( fill="blue", col="red",alpha=0.7),method="loess", span=0.7,formula=y~x, se=T,  lty=1)

algae1+facet_wrap(.~Pond, ncol=4)

names(dataAgg)

gls1 <- gls(Algae ~ Pond*Cat_Y+
              Pond*DaysFromBeg ,
            #random=~1|newdate, 
            method="REML",
            data=dataAgg,
            na.action = na.omit)
anova(gls1)

boxplot(Algae~Pond*Cat_Y, data=dataAgg)

Resid<-resid(gls1, type="normalized")
par(mfrow=c(3,2), mar=c(3,3,3,1))
qqnorm(Resid)
qqline(Resid, col="red")
plot(fitted(gls1), Resid)
#hist(Resid)
boxplot(Resid~Treat, data=dataAgg)
boxplot(Resid~as.factor(Pond), data=dataAgg)
#boxplot(Resid~Pond, data=dataAgg)
plot(DaysFromBeg~Resid, data=dataAgg)
acf(Resid, lag.max=20,main = "")
pacf(Resid, lag.max=20,main = "")



lm1 <- update(gls1, weights=varIdent(form=~1|Pond))
lm2 <- update(gls1, weights=varFixed(~DaysFromBeg))
lm3 <- update(gls1, weights=varComb(varIdent(form=~1|Pond),
                                    varFixed(~DaysFromBeg)))
lm4 <- update(gls1, weights=varComb(varIdent(form=~1|Pond),
                                    varIdent(form=~1|Cat_Y)))
lm4a <- update(gls1, weights=varComb(varIdent(form=~1|Pond),
                                     varPower(form=~DaysFromBeg)))
lm5 <- update(gls1, weights=varComb(varIdent(form=~1|Pond),
                                    varPower(form=~DaysFromBeg|Pond)))

AIC(gls1, lm1, lm2, lm3,  lm4, lm4a,lm5)





###testing for serial autocorrelation (corARMA)

M2 <- update(lm4, cor=corARMA(p=1, form=~1|Pond))
M3 <- update(lm4, cor=corARMA(p=2, form=~1|Pond))
M4 <- update(lm4, cor=corARMA(p=3, form=~1|Pond))
M5  <- update(lm4, cor=corARMA(p=4, form=~1|Pond))
M6  <- update(lm4, cor=corARMA(p=5, form=~1|Pond))

M7  <- update(lm4, cor=corARMA(q=1, form=~1|Pond))
M8  <- update(lm4, cor=corARMA(q=2, form=~1|Pond))
M9  <- update(lm4, cor=corARMA(q=3, form=~1|Pond))
M10  <- update(lm4, cor=corARMA(q=4, form=~1|Pond))
M11 <- update(lm4, cor=corARMA(q=5, form=~1|Pond))

M12  <- update(lm4, cor=corARMA(p=1,q=1, form=~1|Pond))
M13  <- update(lm4, cor=corARMA(p=2,q=1, form=~1|Pond))
M14  <- update(lm4, cor=corARMA(p=3,q=1, form=~1|Pond))
M15  <- update(lm4, cor=corARMA(p=4,q=1, form=~1|Pond))
M16  <- update(lm4, cor=corARMA(p=5, q=1,form=~1|Pond))
M17  <- update(lm4, cor=corARMA(p=1, q=2,form=~1|Pond))
M18  <- update(lm4, cor=corARMA(p=2,q=2, form=~1|Pond))
M19  <- update(lm4, cor=corARMA(p=3,q=2, form=~1|Pond))
#M20 <- update(lm4, cor=corARMA(p=4,q=2, form=~1|Pond))
M21 <- update(lm4, cor=corARMA(p=5,q=2, form=~1|Pond))
M22 <- update(lm4, cor=corARMA(p=1,q=3, form=~1|Pond))
M23 <- update(lm4, cor=corARMA(p=2,q=3, form=~1|Pond))
M24 <- update(lm4, cor=corARMA(p=3,q=3, form=~1|Pond))
M25 <- update(lm4, cor=corARMA(p=4,q=3, form=~1|Pond))
M26 <- update(lm4, cor=corARMA(p=5,q=3, form=~1|Pond))
M27 <- update(lm4, cor=corARMA(p=1,q=4, form=~1|Pond))
M28 <- update(lm4, cor=corARMA(p=2,q=4, form=~1|Pond))
#M29 <- update(lm4, cor=corARMA(p=3,q=4, form=~1|Pond))
M30 <- update(lm4, cor=corARMA(p=4,q=4, form=~1|Pond))
M31 <- update(lm4, cor=corARMA(p=5,q=4, form=~1|Pond))
M32 <- update(lm4, cor=corARMA(p=1,q=5, form=~1|Pond))
#M33 <- update(lm4, cor=corARMA(p=2,q=5, form=~1|Pond))
#M34 <- update(lm4, cor=corARMA(p=3,q=5, form=~1|Pond))
#M35 <- update(lm4, cor=corARMA(p=4,q=5, form=~1|Pond))
#M36 <- update(lm4, cor=corARMA(p=5,q=5, form=~1|Pond))

AIC(lm1,M2,M3,M4, M5, M6,M7,M8,M9, M10,
    M11, M12, M13, M14,M15, M16, M17,M18,M19, 
    M21,M22, M23, M24, M25,M26, M27,  M28, M30,
    M31,M32)


M2 <- gls(Algae ~ Pond*Cat_Y+
            Pond*DaysFromBeg,
          weights=varComb(varIdent(form=~1|Pond),
                          varIdent(form=~1|Cat_Y)),
          #random=~1|newdate, 
          method="REML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)

anova(M2)
qqnorm(M2, abline = c(0,1))
plot(fitted(M2), Resid.gls)
Resid.gls<-resid(M2, type="normalized")
boxplot(Resid.gls~Pond, data=dataAgg)
boxplot(Resid.gls~Treat, data=dataAgg)
plot(Resid.gls~DaysFromBeg, data=dataAgg)
acf(Resid.gls, lag.max=20,main = "")
pacf(Resid.gls, lag.max=20,main = "")


#select factors
M2 <- gls(Algae ~ Pond*Cat_Y+
            Pond*DaysFromBeg,
          weights=varComb(varIdent(form=~1|Pond),
                          varIdent(form=~1|Cat_Y)),
          #random=~1|newdate, 
          method="ML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)



m1<-update(M2,.~.-Pond:Cat_Y)
m2<-update(M2,.~.-DaysFromBeg:Pond)

anova(M2,m1)
anova(M2,m2)


M2 <- gls(Algae ~ Pond*Cat_Y+
            DaysFromBeg,
          weights=varComb(varIdent(form=~1|Pond),
                          varPower(form=~DaysFromBeg)),
          #random=~1|newdate, 
          method="ML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)



m1<-update(M2,.~.-Cat_Y:Pond)
anova(M2,m1)

M2 <- gls(Algae ~ Pond+Cat_Y+
            DaysFromBeg,
          weights=varComb(varIdent(form=~1|Pond),
                          varPower(form=~DaysFromBeg)),
          #random=~1|newdate, 
          method="ML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)

m3<-update(M2,.~.- DaysFromBeg)
m4<-update(M2,.~.-Pond)
m5<-update(M2,.~.-Cat_Y)
anova(M2,m3)
anova(M2,m4)
anova(M2,m5)


M2 <- gls(Algae ~ Pond+
            DaysFromBeg+Cat_Y ,
          #random=~1|newdate, 
          method="REML",
          data=dataAgg,
          cor=corARMA(p=1,q=0, form=~1|Pond),
          na.action = na.omit)


plot(M2)
plot(fitted(M2), Resid.lme)
Resid.lme<-resid(M2, type="normalized")
boxplot(Resid.lme~Pond, data=dataAgg)
boxplot(Resid.lme~Treat, data=dataAgg)
plot(Resid.lme~DaysFromBeg, data=dataAgg)

summary(M2)
anova(M2)
r.squaredGLMM(M2)

Fit<- fitted(M2)
names(Fit)

fit1<-ggplot(dataAgg, aes(DaysFromBeg, Fit))+
  geom_point(aes(DaysFromBeg, Algae))+
  geom_smooth(aes(col="gray"),col="red",method = "lm", lty=1)
#geom_smooth(aes(DaysFromBeg, Abund, col="gray"),method="loess", span=0.3,formula=y~x, se=FALSE, col="blue", lty=1)
fit1+facet_wrap(~Pond, ncol=2, nrow=6)

fit1<-ggplot(dataAgg, aes(DaysFromBeg, Algae))+
  geom_boxplot(aes(x=Pond))



########gam

g2<-gamm(Cond~ s(DaysFromBeg, by=as.numeric(Pond=="AP1"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP2"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP3"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP4"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP5"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP6"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP7"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP8"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP9"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP10"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP11"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP12"))+
           Pond,
         method="REML",
         #correlation=corARMA(form=~DaysFromBeg|Pond, p=1,q=1),
         data=SubFactorsAgg,
         na.action=na.omit)


g3<-gamm(Cond~   s(DaysFromBeg, by=as.numeric(Pond=="AP1"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP2"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP3"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP4"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP5"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP6"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP7"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP8"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP9"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP10"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP11"))+
           s(DaysFromBeg, by=as.numeric(Pond=="AP12"))+
           Pond,
         correlation=corARMA(form=~DaysFromBeg|Pond, p=1,q=1),
         data=SubFactorsAgg,
         na.action=na.omit)


AIC(g, g1, g2, g3)


plot(g1$gam)
plot(g,  rug=TRUE, pages=3,residuals = TRUE,  pch = 1, cex = 1,
     shade = TRUE, shade.col = "lightblue",all.terms = TRUE,
     seWithMean = TRUE, shift = coef(m)[1])

vis.gam(g, theta =120)
summary(g3)
anova(g1$lme)
anova(g1$gam)
par(mfrow=c(2,2), mar=c(2,2,2,2))
gam.check(g)

E2<-resid(g)
Eall<-vector(length=length(SubFactorsAgg))
Eall[]<-NA
I1<-!is.na(SubFactorsAgg)
Eall[I1]<-E2
library(lattice)
xyplot(Eall~SubFactorsAgg$DaysFromBeg|SubFactorsAgg$Pond, col=1)

E2<-resid(g4)
Eall<-vector(length=length(SubFactorsAgg))
Eall[]<-NA
I1<-!is.na(SubFactorsAgg)
Eall[I1]<-E2
library(lattice)
xyplot(Eall~SubFactorsAgg$DaysFromBeg|SubFactorsAgg$Pond, col=1)


E1<-Eall[data$Pond=="AP1"]
E2<-Eall[data$Pond=="AP2"]
E3<-Eall[data$Pond=="AP3"]
E4<-Eall[data$Pond=="AP4"]
E5<-Eall[data$Pond=="AP5"]
E6<-Eall[data$Pond=="AP6"]
E7<-Eall[data$Pond=="AP7"]
E8<-Eall[data$Pond=="AP8"]
E9<-Eall[data$Pond=="AP9"]
E10<-Eall[data$Pond=="AP10"]
E11<-Eall[data$Pond=="AP11"]
E12<-Eall[data$Pond=="AP12"]


par(mfrow=c(4,3))
acf(E1, na.action=na.pass)
acf(E2, na.action=na.pass)
acf(E3, na.action=na.pass)
acf(E4, na.action=na.pass)
acf(E5, na.action=na.pass)
acf(E6, na.action=na.pass)
acf(E7, na.action=na.pass)
acf(E8, na.action=na.pass)
acf(E9, na.action=na.pass)
acf(E10, na.action=na.pass)
acf(E11, na.action=na.pass)
acf(E12, na.action=na.pass)


D<-cbind(E1,E2, E3, E4, E5, E6,
         E7, E8, E9, E10, E11, E12)

cor<-cor(D, use="pairwise.complete.obs")
