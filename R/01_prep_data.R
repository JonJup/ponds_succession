### -------------------------- ###
### --- Prepare Pons data  --- ### 
### -------------------------- ###

# --------------- #
# date:  10.03.21
# files in 
	#-> ENV_MACRO_08_03_21.csv
# files out
	#<- prep_data.rds
# Ponds succession EERES  
# Inspect and prepare data for futher analysis  
# --------------- #

# setup -----------------------------------------------
source("R/setup.R")

# load data -------------------------------------------
data = fread("data/ENV_MACRO_08_03_21.csv")

# prepare data  -----------------------------------------------------------
names(data)[which(names(data) == "Algae")] <- "SubVeg"
data[, date := lubridate::dmy(Date)]
setorderv(data, "date")

#convert in factor (from character)
data$sampler_ID   %<>% as.factor()
data$sampler_ID2  %<>% as.factor()
data$SeasonUnique %<>% as.factor()
data$Year         %<>% as.factor()
data$Date2        %<>% as.factor()
data$Date3        %<>% as.factor()
data$Season       %<>% as.factor()
data$Pond         %<>% as.factor()
data$Treat        %<>% as.factor()

#re-order ponds from 1 to 12
data$Pond <- factor(data$Pond, levels = c("P1", "P2", "P3", "P4",
		  "P5", "P6", "P7", "P8", 
		  "P9", "P10", "P11", "P12"))
#re-order seasons over the years
data$Season_Year  %<>% ordered(levels = c("W_17", "SP_18", "SU_18", "A_18", "W_18", "SP_19", "SU_19"))
data$Season_Year2 %<>% ordered(levels = c("Win1", "Sp1", "Su1", "Au1", "Win2", "Sp2", "Su2"))

###removed variables with too many NAs (not representative) and variables not related to water-environment
datared <-
	data[, c(
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
	) := NULL
	]

saveRDS(datared, "data/prep_data.rds")
