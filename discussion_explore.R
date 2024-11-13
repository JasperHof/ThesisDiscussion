
library(readxl)
library(survminer)
library(survival)

####################
### introduction - sample sizes?
###################

data <- read.csv2("Q:/PL Sita Vermeulen/P Recurrent events/Recurrences_UroLife_NBCS.csv")
data_nbcs  <- read.csv2("Q:/PL Sita Vermeulen/P Recurrent events/Old_NBCS.csv")

length(unique(data$key_nkr_UroLife[!is.na(data$key_nkr_UroLife)]))
length(unique(data_nbcs$GenoID[!is.na(data_nbcs$GenoID)])) + length(unique(data$key_nkr_NBCS[!is.na(data$key_nkr_NBCS)]))

####################
# How are patients with recurrent LR tumours treated?
####################

library(survival)
library(survminer)
library(ggplot2)
library(grid)
library(gridExtra)

### Prepare data

data <- read.csv2("Q:/PL Sita Vermeulen/P Recurrent events/Recurrences_UroLife_NBCS.csv")
if(any(data$ID_UroLife %in% c('UL117034','UL105050','UL104144'))) data = data[-which(data$ID_UroLife %in% c('UL117034','UL105050','UL104144')),]
if(any(data$T_N == 1 & data$Diag_time == 0)) data = data[-which(data$T_N == 1 & data$Diag_time == 0),]
data$EAU[which(data$EAU_old == 'MIBC')] = 'MIBC'
data$EAU[which(data$ID_NBCS == 6221 & data$Diag_time == 0)] = 'High Risk'
library(readxl)
bron <- read_xlsx('Q:/PL Sita Vermeulen/P Recurrent events/Cleaning/Cleaning Nieuwe Urolife/Urolife_20220719.xlsx')
data$lymp = bron$lymfangio[match(data$key_eid_UroLife, bron$key_eid)]; data$lymp[is.na(data$lymp)] = 9
data$diff = as.numeric(data$Morfology %in% c(8041, 8045, 8122, 8131, 8480) | 
                         data$lymp == 1)
data$EAU[data$diff == 1 & data$Diag_time == 0] = 'Very High Risk'
data$T_N[which(data$T_N %in% c('','X'))] = 0
data$T_M[which(data$T_M %in% c('','X'))] = 0
data$T_N_next[which(data$T_N_next %in% c('','X'))] = 0
data$T_M_next[which(data$T_M_next %in% c('','X'))] = 0
data = data[-which(data$EAU == 'MIBC'),]
if(any(data$T_N == 1 & data$Diag_time == 0)) data = data[-which(data$T_N == 1 & data$Diag_time == 0),]
data$subject = paste0(data$ID_UroLife, data$ID_NBCS)
utuc_sync_subs = c("6323","UL112053NA","6018","UL120017NA","UL119052NA","UL103069NA","UL105095NA","6445","6209","6493","6017","6502","UL104042NA","6045","UL107005NA",
                   "UL120024NA","6292","UL123040NA","UL107017NA","UL113008NA","UL105045NA","6426","6399","UL102049NA","UL109025NA","UL118062NA","UL102048NA","UL101029NA","6070","UL113067NA"
                   ,"UL102025NA","UL109004NA","6294","6056","6345","UL102086NA")
data = data[!data$subject %in% utuc_sync_subs,]

#####



