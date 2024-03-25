#Author: Ashleigh Macdonald

#set working directory to the file currently shown in the source pane (i.e. this R script)
#This line of code can be found on stack overflow
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#load required packages
library("tidyverse")
library("dplyr")
library("ggplot2")
library("cowplot")

#From the experiment, will receive two excel files - one containing data before treatment application, and one containing data after treatment application
#In each of these files, there is data collected at the 405 nm excitation wavelength, and data collected at the 480 nm excitation wavelength.
#Manual changes need to first be made to these files which are described below.
#These manual changes were made so that one '405' file and one '480' file were generated - each file would contain either the data collected at the 405 or 480 nm excitation wavelength, that was collected before and after treatment application.

#In each file change the time to time + (time after midnight). This required to be done to merge the two different datasets that were generated before and after treatment application.

#Copy '405' data from the 'before treatment' spreadsheet, including time plus midnight, temp and row labels - paste into new excel file using paste special transpose and tick values box.
#Repeat this with the '405' data from the 'after treatment' spreadsheet - paste into the same excel file below the 'before treatment' data, again using the paste special transpose function.
#Delete temperature data
#Change the heading of 'Time' to 'time'
#Save new dataset as a csv file

#Repeat the above 4 steps for the '480' data.

#####################################################################################
#Initial data processing
#####################################################################################
#Load in relevant datasets
data_405 = read_csv("5.3.23_H2O21_405.csv")
data_480 = read_csv("5.3.23_H2O21_480.csv")

#Convert to long formats
data_405_long = data_405 %>% 
  mutate(time = time - 59152) %>% #sets time to zero (this number will change for each plate that is run)
  mutate(time=round(time/60,2)) %>% #convert time to mins
  gather(sample,gfp405,-time) %>%
  select(sample,time,gfp405) %>%
  extract(sample,"row","([[:alpha:]]+)",remove=FALSE) %>%
  extract(sample,"column","([[:digit:]]+)",remove=FALSE) %>%
  mutate_at("column",funs(as.numeric))

data_480_long = data_480 %>% 
  mutate(time = time - 59152) %>%
  mutate(time=round(time/60,2)) %>% 
  gather(sample,gfp480,-time) %>% 
  select(sample,time,gfp480) %>% 
  extract(sample,"row","([[:alpha:]]+)",remove=FALSE) %>% 
  extract(sample,"column","([[:digit:]]+)",remove=FALSE) %>%
  mutate_at("column",funs(as.numeric))

#Check each measurement overtime
ggplot(data_405_long, aes(x= time))+
  facet_grid(row~column, scale = 'fixed')+
  geom_line(aes(y=gfp405), col = "black")

ggplot(data_480_long, aes(x= time))+
  facet_grid(row~column, scale = 'fixed')+
  geom_line(aes(y=gfp480), col = "black")

#######################################################################
#Make the sample index
#######################################################################
#Label wells by treatment (using an index)
sample <- c('A1', 'B1', 'C1', 'D1', 'E1','F1', 'A2', 'B2', 'C2', 'D2', 'E2','F2','A3','B3','C3','D3','E3','F3','A4','B4','C4','D4','E4','F4','A5','B5','C5','D5','E5','F5','A6','B6','C6','D6','E6','F6')
print(sample)

#Now make a new data frame containing all sample info
sample_pair_ID = c(rep(c('1','2','3','4','5','6'),2), rep(c('7','8','9','10','11','12'),2), rep(c('13','14','15','16','17','18'),2))
print(sample_pair_ID)
plasmid = c(rep(c('roGFP', 'pDSK'), each = 6), rep(c('roGFP', 'pDSK'), each = 6), rep(c('roGFP', 'pDSK'), each = 6))
print(plasmid)
strain = c(rep(c('WT'), 12), rep(c('KO'), 12), rep(c('KO_katG'), 12))
print(strain)
H2O2 = c(rep(c('NA','100','0','2.5','1','0.01'),6))
print(H2O2)
DTT = c(rep(c('5','NA','NA','NA','NA','NA'),6))
print(DTT)

#Add all sample info into an index
data_index = data.frame(sample = sample, strain = strain, sample_pair_ID = sample_pair_ID, plasmid = plasmid, DTT = DTT, H2O2 = H2O2)

#Join the sample info index to the sample data
data_final_data = data_405_long %>%
  left_join(data_index) %>% #add in the data index
  left_join(data_480_long) %>% #add in the 480 data
  mutate(H2O2 = factor(H2O2, levels =c(0,0.01,1,2.5,100))) %>%
  mutate(sample_pair_ID = factor(sample_pair_ID, levels =c(1:18)))

###################################################################################
#Complete roGFP calculations
###################################################################################

#Blank data (using pDSK control)
#Do roGFP405 - pDSK405 and roGFP480 - pDSK480 at each time point and at each condition
#Do roGFP405 blanked /roGFP480 blanked for R (ratio) value
#Normalise R to the max oxidised and max reduced values
#Need to normalise R for 100 mM H2O2 to 0.1
#Need to normalise R for DTT to 1

data_final_processed = data_final_data %>%
  group_by(time, sample_pair_ID) %>%
  mutate(blanked_roGFP405 = gfp405[which(plasmid == "roGFP")] - gfp405[which(plasmid == "pDSK")]) %>%
  mutate(blanked_roGFP480 = gfp480[which(plasmid == "roGFP")] - gfp480[which(plasmid == "pDSK")]) %>%
  mutate(R = blanked_roGFP405/blanked_roGFP480) %>%
  ungroup() %>%
  group_by(time, strain) %>%
  mutate(normalised_R = 0.9*((R - R[which(DTT == 5)])/(R[which(H2O2 == 100)] - R[which(DTT == 5)]))+0.1)

###############################################################
#Plot the charts
###############################################################

#Arrow position will change for each dataset = point at which treatments were added.
#Arrow position = 60311 (first after midnight time for after on raw data) - 59152 (was previously used to zero all data)
#60311 - 59152 = 1159. 1159/60 = 19.32

KO_chart_data <- filter(data_final_processed, !is.na(H2O2))
normalised_KO_chart_data <- filter(data_final_processed, !is.na(H2O2)&H2O2 != 100, time >19.3)

#Raw R KO
ggplot(data = KO_chart_data, aes(x= time, y= R, col=H2O2))+
  facet_grid(.~strain)+
  geom_line(size = 2)+
  xlab("Time (mins)")+
  ylab("405/480 ratio")+
  annotate("segment", x = 19.32, xend = 19.32, y = 0.4, yend = 0.35, arrow = arrow(),size=3) + 
  theme_grey() +
  theme(text = element_text(size = 45)) +
  theme(legend.position="right",axis.text=element_text(size=45),axis.title=element_text(size=45),legend.text = element_text(size=45),legend.title = element_text(size=45), legend.key.size=unit(2.5,"cm"))

ggsave2("5.3.23_catalase_KOs_rawR.png", width=20, height = 20, path = "C:/Users/Ashleigh/OneDrive - Nexus365/Bioreporter_Tests/roGFP_Constructs/roGFP_robust_results_knockouts_H2O2/Charts")

#Normalised R KO
ggplot(data = normalised_KO_chart_data, aes(x= time, y= normalised_R, col=H2O2))+
  facet_grid(.~strain)+
  geom_line(size = 2)+
  xlab("Time (mins)")+
  ylab("405/480 ratio")+
  theme_grey() +
  theme(text = element_text(size = 45)) +
  ylim(c(-0.4,1.2))+
  theme(legend.position="right",axis.text=element_text(size=45),axis.title=element_text(size=45),legend.text = element_text(size=45),legend.title = element_text(size=45), legend.key.size=unit(2.5,"cm"))

ggsave2("5.3.23_catalase_KOs_normalisedR.png", width=20, height = 20, path = "C:/Users/Ashleigh/OneDrive - Nexus365/Bioreporter_Tests/roGFP_Constructs/roGFP_robust_results_knockouts_H2O2/Charts")
