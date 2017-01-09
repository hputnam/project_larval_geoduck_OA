#project-geoduck-oa
#Data from Roberts et al NOAA OA
#last modified 20160320 H Putnam


rm(list=ls()) # removes all prior objects

#Read in required libraries
library("seacarb") #seawater carbonate chemistry
library("reshape") #reshape data
library("plotrix") #functions in tapply
library(ggplot2)
library(gridExtra)

#############################################################
setwd("/Users/hputnam/MyProjects/Geoduck_Epi/project_larval_geoduck_OA/RAnalysis/Data/") #set working directory

#Required Data files
#pH_Calibration_Files/
#SW_Chem_Trial3.csv
#Cell_Counts_Trial3.csv
#Larval_Counts_Trial3.csv

#SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS

##### pH Tris Calibration Curves #####
#For conversion equations for pH from mV to total scale using tris standard

path <-("/Users/hputnam/MyProjects/Geoduck_Epi/project-geoduck-oa/RAnalysis/Data/pH_Calibration_Files/")

#list all the file names in the folder to get only get the csv files
file.names<-list.files(path = path, pattern = "csv$")

pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=4, dimnames=list(file.names,c("Date", "Intercept", "Slope","R2")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  R <- summary(model)$r.squared #extracts the R2
  pH.cals[i,2:3] <- coe #inserts coef in the dataframe
  pH.cals[i,4] <- R #inserts R2 in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}

colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope", "R2")
pH.cals

# read in total alkalinity, temperature, and salinity
SW.chem <- read.csv("SW_Chem_Trial3.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, SW.chem, by="Calib.Date")

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-27.5 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

##### Seacarb Calculations #####

#Calculate CO2 parameters using seacarb
carb.output <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb

carb.output$ALK <- carb.output$ALK*1000000 #convert to µmol kg-1
carb.output$CO2 <- carb.output$CO2*1000000 #convert to µmol kg-1
carb.output$HCO3 <- carb.output$HCO3*1000000 #convert to µmol kg-1
carb.output$CO3 <- carb.output$CO3*1000000 #convert to µmol kg-1
carb.output$DIC <- carb.output$DIC*1000000 #convert to µmol kg-1

carb.output <- cbind(SW.chem$Measure.Date,  SW.chem$Tank,  SW.chem$Treatment, carb.output) #combine the sample information with the seacarb output
colnames(carb.output) <- c("Date",  "Tank",  "Treatment",	"flag",	"Salinity",	"Temperature",	"Pressure",	"pH",	"CO2",	"pCO2",	"fCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat", 	"Calcite.Sat") #Rename columns to describe contents
carb.output <- subset(carb.output, select= c("Date",  "Tank",  "Treatment",	"Salinity",	"Temperature",		"pH",	"CO2",	"pCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat"))

##### Descriptive Statistics #####
#Tanks
mean_pCO2=tapply(carb.output$pCO2, carb.output$Tank, mean)
se_pCO2=tapply(carb.output$pCO2, carb.output$Tank, std.error)
mean_Temp=tapply(carb.output$Temperature, carb.output$Tank, mean)
se_Temp=tapply(carb.output$Temperature, carb.output$Tank, std.error)
mean_Sal=tapply(carb.output$Salinity, carb.output$Tank, mean)
se_Sal=tapply(carb.output$Salinity, carb.output$Tank, std.error)
mean_TA=tapply(carb.output$TA, carb.output$Tank, mean)
se_TA=tapply(carb.output$TA, carb.output$Tank, std.error)
mean_pH=tapply(carb.output$pH, carb.output$Tank, mean)
se_pH=tapply(carb.output$pH, carb.output$Tank, std.error)
mean_DIC=tapply(carb.output$DIC, carb.output$Tank, mean)
se_DIC=tapply(carb.output$DIC, carb.output$Tank, std.error)

#Treatments
gmean_pCO2 <- tapply(carb.output$pCO2, carb.output$Treatment, mean)
gse_pCO2 <- tapply(carb.output$pCO2, carb.output$Treatment, std.error)
gmean_Temp <- tapply(carb.output$Temperature, carb.output$Treatment, mean)
gse_Temp <- tapply(carb.output$Temperature, carb.output$Treatment, std.error)
gmean_Sal <- tapply(carb.output$Salinity, carb.output$Treatment, mean)
gse_Sal <- tapply(carb.output$Salinity, carb.output$Treatment, std.error)
gmean_TA <- tapply(carb.output$TA, carb.output$Treatment, mean)
gse_TA <- tapply(carb.output$TA, carb.output$Treatment, std.error)
gmean_pH <- tapply(carb.output$pH, carb.output$Treatment, mean)
gse_pH <- tapply(carb.output$pH, carb.output$Treatment, std.error)
gmean_DIC <- tapply(carb.output$DIC, carb.output$Treatment, mean)
gse_DIC <- tapply(carb.output$DIC, carb.output$Treatment, std.error)

mean.carb.output <- rbind(gmean_pCO2, gse_pCO2, gmean_pH, gse_pH, gmean_Temp, gse_Temp, gmean_Sal, gse_Sal, gmean_TA, gse_TA, gmean_DIC, gse_DIC)
mean.carb.output <- as.data.frame(mean.carb.output)
row.names(mean.carb.output) <- c("mean pCO2", "SE pCO2", "mean pH", "SE pH", "mean Temperature", "SE Temperature", "mean Salinity", "SE Salinity", "mean Total Alkalinity", "SE Total Alkalinity", "mean DIC", "SE DIC")
mean.carb.output$Variables <- row.names(mean.carb.output)
write.table (mean.carb.output, file="/Users/hputnam/MyProjects/Geoduck_Epi/project-geoduck-oa/RAnalysis/Output/Seawater_chemistry_table_Output_Trial3.csv", sep=",", row.names = FALSE)

##### Cell Counts #####

cell.counts <- read.csv("Cell_Counts_Trial3.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cell.counts$Avg.Cells <- rowMeans(cell.counts[,c("Count1",  "Count2")], na.rm = TRUE) #calculate average of counts
cell.counts$cells.ml <- cell.counts$Avg.Cells/cell.counts$Volume.Counted #calculate density

#Tanks
mean_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, mean, na.rm = TRUE)
se_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, std.error, na.rm = TRUE)

#Treatments
gmean_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, mean, na.rm = TRUE)
gse_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, std.error, na.rm = TRUE)

##### Larval Counts #####
##### Larval Counts #####
larval.counts <- read.csv("Larval_Counts_Trial3.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
larval.counts$Avg.Live <- rowMeans(larval.counts[,c("Live1",  "Live2",  "Live3")], na.rm = TRUE) #calculate average of counts
larval.counts$Avg.Dead <- rowMeans(larval.counts[,c("Dead1",  "Dead2",  "Dead3")], na.rm = TRUE) #calculate average of counts
larval.counts$Live.cells.ml <- larval.counts$Avg.Live/larval.counts$Volume.Counted.ml #calculate density
larval.counts$Dead.cells.ml <- larval.counts$Avg.Dead/larval.counts$Volume.Counted.ml #calculate density
larval.counts$total.live.larvae <- larval.counts$Live.cells.ml*larval.counts$Vol.Tripour
larval.counts$total.dead.larvae <- larval.counts$Dead.cells.ml*larval.counts$Vol.Tripour
larval.counts$per.sur <- ((larval.counts$total.live.larvae/(larval.counts$total.live.larvae[1]))*100)

#calculate mean and se of live larvae
larval.counts.means <- aggregate(total.live.larvae ~ Day*Treatment, data=larval.counts, mean)
larval.counts.ses <- aggregate(total.live.larvae ~ Day*Treatment, data=larval.counts, std.error)
means <- cbind(larval.counts.means, larval.counts.ses$total.live.larvae) #create dataframe
colnames(means) <- c("Day", "Treatment","mean","se") #rename columns

#calculate mean and se of percent survival
larval.survival.means <- aggregate(per.sur ~ Day*Treatment, data=larval.counts, mean)
larval.survival.ses <- aggregate(per.sur ~ Day*Treatment, data=larval.counts, std.error)
sur.means <- cbind(larval.survival.means, larval.survival.ses$per.sur) #create dataframe
colnames(sur.means) <- c("Day", "Treatment","mean","se") #rename columns


#plot survivorship
Fig.Survivorship <- ggplot(sur.means, aes(x=Day, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Treatment), position = position_dodge(width = 0.2), size=2) + #plot points
  scale_shape_manual(values=c(16,17,15,17)) +
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  #scale_x_discrete(limits=c("Time0", "Time1", "Time2", "Time3", "Time4", "Time5"))+
  xlab("Time") + #Label the X Axis
  ylab("% Survivorship") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        panel.border = element_blank(), #Set the border
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position=c(0.8,0.6)) #set legend position
Fig.Survivorship

#####Larval Size#####
size <- read.csv("OA_Size_Larval_Geoduck_Trial3.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
info <- read.csv("Larval_Sample.Info.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#check the variation in the scale between pictures
CV.mean <- aggregate(Scale ~ Tube.ID, data=size, mean, na.rm=TRUE) #calculate mean 
CV.sd <- aggregate(Scale ~ Tube.ID, data=size, sd, na.rm=TRUE) #calculate standard deviation
CVs <- cbind(CV.mean, CV.sd$Scale) # merge mean and sd into dataframe
colnames(CVs) <- c("Tube.ID",  "Mean.Scale",  "SD.Scale") #rename the columns of the dataframe
CVs$CV <- (CVs$SD.Scale/CVs$Mean.Scale)*100 #calculate CV of the scale bar measurements
range(CVs$CV) #display min and max of scale CV within a tube
Grand.CV <- (sd(na.omit(size$Scale))/mean(na.omit(size$Scale)))*100#calcualte the overall CV
Grand.CV #coefficient of variation in scale bar between all samples ~2%

data <- merge(size, info, by="Tube.ID") #merge data and info
data$Ratio <- data$Width/data$Length #calcualte width to length ratio

means <- aggregate(cbind(Length, Width,Area, Ratio)  ~ Treatment*Time*Trial, data=data, mean, na.rm=TRUE) #calculate mean
ses <- aggregate(cbind(Length, Width,Area, Ratio) ~ Treatment*Time*Trial, data=data, std.error, na.rm=TRUE) #calculate standard error
sizes <- cbind(means, ses$Length, ses$Width, ses$Area, ses$Ratio)
colnames(sizes) <- c("Treatment","Time", "Trial", "Length", "Width", "Area", "Ratio", "Length.se", "Width.se", "Area.se", "Ratio.se")

Fig.Length <- ggplot(sizes, aes(x=Time, y=Length, group=Treatment)) + 
  geom_errorbar(aes(ymin=Length-Length.se, ymax=Length+Length.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Treatment), position = position_dodge(width = 0.2), size=2) + #plot points
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Shell Length mm") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        panel.border = element_blank(), #Set the border
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position=c(0.8,0.3)) #set legend position
Fig.Length

Fig.Width <- ggplot(sizes, aes(x=Time, y=Width, group=Treatment)) + 
  geom_errorbar(aes(ymin=Width-Width.se, ymax=Width+Width.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Treatment), position = position_dodge(width = 0.2), size=2) + #plot points
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Shell Width mm") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        panel.border = element_blank(), #Set the border
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position='none') #remove legend background
Fig.Width

Fig.Area <- ggplot(sizes, aes(x=Time, y=Area, group=Treatment)) + 
  geom_errorbar(aes(ymin=Area-Area.se, ymax=Area+Area.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Treatment), position = position_dodge(width = 0.2), size=2) + #plot points
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Shell Area mm2") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        panel.border = element_blank(), #Set the border
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position='none') #remove legend background
Fig.Area

Fig.Ratio <- ggplot(sizes, aes(x=Time, y=Ratio, group=Treatment)) + 
  geom_errorbar(aes(ymin=Ratio-Ratio.se, ymax=Ratio+Ratio.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Treatment), position = position_dodge(width = 0.2), size=2) + #plot points
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Shell Ratio (W:L)") + #Label the Y Axis
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        axis.line.x = element_line(color = 'black'), #Set the axes color
        axis.line.y = element_line(color = 'black'), #Set the axes color
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(), #Set plot legend key
        legend.position='none') #remove legend 
Fig.Ratio

setwd("/Users/hputnam/MyProjects/Geoduck_Epi/project_larval_geoduck_OA/RAnalysis/Output")
Larval.Size <- arrangeGrob(Fig.Length, Fig.Width, Fig.Area, Fig.Survivorship, ncol=2)
ggsave(file="Larval_Size_Trial3.pdf", Larval.Size, width =8.5, height = 11, units = c("in"))
