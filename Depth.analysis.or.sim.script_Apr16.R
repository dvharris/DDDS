#Written by Danielle Harris, updated 20 April 2016

######################################################
#AIM: 
#This code that will run various acoustic survey scenarios either as simulation or analysis, given certain parameters

######################################################
#INPUTS:
#The following arguments will be required for either simulation or analysis:

#Type = fixed or towed survey = this will determine which functions are called.
										
#Data = "bearing.elev" (elevation angles)  or "bearing.inc" (incidence angles) or "range" (slant ranges) = this will determine which functions are called.
										
#Dist = depth distribution, "beta" or "uniform"
										
#Params = maxdetrad (maximum slant range given scenario geometry), maxdepth, maxbearing, linelength (for towed surveys only)

#NOTE: maxdepth for the fixed scenarios is the elevation from the seafloor

#NOTE: maxbearing for elevation angle data - if an instrument can monitor the full circle quadrant, its maxbearing elevation angle is 0 i.e., it can monitor horizontally.  This is converted to pi/2 in the code, which works in incidence angles.
										
#Popn = number or density of animals in the simulated population, name ("abundance" or "density") and value
										
#Iterations = number of iterations of the simulation 
									  
#Distparams = depth distribution parameters.

#NOTE: Be careful of whether the maxdepth is the maximum depth or altitude, this will determine the distribution parameters to use.  

#NOTE: Three or five inputs are given - lower and upper limit and spread (here, spread is defined as the maximum range of depths that animals are distributed) are always required.  However, only two out of three values are required so that the third can be calculated.  For example, the upper limit can be the max depth/elevation so that the lower limit can be altered using the spread if different depths are used (in the multiple instrument scenario).

#NOTE: ll can take a value or "NULL", ul can take a value, "NULL" or "maxdepth", spread can take a value or "NULL"

#NOTE: The code does not have sufficient checks to stop invalid combinations of limits of the depth distribution in relation to the specified maximum depth/elevation.  The code might still run with "bad" combinations of values.  

#NOTE: If a beta distribution is also given, then two shape parameters are also set).
                    
#Det.func.params= value of sigma to simulate from, and range of values to maximise over
                    
#Label is a name that will be appended to the results text files that are saved.

#A note on class types of inputs (for information only):
#Type, Data and Dist will be characters
#Popn will also be treated as a character vector as the first term will be a character, so use as.numeric in code to read value correctly
#Params and distparams will be vectors- all distances are given in m.  Distparams will include characters so use as.numeric too.
#Iterations will be numeric

##########################################################
#OUTPUTS - in the working directory

#Simulations will produce:
#(1) A text file with useful information and results from the simulation
#(2) An R workspace containing a single object: results.list, which contains detailed results from the simulation
#(3) A figure of the density estimates compared to the sigma estimates - this plot can be useful for seeing if there were any optimisation problems, and if extreme density values were estimated.
#(4) A figure of the last simulated population in the simulation, as an example of the survey geometry.  Purple dots show the animals that were available for detection (once max. bearing and max. detection radius were considered).  Black solid dots show the detected animals.  The vertical black lines denote the maximum horizontal distance that could be surveyed out to, given any restrictions of the maximum bearing.

#Analyses will produce:
#(1) A text file with useful information and results from the simulation.  Note that the sigma and probability estimates are printed with more decimal places at the end of the text file.
#(2) An R workspace containing a single object: results.list, which contains detailed results from the simulation
#(3) A figure of the estimated detection function.

###########################################################################
#
#STEP 0 - Set up preliminary values

#load likelihood functions
#NB: make sure that this file is in the same working directory as the session
source("Depth.LH.functions.vs3.R")

#load distribution functions:
#NB: make sure that this file is in the same working directory as the session
source("Depth.distributions.vs1.R")

#load the depth.simulation function
source("Depth.simulation.vs2.R")

#load the depth.simulation function
source("Depth.analysis.vs2.R")
  
#########################################################################
#SIMULATIONS IN PAPER - only elevation bearing results were reported in the paper, but incidence angle code also included here.

#Using the parameters for the paper:
#depth, maxdetrad - 2000m, linelength = 4000m
#maxbearing = pi/2, or 0 for bearing.elev data
#beta shape: 1.5, 70 (and vice versa for fixed**) - shallow divers
#beta shape: 5, 10 (and vice versa for fixed**) - deep divers 

#**Note that for fixed sensors, the shape parameters of the beta distribution must be interchanged to create a comparable depth distribution with a simulation using towed instruments, because the calculations using fixed sensors work with elevation, not depth.

#a reminder of the arguments that the function expects
#depth.simulation.vs1<-function(type, data, dist, params, popn, iterations, distparams,detfuncparams,label)

#Five iterations specifed here - each simulation will take a few minutes to run

depth.simulation.vs2("towed","range","uniform",c(2000,2000,pi/2,4000),c("abundance",2500),5,c(0,"maxdepth","NULL",5,10),c(500,0,10000),"towed.range.uniform")

depth.simulation.vs2("fixed","range","uniform",c(2000,2000,pi/2,2000),c("abundance",2500),5,c(0,"maxdepth","NULL",5,10),c(800,0,10000),"fixed.range.uniform")

depth.simulation.vs2("fixed","bearing.elev","uniform",c(2000,2000,0,2000),c("abundance",2500),5,c(0,"maxdepth","NULL",10,5),c(800,0,10000),"fixed.bearing.elev.uniform") #poor results

depth.simulation.vs2("fixed","bearing.inc","uniform",c(2000,2000,pi/2,2000),c("abundance",2500),5,c(0,"maxdepth","NULL",10,5),c(800,0,10000),"fixed.bearing.inc.uniform") #poor results

depth.simulation.vs2("towed","range","beta",c(2000,2000,pi/2,4000),c("abundance",2000),5,c(0,"maxdepth","NULL",5,10),c(500,0,10000),"towed.range.beta.deep")

depth.simulation.vs2("towed","bearing.elev","beta",c(2000,2000,0,4000),c("abundance",2000),5,c(0,"maxdepth","NULL",5,10),c(500,0,10000),"towed.bear.elev.beta.deep")

depth.simulation.vs2("towed","bearing.inc","beta",c(2000,2000,pi/2,4000),c("abundance",2000),5,c(0,"maxdepth","NULL",5,10),c(500,0,10000),"towed.bear.inc.beta.deep")

depth.simulation.vs2("fixed","range","beta",c(2000,2000,pi/2,2000),c("abundance",4000),5,c(0,"maxdepth","NULL",10,5),c(800,0,10000),"fixed.range.beta.deep")

depth.simulation.vs2("fixed","bearing.elev","beta",c(2000,2000,0,2000),c("abundance",4000),5,c(0,"maxdepth","NULL",10,5),c(800,0,1000),"fixed.bearing.elev.beta.deep")

depth.simulation.vs2("fixed","bearing.inc","beta",c(2000,2000,pi/2,2000),c("abundance",4000),5,c(0,"maxdepth","NULL",10,5),c(800,0,1000),"fixed.bearing.inc.beta.deep")

depth.simulation.vs2("towed","range","beta",c(2000,2000,pi/2,4000),c("abundance",1000),5,c(0,"maxdepth","NULL",10,70),c(500,0,10000),"towed.range.beta.shallow")

depth.simulation.vs2("towed","bearing.elev","beta",c(2000,2000,0,4000),c("abundance",1000),5,c(0,"maxdepth","NULL",10,70),c(500,0,10000),"towed.bear.elev.beta.shallow")

depth.simulation.vs2("towed","bearing.inc","beta",c(2000,2000,pi/2,4000),c("abundance",1000),5,c(0,"maxdepth","NULL",10,70),c(500,0,10000),"towed.bear.inc.beta.shallow")

depth.simulation.vs2("fixed","range","beta",c(2000,2000,pi/2,2000),c("abundance",20000),5,c(0,"maxdepth","NULL",70,10),c(800,0,10000),"fixed.range.beta.shallow")

depth.simulation.vs2("fixed","bearing.elev","beta",c(2000,2000,0,2000),c("abundance",20000),5,c(0,"maxdepth","NULL",70,10),c(800,0,1000),"fixed.bearing.elev.beta.shallow")

depth.simulation.vs2("fixed","bearing.inc","beta",c(2000,2000,pi/2,2000),c("abundance",20000),5,c(0,"maxdepth","NULL",70,10),c(800,0,1000),"fixed.bearing.inc.beta.shallow")

############################################################################
#GOALSII DDDS analysis from the paper

###################################
#STEP 1 - Set up data for analysis
###################################

#read in distance data:
obs.data<-read.csv("GOALSII_data_for_analysis.csv",header=TRUE)

#set up the data for analysis by converting distances to km and truncating at 5 km
obs.km<-as.vector(na.omit(obs.data$Observation.Perp.Distance))
obs.m<-obs.km*1000
obs.m.trunc5km<-obs.m[obs.m<=5000]

#make a histogram of the slant ranges
hist(obs.m,xlab="Slant range (m)",main="",breaks=16,xlim=c(0,8000))

###################################
#STEP 2 - Run analysis
###################################

#treat depth as 3800m
#beta distribution will be s1=5,s2=15

depth.analysis.vs2("towed","range","beta",c(5000,3800,pi/2,3441960),c(0,"maxdepth","NULL",5,15),c(0,10000),obs.m.trunc5km,"GOALSII_analysis")

#The estimate of sigma was 3076.2 and the estimate of the detection probability was 0.65

########################################
#STEP 2 - Estimate density and abundance
########################################

#The density estimation in the analysis functions cannot account for a stratified design, so we must estimate density by hand.

offshore.n<-8     #no. encounters
offshore.k<-7     #no. lines
offshore.L<-1085.65 #line length in km 
offshore.w<-5     #truncation distance in km
offshore.A<-60051 #area in km2

seamount.n<-26    #no. encounters
seamount.k<-25    #no. lines
seamount.L<-2356.310 #line length in km 
seamount.w<-5     #truncation distance in km
seamount.A<-45377 #area in km2

global_slant.p<-0.6493449

#Estimate point estimates of abundance and density

offshore.N.slant<-(offshore.n*offshore.A)/(2*offshore.w*offshore.L*global_slant.p)
#[1] 68.14672
seamount.N.slant<-(seamount.n*seamount.A)/(2*seamount.w*seamount.L*global_slant.p)
#[1] 77.10832 
total.N.slant<-offshore.N.slant+seamount.N.slant
#145.255

offshore.D.slant<-(offshore.n)/(2*offshore.w*offshore.L*global_slant.p)
#0.001134814 animals/km2 = 1.13 animals/1000km2
seamount.D.slant<-(seamount.n)/(2*seamount.w*seamount.L*global_slant.p)
#0.001699282 animals/km2 = 1.70 animals/1000km2 
w.avg.D.slant<-((offshore.A*offshore.D.slant)+(seamount.A*seamount.D.slant))/(offshore.A+seamount.A)
#[1] 0.001377765 animals/km2 = 1.38 animals/1000km2 

#############################################
#STEP 3 - Run bootstrap for confidence limits
#############################################

#The analysis function will return sigma and detection probability values, which can be saved.  The number of observations in each strata will also be saved. The line length of each stratum will also be saved.  Abundance and density can then be calculated.

#To save running the bootstrap, the workspace is provided that contains the bootstrap results - this can be read in (skip from here to the load("GOALSII_999bs_stratsamp.RData") code, currently commented out) and used to explore the density and abundance calculations below

GOALSII_bootstrap<-array(NA,c(999,6))
no.strata<-2
offshore.lines<-as.vector(unique(obs.data$Line.Transect.Label[obs.data$Region.Label=="offshore"]))  
seamount.lines<-as.vector(unique(obs.data$Line.Transect.Label[obs.data$Region.Label=="seamount"]))  

for (i in 1:999){
  obs.m.temp<-numeric()
  obs.line.length.temp<-numeric()
  obs.data.bs<-numeric()
  
  offshore.lines.temp<-as.vector(sample(offshore.lines,length(offshore.lines),replace=TRUE))
  
  for (k in 1:length(offshore.lines.temp)){
    obs.m.bs.temp<-obs.data$Observation.Perp.Distance[obs.data$Line.Transect.Label==offshore.lines.temp[k]]*1000
    obs.m.temp<-c(obs.m.temp,obs.m.bs.temp)
    obs.line.length.temp<-c(obs.line.length.temp,as.vector(unique(obs.data$Line.Transect.Line.Length[obs.data$Line.Transect.Label==offshore.lines.temp[k]])))
  }
  
  #truncate at 5 km
  obs.data.bs<-obs.m.temp[obs.m.temp<=5000]
  GOALSII_bootstrap[i,1]<-sum(na.omit(obs.line.length.temp)) #line length
  GOALSII_bootstrap[i,2]<-length(na.omit(obs.data.bs)) #no. obs
  
  seamount.lines.temp<-as.vector(sample(seamount.lines,length(seamount.lines),replace=TRUE))
  
  obs.m.temp<-numeric()
  obs.line.length.temp<-numeric()
  
  for (l in 1:length(seamount.lines.temp)){
    obs.m.bs.temp<-obs.data$Observation.Perp.Distance[obs.data$Line.Transect.Label==seamount.lines.temp[l]]*1000
    obs.m.temp<-c(obs.m.temp,obs.m.bs.temp)
    obs.line.length.temp<-c(obs.line.length.temp,as.vector(unique(obs.data$Line.Transect.Line.Length[obs.data$Line.Transect.Label==seamount.lines.temp[l]])))
  }
  
  GOALSII_bootstrap[i,3]<-sum(na.omit(obs.line.length.temp))
  GOALSII_bootstrap[i,4]<-length(na.omit(obs.m.temp))
  obs.data.bs<-c(obs.data.bs,obs.m.temp[obs.m.temp<=5000])
  obs.m.trunc<-as.vector(na.omit(obs.data.bs))
  
  if (sum(is.na(GOALSII_bootstrap[i,1:4]))==0){
    GOALSII_bootstrap[i,5:6]<-depth.analysis.vs2("towed","range","beta",c(5000,3800,pi/2,3441960),c(0,"maxdepth","NULL",5,15),c(0,10000),obs.m.trunc,"GOALSII_bootstrap")}
  
  print(i)}

#here the pre-saved workspace can be read in
#load("GOALSII_999bs_stratsamp.RData")

#Get 95% confidence intervals
GOALSII_bootstrap_sigma<-GOALSII_bootstrap[,5]
GOALSII_bootstrap_sigma_sort<-sort(GOALSII_bootstrap_sigma)
GOALSII_bootstrap_sigma_sort[25]
#[1]  1785.639
GOALSII_bootstrap_sigma_sort[975]
#[1] 10000

GOALSII_bootstrap_pdet<-GOALSII_bootstrap[,6]
GOALSII_bootstrap_pdet_sort<-sort(GOALSII_bootstrap_pdet)
GOALSII_bootstrap_pdet_sort[25]
#[1] 0.3810148
GOALSII_bootstrap_pdet_sort[975]
#[1] 0.9364848

#################################################################
#STEP 4 - Use the bootstrap results for density and abundance CIs
#################################################################

#Density for each stratum
offshore.D.slant.bs<-(GOALSII_bootstrap[,2])/(2*offshore.w*GOALSII_bootstrap[,1]*GOALSII_bootstrap_pdet)
seamount.D.slant.bs<-(GOALSII_bootstrap[,4])/(2*seamount.w*GOALSII_bootstrap[,3]*GOALSII_bootstrap_pdet)

#Weighted average density using area as weight
w.avg.D.slant_bs<-((offshore.A*offshore.D.slant.bs)+(seamount.A*seamount.D.slant.bs))/(offshore.A+seamount.A)

#Sort and take 25th and 975th values
w.avg.D.slant_bs_sort<-sort(w.avg.D.slant_bs)
w.avg.D.slant_bs_sort[25]
#[1] 0.0006770153
w.avg.D.slant_bs_sort[975]
#[1]  0.002685653

#Abundance
offshore.N.slant.bs<-(GOALSII_bootstrap[,2]*offshore.A)/(2*offshore.w*GOALSII_bootstrap[,1]*GOALSII_bootstrap_pdet)
seamount.N.slant.bs<-(GOALSII_bootstrap[,4]*seamount.A)/(2*seamount.w*GOALSII_bootstrap[,3]*GOALSII_bootstrap_pdet)

#Total abundance
w.avg.N.slant_bs<-offshore.N.slant.bs+seamount.N.slant.bs

#Sort and take 25th and 975th values
w.avg.N.slant_bs_sort<-sort(w.avg.N.slant_bs)
w.avg.N.slant_bs_sort[25]
#[1] 71.37637
w.avg.N.slant_bs_sort[975]
#[1] 283.1431










