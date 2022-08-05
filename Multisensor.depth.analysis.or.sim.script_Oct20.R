#Written by Danielle Harris, updated Oct 2020 (a small code typo was fixed from the Apr 2016 version)

######################################################
#AIM: 
#This code that will run various acoustic survey scenarios either as simulation or analysis, given certain parameters
#This code has been specifically designed for running analyses with data from several different sites that have different deployment depths and other properties
######################################################
#INPUTS:
#The following arguments will be required for either simulation or analysis:

#Type = fixed or towed survey = this will determine which functions are called.

#Data = "bearing.elev" (elevation angles)  or "bearing.inc" (incidence angles) or "range" (slant ranges) = this will determine which functions are called.

#Dist = depth distribution, "beta" or "uniform"

#Params1 - Params 4 = four separate vectors for maxdetrad (maximum slant range given scenario geometry), maxdepth, maxbearing, and linelength (for towed surveys only) values

#NOTE: each vector will be the length of the number of different monitoring sites 

#NOTE: maxdepth for the fixed scenarios is the elevation from the seafloor

#NOTE: maxbearing for elevation angle - if an instrument can monitor the full circle quadrant, its maxbearing elevation angle is 0 i.e., it can monitor horizontally.  This is converted to pi/2 in the code, which works in incidence angles.
										
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
#Params and distparams will be vectors - all distances are given in m.  Distparams will include characters so uses as.numeric too.
#Iterations will be numeric

##########################################################
#OUTPUTS - in the working directory

#Simulations will produce:
#(1) A text file with useful information and results from the simulation
#(2) An R workspace containing a single object: results.list, which contains detailed results from the simulation
#(3) A histogram of the estimated values of the detection function scale parameter - the scale parameter is the same for each sensor, so only one plot is required.
#(4) Figures of each scenario using data from the last iteration.  This will give an example of the geometry.  Purple dots show the animals that were available for detection (once max. bearing and max. detection radius were considered).  Black solid dots show the detected animals.  The vertical black lines denote the maximum horizontal distance that could be surveyed out to, given any restrictions of the maximum bearing.

#Analyses will produce:
#(1) A text file with useful information and results from the simulation.  Note that the sigma and probability estimates are printed with more decimal places at the end of the text file.
#(2) An R workspace containing a single object: results.list, which contains detailed results from the simulation
#(3) A figure of the estimated detection function.

###########################################################################
#STEP 0 - Source all functions that I need
#NB: make sure that these files are all in the same working directory 
##########################################################################

#load likelihood functions:
source("Depth.LH.functions.vs3.R")

#load distribution functions:
source("Depth.distributions.vs1.R")

#load the optimization functions
source("Multisensor.optimization.vs1.R")

#load the depth.simulation function
source("Multisensor.depth.simulation.vs2.R")

#load the depth.analysis function
source("Multisensor.depth.analysis.vs2.R")

#########################################################################
#OBS-BASED ANALYSIS IN THE PAPER 

#OBS dataset information about depth:
#OBS04: 1993
#OBS10: 2067
#OBS16: 2069

#The critical indidence angle is asin(1.5/1.8)

#These data are needed for working out (a) the maximum detection range (slant range) of the scenarios and (b) the maximum horizontal radius (critical range).

###################################
#STEP 1 - Set up data for analysis
###################################
datain<-read.table("OBSdata_3OBS.txt",header=TRUE)

#convert incident angles to radians
inc.angle.all<-(pi/180)*datain$A.Inc

#split into OBSs

obs<-list(
  "OBS04"=inc.angle.all[datain$OBS=="4"],
  "OBS10"=inc.angle.all[datain$OBS=="10"],
  "OBS16"=inc.angle.all[datain$OBS=="16"])

#draw histograms used in paper

#Inc angles in degrees
hist1<-hist(unlist(obs)/(pi/180),main="",xlab="Incidence angle (degrees)", xlim=c(0,60),ylim=c(0,300))

#Estimated perpendicular range (assuming whale is at the surface)
range.all<-datain$Range

obs.range<-list(
  "OBS04"=range.all[datain$OBS=="4"],
  "OBS10"=range.all[datain$OBS=="10"],
  "OBS16"=range.all[datain$OBS=="16"])

hist2<-hist(unlist(obs.range),main="",xlab="Perpendicular range (m)", ylim=c(0,300),xlim=c(0,3000))

#For the three OBS of interest - their maxdetrads (maximum slant range) are:

#for OBS04
mdrOBS04<-1993/cos(asin(1.5/1.7))
#4235.125

#for OBS10
mdrOBS10<-2067/cos(asin(1.5/1.7))
#4392.375

#for OBS16
mdrOBS16<-2069/cos(asin(1.5/1.7))
#4396.625

###################################
#STEP 2 - Run analysis
###################################

type = "fixed"
data = "bearing.inc"
dist = "beta"
params1<-c(mdrOBS04,mdrOBS10,mdrOBS16)
params2<-c(1993,2067,2069)
params3<-c(asin(1.5/1.7))
params4<-c(0,0,0)
distparams<-c("NULL","maxdepth",100,10,5) #here, the animals will be distributed in the top 100 m of the water column, regardless of instrument depth
det.func.params<-c(0,10000)
label<-"3OBS_analysis"

multisensor.depth.analysis.vs2(type, data, dist, params1, params2, params3, params4, distparams,det.func.params,obs,label)

########################################
#STEP 2 - Estimate density and abundance
########################################

#calculate the horizontal critical ranges for each OBS
maxxy_OBS04<-tan(asin(1.5/1.7))*1.993
maxxy_OBS10<-tan(asin(1.5/1.7))*2.067
maxxy_OBS16<-tan(asin(1.5/1.7))*2.069

#need to weight by area monitored - this uses the horizontal critical ranges
(((156/5.064626e-02 )/(24*pi*maxxy_OBS04^2)*(pi*maxxy_OBS04^2))+((975/4.248443e-02 )/(24*pi*maxxy_OBS10^2)*(pi*maxxy_OBS10^2))+((196/4.228191e-02)/(24*pi*maxxy_OBS16^2)*(pi*maxxy_OBS16^2)))/((pi*maxxy_OBS04^2)+(pi*maxxy_OBS10^2)+(pi*maxxy_OBS16^2))
#[1] 9.236253

#############################################
#STEP 3 - Run bootstrap for confidence limits
#############################################

#The analysis function will return sigma and detection probability values, which can be saved.  The number of observations in each sample will also be saved. Abundance and density can then be calculated.

#To save running the bootstrap, the workspace is provided that contains the bootstrap results - this can be read in (skip from here to the load("OBS_bootstrap.RData") code, currently commented out) and used to explore the density and abundance calculations below

#Note: the bootstrap resamples the OBSs with replacement - if 3 x OBS04 are chosen, then the analysis does not converge.  Therefore, need to run more than 999 iterations so that non-converged results can be omitted.

OBS_bootstrap<-array(NA,c(1200,7))  

for (i in 1:1200){
  
  obs.temp<-sample(obs,3,replace=T)
  params1<-c()
  params2<-c()
  
  OBS_bootstrap[i,1]<-as.numeric(summary(obs.temp)[1]) #no. obs in OBS 1
  OBS_bootstrap[i,2]<-as.numeric(summary(obs.temp)[2]) #no. obs in OBS 2
  OBS_bootstrap[i,3]<-as.numeric(summary(obs.temp)[3]) #no. obs in OBS 3
  
  #set maxdetrads and maxdepths depending on which OBSs have been selected.  Use the number of observations in each OBS to identify which OBS was selected.
  
  if (OBS_bootstrap[i,1]==156){params1[1]<-mdrOBS04}
  if (OBS_bootstrap[i,1]==975){params1[1]<-mdrOBS10}
  if (OBS_bootstrap[i,1]==196){params1[1]<-mdrOBS16}
  if (OBS_bootstrap[i,2]==156){params1[2]<-mdrOBS04}
  if (OBS_bootstrap[i,2]==975){params1[2]<-mdrOBS10}
  if (OBS_bootstrap[i,2]==196){params1[2]<-mdrOBS16}
  if (OBS_bootstrap[i,3]==156){params1[3]<-mdrOBS04}
  if (OBS_bootstrap[i,3]==975){params1[3]<-mdrOBS10}
  if (OBS_bootstrap[i,3]==196){params1[3]<-mdrOBS16}
  
  if (OBS_bootstrap[i,1]==156){params2[1]<-1993}
  if (OBS_bootstrap[i,1]==975){params2[1]<-2067}
  if (OBS_bootstrap[i,1]==196){params2[1]<-2069}
  if (OBS_bootstrap[i,2]==156){params2[2]<-1993}
  if (OBS_bootstrap[i,2]==975){params2[2]<-2067}
  if (OBS_bootstrap[i,2]==196){params2[2]<-2069}
  if (OBS_bootstrap[i,3]==156){params2[3]<-1993}
  if (OBS_bootstrap[i,3]==975){params2[3]<-2067}
  if (OBS_bootstrap[i,3]==196){params2[3]<-2069}
  
  OBS_bootstrap[i,4:7]<-multisensor.depth.analysis.vs2("fixed","bearing.inc","beta",params1,params2,c(asin(1.5/1.7)),c(0,0,0),c("NULL","maxdepth",100,10,5),c(0,10000),obs.temp,"bs_3OBS_stratsamp_Sept18")
  
  print(i)}

save.image(file="OBS_bootstrap.RData")

#here the pre-saved workspace can be read in
#load("OBS_bootstrap.RData")

#take away the nonconverged results
find.non.converge<-rowSums(OBS_bootstrap[,1:3])
OBS_bootstrap_converge<-OBS_bootstrap[!find.non.converge==606,]
OBS_bootstrap_converge_999<-OBS_bootstrap_converge[1:999,]

OBS_bootstrap_sigma<-OBS_bootstrap_converge_999[,4]
OBS_bootstrap_sigma_sort<-sort(OBS_bootstrap_sigma)
OBS_bootstrap_sigma_sort[25]
#[1] 769.6891
OBS_bootstrap_sigma_sort[975]
#[1] 2198.783

#can use the sigma UCI and LCI to estimate the CIs of each pdet
type = "fixed"
data = "bearing.inc"
dist = "beta"
params1<-c(mdrOBS04,mdrOBS10,mdrOBS16)
params2<-c(1993,2067,2069)
params3<-c(asin(1.5/1.7))
params4<-c(0,0,0)
distparams<-c("NULL","maxdepth",100,10,5)
det.func.params<-c(0,10000)

#rerunning using Multisensor.depth.analysis code (run Steps 1-4 in Multisensor.depth.analysis.vs2.R, with sigmaest in Step 4 set to LCI and UCL values for sigma)
#LCI
sigmaest<-769.6891
#results.pdet
#[1] 0.003321591 0.002406643 0.002385557

#UCI
sigmaest<-2198.783
#results.pdet
#[1] 0.3502738 0.3264349 0.3258053

#################################################################
#STEP 4 - Use the bootstrap results for call density CIs
#################################################################

#now need to know which OBSs were picked so that density can be estimated and weighted by monitored area

max_xys<-c(1.993,2.067,2.069)*tan(asin(1.5/1.7)) 
#[1] 3.736875 3.875625 3.879375

#set up the coorect critical ranges for each bootstrap - this can be done by using the number of observations to identify which OBS was selected.
OBS_max_xy<-OBS_bootstrap_converge_999[,1:3]
OBS_max_xy<-replace(OBS_max_xy,OBS_max_xy==156,max_xys[1])
OBS_max_xy<-replace(OBS_max_xy,OBS_max_xy==975,max_xys[2])
OBS_max_xy<-replace(OBS_max_xy,OBS_max_xy==196,max_xys[3])

OBS_bootstrap_density<-(((OBS_bootstrap_converge_999[,1]/OBS_bootstrap_converge_999[,5])/(24*pi*OBS_max_xy[,1]^2)*(pi*OBS_max_xy[,1]^2))+((OBS_bootstrap_converge_999[,2]/OBS_bootstrap_converge_999[,6])/(24*pi*OBS_max_xy[,2]^2)*(pi*OBS_max_xy[,2]^2))+((OBS_bootstrap_converge_999[,3]/OBS_bootstrap_converge_999[,7])/(24*pi*OBS_max_xy[,3]^2)*(pi*OBS_max_xy[,3]^2)))/((pi*OBS_max_xy[,1]^2)+(pi*OBS_max_xy[,2]^2)+(pi*OBS_max_xy[,3]^2))

OBS_bootstrap_density_sort<-sort(OBS_bootstrap_density)
OBS_bootstrap_density_sort[25]
#[1] 0.4229983
OBS_bootstrap_density_sort[975]
#[1] 72.40714

#########################################
#OBS SIMULATION (not presented in paper).
##########################################
#Using the estimated sigma from the analysis.

type = "fixed"
data = "bearing.inc"
dist = "beta"
params1<-c(mdrOBS04,mdrOBS10,mdrOBS16)
params2<-c(1993,2067,2069)
params3<-c(asin(1.5/1.7))
params4<-c(0,0,0)
distparams<-c("NULL","maxdepth",100,10,5)  
det.func.params<-c(1184,0,10000)
popn<-c("density",100)
iterations<-5                              
#5 iterations set here but results of 999 simulations are given in a pre-saved workspace (below)
label<-"OBS_3instr_sim"

multisensor.depth.simulation.vs2(type, data, dist, params1, params2, params3, params4, popn, iterations, distparams,det.func.params,label)

#load("OBS_3instr_sim_999.RData")
