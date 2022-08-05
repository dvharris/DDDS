depth.analysis.vs2<-function(type, data, dist, params, distparams,det.func.params,obs,label){

#NB: det.func.params is only 2 values - range to look over for optimization. 
  ################################################################################
  #STEP 1: Set up the analysis
  
  #Results vectors
  
  if (type=="fixed") {

    results.sig<-numeric()
    results.pdet<-numeric()
    results.est.cyl.max_xy<-numeric()
    results.dens<-numeric()
    results.obs<-numeric() } else 
      
  if (type=="towed") {

    results.sig<-numeric()
    results.pdet<-numeric()
    results.est.rect.max_xy<-numeric()
    results.dens<-numeric()
    results.obs<-numeric() }

  #Enter the parameters of the survey area and simulation 

  #Define the 3D survey area

  #Maximum radial distance of the detection process - this is the maximum slant detection range for either a towed or fixed transect:
  maxdetrad<-params[1]
  #Maximum depth
  #Note that the calculations below work in exactly the same way for towed and fixed scenarios. 
  #The only difference is that depth in the towed scenario is actual depth i.e., 2000m = the seafloor
  #In contrast, in the fixed scenario, depth is actually elevation i.e., 2000m = the sea surface.
  #The metric is actually the distance between the instrument and a boundary, either the sea surface or sea floor.
   maxdepth<-params[2]
  #Maximum vertical bearing (folding observations over vertical axis)
  #Note that all calculations are done using incident angles - if bearing data are used, then all the data need to be converted into incident angles.
  if (data=="bearing.inc" | data=="range") { 
  maxbearing<-params[3]}
  if (data=="bearing.elev") { 
    maxbearing<-(pi/2)-params[3]}

  #Line length - for towed scenarios only
  if (type=="towed") {
    linelength<-params[4] }
  
  #distribution parameters
  #for either uniform or beta - this allows distributions that do not range between 0 and maxdepth

  #Spread can be used with one of the two limits instead, particularly if the depth varies in the case of multiple sensors.

  if (distparams[1]=="NULL"){
    if (distparams[2]=="maxdepth"){ 
    ul<-maxdepth} else{
    ul<-as.numeric(distparams[2])}
    spread<-as.numeric(distparams[3])
    ll<-ul-spread}
   
  if (distparams[2]=="NULL"){
    ll<-as.numeric(distparams[1])
    spread<-as.numeric(distparams[3])
    ul<-ll+spread}
    
  if (distparams[3]=="NULL"){
    ll<-as.numeric(distparams[1])
    if (distparams[2]=="maxdepth"){ 
      ul<-maxdepth} else{
        ul<-as.numeric(distparams[2])}
    spread<-ul-ll}
 
  #check that the ll, ul and spread have been generated
  if(exists("ll")==FALSE){
  stop("Lower limit, upper limit and spread have not been created - check that input values make sense")}
  
  #Enter shape parameters if a beta distribution is being used

  if (dist=="beta"){
    
    shape1<-as.numeric(distparams[4])  
    shape2<-as.numeric(distparams[5])
    
  }

  #NB: The maximum vertical bearing is an incidence angle i.e., vertical is 0 rad (0 deg), horizontal is pi/2 rad (90 deg).
  #If maxbearing is less than pi/2, then will need to calculate max_xy, as this will be different from maxdetrad
    
      max_xy<-maxdetrad*sin(maxbearing) 
  
  ################################################################################
  #STEP 2: Enter the observations

    if (data=="range"){
      #make the vertical radial distances of these detections the data for the likelihood
      vertr_obs<-obs
      if (max(vertr_obs)>maxdetrad){stop("Error: observed ranges exceeds maximum detection radius")}
      results.obs<-length(vertr_obs)} else {
        #note that in either bearing case, the code is written using the incidence angle
        #if in an analysis, a user specified that they had bearing.elev data, then these would simply be converted to incidence angles.
       
        #In either case, it would be assumed that the tilt of the instrument was known (assumed to be zero here).
        
        #now estimate the vertical bearing i.e., the elevation angle (we ignore the horizontal bearing - the detection function is symmetrical around the vertical axis)      
        if (data=="bearing.inc"){
        
        verttheta_obs<-obs
        if (max(verttheta_obs)>maxbearing){stop("Error: observed bearings exceeds maximum bearing")}
        results.obs<-length(verttheta_obs) } else {
        
          if (data=="bearing.elev"){
            
            verttheta_obs<-(pi/2)-obs   
            if (max(verttheta_obs)>maxbearing){stop("Error: observed bearings exceeds maximum bearing")}
        results.obs<-length(verttheta_obs) }}}
    
    ############################################################################
    #STEP 3: Likelihood optimisation - using functions which need to be read in from "Depth.LH.functions.vs2.r"
    
    if (dist=="beta" & data=="range"){
      maxLH<-optimize(model.log.lik.fct.beta.r_obs,c(det.func.params[1],det.func.params[2]),shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,vertr_obs)}
    
    if (dist=="beta" && (data=="bearing.inc" | data=="bearing.elev")){
      maxLH<-optimize(model.log.lik.fct.beta.theta_obs,c(det.func.params[1],det.func.params[2]),shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,verttheta_obs)}
    
    if (dist=="uniform" & data=="range"){
      maxLH<-optimize(model.log.lik.fct.uniform.r_obs,c(det.func.params[1],det.func.params[2]),max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,vertr_obs)}
    
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="towed"){
      stop("A towed scenario with a uniform depth distribution and bearing data cannot be analysed")}
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="fixed"){
      maxLH<-optimize(model.log.lik.fct.uniform.theta_obs,c(det.func.params[1],det.func.params[2]),max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,verttheta_obs)}
    
    results.sig<-maxLH$minimum
    
  #########################################################################
    #
    #STEP 4 - Estimate the average probability of detection.  This is the denominator of the LH function using the maximum LH value.  
    
    sigmaest<-as.numeric(maxLH[1])
    
 #######################################################################
#Define the range of r and theta
#Ideally want small cutpoints for gridded integration.  Initially had 50 but these were too coarse so using 200.   However, this has made the routines a lot slower.
  
#In this new version of the LH functions, have tried to speed up the code

rrange<-seq(0,maxdetrad,length.out=201)
#get the midpoint between the range cuts
meanr.vec<-seq(rrange[2]/2,maxdetrad-(rrange[2]/2),length.out=200)
thetarange<-seq(0,maxbearing,length.out=201)
#get the midpoint between the theta cuts
meant.vec<-seq(thetarange[2]/2,maxbearing-(thetarange[2]/2),length.out=200)
#expand all combinations of r and theta midpoints so that can evaluate the functions without for loops
#transform so that have a lot of columns and fewer rows - R is quicker at column computation.
meanrthetagrid<-t(expand.grid(meanr.vec,meant.vec))
   
    ###############################################################################
    #calculate the denominator for either scenario with a beta distribution

if (dist=="beta"){ 
  denom1pdet<-array(0,dim(meanrthetagrid)[2])
      for (i in 1:dim(meanrthetagrid)[2]){  
        #evaluate the functions
        if (type=="fixed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1pdet[i]<-0} else
        {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
        if (type=="towed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1pdet[i]<-0} else
        {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_towed(meanrthetagrid[1,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
      }
       
#Multiply the results by the length of each range and theta cell.
denom2pdet<-denom1pdet*rrange[2]*thetarange[2]
    
#Sum the results of denom2pdet to get the final integral.
 denom3pdet<-sum(denom2pdet)
 results.pdet<-denom3pdet
}    
    
    ###############################################################################
    #calculate the denominator for either scenario with a uniform distribution

if (dist=="uniform"){   
  denom1pdet<-array(0,dim(meanrthetagrid)[2])
  for (i in 1:dim(meanrthetagrid)[2]){  
    #evaluate the functions
    if (type=="fixed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1pdet[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1pdet[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1pdet[i]<-0} else
        {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*uniformpdf_rtheta(spread)}}
    if (type=="towed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1pdet[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1pdet[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1pdet[i]<-0} else
        {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_towed(meanrthetagrid[1,i],max_xy)*uniformpdf_rtheta(spread)}}
      }
    
    
#Multiply the results by the length of each range and theta cell.
denom2pdet<-denom1pdet*rrange[2]*thetarange[2]
    
#Sum the results of denom2 to get the final integral.
denom3pdet<-sum(denom2pdet)
results.pdet<-denom3pdet
}
       
    #########################################################################
    #STEP 5 - work up to estimating density
    #
    #corrected number of animals in the cylinder or rectangle.
    n.corrected<-results.obs/results.pdet
    if (type=="fixed"){
      results.est.cyl.max_xy<-n.corrected}
    if (type=="towed"){
      results.est.rect.max_xy<-n.corrected}
   

  
  #this is abundance, want density in the horizontal plane i.e., per km2
  if (type=="fixed"){
    results.dens<-results.est.cyl.max_xy/(pi*((max_xy/1000)^2))}
  if (type=="towed"){
    results.dens<-results.est.rect.max_xy/(2*(max_xy/1000)*(linelength/1000))}

  ############################################################ 
  #STEP 6 - Return useful output
  
  if (type=="fixed") {
    
    results.list<-list("results.sig"=results.sig,
                       "results.pdet"=results.pdet,
                       "results.est.cyl.max_xy"=results.est.cyl.max_xy,
                       "results.dens"=results.dens,
                       "results.obs"=results.obs)} else
                         
  if (type=="towed") {
                           
        results.list<-list("results.sig"=results.sig,
                           "results.pdet"=results.pdet,
                           "results.est.rect.max_xy"=results.est.rect.max_xy,
                            "results.dens"=results.dens,
                            "results.obs"=results.obs)} 
 
  #Print useful output
  
  #Print information about the simulation set up
  print1<-print(paste("The survey is a", type, "survey"))
  cat(print1,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  if (data=="range"){
    datatype<-"ranges."}
  if (data=="bearing.inc"){
    datatype<-"incident angles i.e., vertical is zero."}
  if (data=="bearing.elev"){
    datatype<-"elevation angles i.e., horizontal is zero."}
  print2<-print(paste("The collected data are", datatype))
  print3<-print(paste("The depth distribution is assumed to be",dist,"."))
  print4<-print(paste("The maximum detection radius was assumed to be", params[1],"m."))
  print5<-print(paste("The maximum depth was assumed to be", params[2],"m."))
  print6<-print(paste("The maximum vertical detection bearing was assumed to be ", round(params[3],2),"rad."))
  cat(c(print2,print3,print4,print5,print6),file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)


  if (type=="towed"){
    print7<-print(paste("The survey line length was assumed to be", params[4],"."))
    cat(print7,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}
  
  if (dist=="beta"){
    print8<-print(paste("The shape parameters of the beta distribution were", shape1, "and", shape2,"."))
  cat(print8,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}
  print9<-print(paste("The range of values used to estimate sigma were", det.func.params[1], "and",det.func.params[2]))
  cat(print9,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  print9b<-print(paste("The lower limit,upper limit and spread of the depth distribution parameters were", ll, ul, "and", spread,"."))
  cat(print9b,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  #Now print some results
  print10<-print(paste("The estimate of sigma was", round(mean(results.sig),2)))
  print11<-print(paste("The estimate of the detection probability was", round(mean(results.pdet),2)))
  print12<-print(paste("The estimated density of animals was", round(mean(results.dens),2),"animals/km2"))
  print13<-print(paste("The number of detections was", round(mean(results.obs),2)))
  cat(c(print10,print11,print12,print13),file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  if (round(mean(results.sig),2)==det.func.params[2]){
  print14<-print("NOTE: The estimate of sigma seems to have hit the upper bound - suggest expanding the range of values in the optimisation.")
  cat(print14,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}

  
  #plot the detection function
  jpeg(paste('Analysis_plot_',label,'.jpg'))
  
plot(c(0:maxdetrad),grtheta(c(0:maxdetrad),results.sig),main=paste("Estimated sigma =",round(results.sig,2)),xlab="Slant range(m)",ylab="Probability of detection")
 
dev.off()

save(results.list, file = paste(label,".RData"))
return(c(results.sig,results.pdet))
}

