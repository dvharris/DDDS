multisensor.depth.analysis.vs2<-function(type, data, dist, params1, params2, params3, params4, distparams,det.func.params,obs,label){

#NB: det.func.params is only 2 values - range to look over for optimization.
#obs needs to be a list of all observations, split by sensor. 
  ################################################################################
  #STEP 1: Set up the analysis
  
  #Params will become a matrix of different detection radii,max depths,bearings and line lengths.

  params<-array(NA,c(length(params1),4))
  params[,1]<-t(params1)
  params[,2]<-t(params2)
  params[,3]<-t(params3)
  params[,4]<-t(params4)
  
  #Results vectors - one entry for each instrument
  
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
   maxdetrads<-params[,1]
  #Maximum depth#Note that the calculations below work in exactly the same way for towed and fixed scenarios. 
  #The only difference is that depth in the towed scenario is actual depth i.e., 2000m = the seafloor
  #In contrast, in the fixed scenario, depth is actually elevation i.e., 2000m = the sea surface.
  #The metric is actually the distance between the instrument and a boundary, either the sea surface or sea floor.
  maxdepths<-params[,2]
  #Maximum vertical bearing (folding observations over vertical axis)
  #Note that all calculations are done using incident angles - if bearing data are used, then all the data need to be converted into incident angles.  
  if (data=="bearing.inc" | data=="range") { 
    maxbearings<-params[,3]}
  if (data=="bearing.elev") { 
    maxbearings<-(pi/2)-params[,3]}
  
  #Line length - for towed scenarios only
  if (type=="towed") {
    linelengths<-params[,4] }
  
  #distribution parameters
  #for either uniform or beta - this allows distributions that do not range between 0 and maxdepth

  #Spread can be used with one of the two limits instead, particularly if the depth varies in the case of multiple sensors.
  if (distparams[1]=="NULL"){
    if (distparams[2]=="maxdepth"){ 
      uls<-maxdepths} else{
        uls<-c(rep(as.numeric(distparams[2]),dim(params)[1]))}
    spreads<-c(rep(as.numeric(distparams[3]),dim(params)[1]))
    lls<-uls-spreads}
  
  if (distparams[2]=="NULL"){
    lls<-c(rep(as.numeric(distparams[1]),dim(params)[1]))
    spreads<-c(rep(as.numeric(distparams[3]),dim(params)[1]))
    uls<-lls+spreads}
  
  if (distparams[3]=="NULL"){
    lls<-c(rep(as.numeric(distparams[1]),dim(params)[1]))
    if (distparams[2]=="maxdepth"){ 
      uls<-maxdepths} else{
        uls<-lls<-c(rep(as.numeric(distparams[2]),dim(params)[1]))}
    spreads<-uls-lls}
 
  #check that the ll, ul and spread have been generated
 if(exists("lls")==FALSE){
    stop("Lower limit, upper limit and spread have not been created - check that input values make sense")}
  
  #Enter shape parameters if a beta distribution is being used  
  if (dist=="beta"){
    
    shape1<-as.numeric(distparams[4])  
    shape2<-as.numeric(distparams[5])
    
  }
  
  #NB: The maximum vertical bearing is an incidence angle i.e., vertical is 0 rad (0 deg), horizontal is pi/2 rad (90 deg).
  #If maxbearing is less than pi/2, then will need to calculate max_xy, as this will be different from maxdetrad
  
  max_xys<-maxdetrads*sin(maxbearings) 
  
    
  ################################################################################
  #STEP 2: Enter the observations
    if (data=="range"){
      #make the vertical radial distances of these detections the data for the likelihood
      results.all.obs<-obs
      for (s in 1:length(results.all.obs)){
      if (max(results.all.obs[[s]])>maxdetrads[s]){stop("Error: observed ranges exceeds maximum detection radius")}}
      for (s in 1:length(results.all.obs)){
      results.obs[s]<-length(results.all.obs[[s]])}} else {
        #note that in either bearing case, the code is written using the incidence angle
        #if in an analysis, a user specified that they had bearing.elev data, then these would simply be converted to incidence angles.
       
        #In either case, it would be assumed that the tilt of the instrument was known (assumed to be zero here).
        
        #now estimate the vertical bearing i.e., the elevation angle (we ignore the horizontal bearing - the detection function is symmetrical around the vertical axis)      
        if (data=="bearing.inc"){
        
        results.all.obs<-obs
        for (s in 1:length(results.all.obs)){
        if (max(results.all.obs[[s]])>maxbearings[s]){stop("Error: observed bearings exceeds maximum bearing")}}
        for (s in 1:length(results.all.obs)){
      results.obs[s]<-length(results.all.obs[[s]])}} else {
        
          if (data=="bearing.elev"){
          
          results.all.obs<-list()
          for (s in 1:length(obs)){  
           results.all.obs[[s]]<-(pi/2)-obs[[s]]   
           for (s in 1:length(results.all.obs)){
        if (max(results.all.obs[[s]])>maxbearings[s]){stop("Error: observed bearings exceeds maximum bearing")}}
        results.obs[s]<-length(results.all.obs[[s]]) }}}}
 
 ###################################################################################   
  #STEP 3: Likelihood optimisation
    
    
maxLHs<-optimize(multisensor_optimize_2,c(det.func.params[1],det.func.params[2]),type, data, dist, maxdetrads, maxdepths, maxbearings, max_xys,linelengths, popn, iterations, lls, uls, spreads,det.func.params,results.all.obs,shape1,shape2)

  
results.sig<-maxLHs$minimum


#########################################################################
    #
    #STEP 4 - Estimate the average probability of detection.  This is the denominator of the LH function using the maximum LH value.  
    
    sigmaest<-as.numeric(maxLHs[1])
    
    ###############################################################################
    #calculate the denominator for either scenario with a beta distribution
    #Note that there is a separate denominator and pdet for each instruments scenario    
    

    for (s in 1:dim(params)[1]){
    
      #Set up the integration for the denominator
      #############################################################################
      #Define the range of r and theta (the two parameters)
      
      rrange<-seq(0,maxdetrads[s],length.out=201)
      #get the midpoint between the range cuts
      meanr.vec<-seq(rrange[2]/2,maxdetrads[s]-(rrange[2]/2),length.out=200)
      thetarange<-seq(0,maxbearings[s],length.out=201)
      #get the midpoint between the theta cuts
      meant.vec<-seq(thetarange[2]/2,maxbearings[s]-(thetarange[2]/2),length.out=200)
      #expand all combinations of r and theta midpoints so that can evaluate the functions without for loops
      #transform so that have a lot of columns and fewer rows - R is quicker at column computation.
      meanrthetagrid<-t(expand.grid(meanr.vec,meant.vec))
      
      if (dist=="beta"){ 
        denom1pdet<-array(0,dim(meanrthetagrid)[2])
        for (i in 1:dim(meanrthetagrid)[2]){  
          #evaluate the functions
          if (type=="fixed"){
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepths[s]) {denom1pdet[i]<-0} else 
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > uls[s]) {denom1pdet[i]<-0} else 
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < lls[s]) {denom1pdet[i]<-0} else
                {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xys[s])*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,lls[s],uls[s])}}
          if (type=="towed"){
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepths[s]) {denom1pdet[i]<-0} else 
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > uls[s]) {denom1pdet[i]<-0} else 
            if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < lls[s]) {denom1pdet[i]<-0} else
                {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_towed(meanrthetagrid[1,i],max_xys[s])*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,lls[s],uls[s])}}
        }
        
        #Multiply the results by the length of each range and theta cell.
        denom2pdet<-denom1pdet*rrange[2]*thetarange[2]
        
        #Sum the results of denom2 to get the final integral.
        denom3pdet<-sum(denom2pdet)
        results.pdet[s]<-denom3pdet
      }
       
    ###############################################################################
    #calculate the denominator for either scenario with a uniform distribution

    if (dist=="uniform"){   
      denom1pdet<-array(0,dim(meanrthetagrid)[2])
      for (i in 1:dim(meanrthetagrid)[2]){  
        #evaluate the functions
        if (type=="fixed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepths[s]) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > uls[s]) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < lls[s]) {denom1pdet[i]<-0} else
              {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xys[s])*uniformpdf_rtheta(spreads[s])}}
        if (type=="towed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepths[s]) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > uls[s]) {denom1pdet[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < lls[s]) {denom1pdet[i]<-0} else
              {denom1pdet[i]<-grtheta(meanrthetagrid[1,i],sigmaest)*rxdist_towed(meanrthetagrid[1,i],max_xys[s])*uniformpdf_rtheta(spreads[s])}}
      }
      
      
      #Multiply the results by the length of each range and theta cell.
      denom2pdet<-denom1pdet*rrange[2]*thetarange[2]
      
      #Sum the results of denom2 to get the final integral.
      denom3pdet<-sum(denom2pdet)
      results.pdet[s]<-denom3pdet
    }
    print(s)
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
    results.dens<-results.est.cyl.max_xy/(pi*((max_xys/1000)^2))}
  if (type=="towed"){
    results.dens<-results.est.rect.max_xy/(2*(max_xys/1000)*(linelengths/1000))}
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
  print1aa<-print(paste("The number of separate monitoring lines or points with different parameters are",dim(params)[1]))
cat(print1aa,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  if (data=="range"){
    datatype<-"ranges."}
  if (data=="bearing.inc"){
    datatype<-"incident angles i.e., vertical is zero."}
  if (data=="bearing.elev"){
    datatype<-"elevation angles i.e., horizontal is zero."}
  print2<-print(paste("The collected data are", datatype))
  print3<-print(paste("The depth distribution is assumed to be",dist,"."))
  print4a<-print("The maximum detection radii were assumed to be:")
  print4b<-print(params[,1])
  print4c<-print("m.")
  print5a<-print("The maximum depths were assumed to be")
  print5b<-print(params[,2])
  print5c<-print("m.")
  print6a<-print("The maximum vertical detection bearings were assumed to be ")
  print6b<-print(round(params[,3],2))
  print6c<-print("rad.")
  cat(c(print2,print3,print4a,print4b,print4c,print5a,print5b,print5c,print6a,print6b,print6c),file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  if (type=="towed"){
    print7a<-print("The survey line lengths were assumed to be")
    print7b<-print(params[,4])
    print7c<-print("m.")
    cat(c(print7a,print7b,print7c),file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}
  
  if (dist=="beta"){
    print8<-print(paste("The shape parameters of the beta distribution were", shape1, "and", shape2,"."))
  cat(print8,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}
  print9<-print(paste("The range of values used to estimate sigma were", det.func.params[1], "and",det.func.params[2]))
  cat(print9,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  print9b<-print(paste("The lower limit,upper limit and spread of the depth distribution parameters were set to be", distparams[1],",", distparams[2], "and", distparams[3],"."))
  cat(print9b,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  #Now print some results
  print10<-print(paste("The estimate of sigma was", round(mean(results.sig),2)))
print11a<-print("The estimated detection probabilities were")
  print11b<-print(round(results.pdet,2))
  print12a<-print("The estimated densities of animals were")
  print12b<-print(round(results.dens,0))
  print12c<-print("animals/km2")
  print13a<-print("The number of detections per line/point were")
  print13b<-print(results.obs)
  cat(c(print10,print11a,print11b,print12a,print12b,print12c,print13a,print13b),file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)
  if (round(results.sig,2)==det.func.params[2]){
  print14<-print("NOTE: The estimate of sigma seems to have hit the upper bound - suggest expanding the range of values in the optimisation.")
  cat(print14,file=paste("Analysis_output_",label,".txt"),sep="\n",append=TRUE)}
  
  #plot the detection function
  
  jpeg(paste('Analysis_plot_',label,'.jpg'))
  
for (s in 1:length(results.all.obs)){ 
plot(c(0:max(maxdetrads)),grtheta(c(0:max(maxdetrads)),results.sig),main=paste("Estimated sigma =",round(results.sig,2)),xlab="Slant range(m)",ylab="Probability of detection")
}
dev.off()

save(results.list, file = paste(label,".RData"))
return(c(results.sig,results.pdet))

}
