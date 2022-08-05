

multisensor.depth.simulation.vs2<-function(type, data, dist, params1, params2, params3, params4, popn, iterations, distparams,det.func.params,label){
 ############################################################################################
  #Step 1: Create results vectors

  
#Params will become a matrix of different detection radii,max depths,bearings and line lengths.

  params<-array(NA,c(length(params1),4))
  params[,1]<-t(params1)
  params[,2]<-t(params2)
  params[,3]<-t(params3)
  params[,4]<-t(params4)
  
#Results vectors - these can all be combined using cbind to create results matrices and saved as one list.  Haven't found a way to save all the same results as the simpler version
  
  if (type=="fixed") {
    
    results.cyl.mdr<-array(NA,c(dim(params)[1],iterations))
    results.cyl.maxxy<-array(NA,c(dim(params)[1],iterations))
    results.cone<-array(NA,c(dim(params)[1],iterations))
    results.sph<-array(NA,c(dim(params)[1],iterations))
    results.sig<-array(NA,iterations)
    results.pdet<-array(NA,c(dim(params)[1],iterations))
    results.est.cyl.max_xy<-array(NA,c(dim(params)[1],iterations))
    results.dens<-array(NA,c(dim(params)[1],iterations))
    results.obs<-array(NA,c(dim(params)[1],iterations))
    meanpercentbias<-numeric()
    } else 
      
    if (type=="towed") {
        
    results.rect.mdr<-array(NA,c(dim(params)[1],iterations))
    results.rect.maxxy<-array(NA,c(dim(params)[1],iterations))
    results.triangle<-array(NA,c(dim(params)[1],iterations))
    results.semicircle<-array(NA,c(dim(params)[1],iterations))
    results.sig<-array(NA,iterations)
    results.pdet<-array(NA,c(dim(params)[1],iterations))
    results.est.rect.max_xy<-array(NA,c(dim(params)[1],iterations))
    results.dens<-array(NA,c(dim(params)[1],iterations))
    results.obs<-array(NA,c(dim(params)[1],iterations))
    meanpercentbias<-numeric()
    } 
  
  ###################################################################################################
  #Step 2: Enter the parameters of the survey area and simulation 
  
  #Create multiple 3D survey areas
  
  #Maximum radial distance of the detection process - this is the maximum slant detection range for either a towed or fixed transect:
  maxdetrads<-params[,1]
  #Maximum depth (elevation in a fixed scenario)
  #Note that the calculations below work in exactly the same way for towed and fixed scenarios. 
  #The only difference is that depth in the towed scenario is actual depth i.e., 2000m = the seafloor
  #In contrast, in the fixed scenario, depth is actually elevation i.e., 2000m = the sea surface.
  #The metric is actually the distance between the instrument and a boundary, either the sea surface or sea floor.
  
  maxdepths<-params[,2]
  #Maximum vertical bearing (folding observations over vertical axis)

  if (data=="bearing.inc" | data=="range") { 
    maxbearings<-params[,3]}
  if (data=="bearing.elev") { 
    maxbearings<-(pi/2)-params[,3]}
  
  #Line length - for towed scenarios only
  if (type=="towed") {
    linelengths<-params[,4] }
  
  #density/number of animals in the population
  if (popn[1]=="abundance"){
  totpops<-c(rep(as.numeric(popn[2]),dim(params)[1]))
  #Calculate the true density and print the result
  if (type=="fixed"){
    truedensities.km2<-totpops/((2*maxdetrads/1000)*(2*maxdetrads/1000))} 
  if  (type=="towed"){
    truedensities.km2<-totpops/((2*maxdetrads/1000)*(linelengths/1000))}
  }
  if (popn[1]=="density"){
    truedensities.km2<-c(rep(as.numeric(popn[2]),dim(params)[1]))
    #Calculate the true abundance 
    if (type=="fixed"){
      totpops<-as.numeric(popn[2])*((2*maxdetrads/1000)*(2*maxdetrads/1000))} 
    if  (type=="towed"){
      totpops<-as.numeric(popn[2])*((2*maxdetrads/1000)*(linelengths/1000))} 
    }
  
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
        uls<-c(rep(as.numeric(distparams[2]),dim(params)[1]))}
    spreads<-uls-lls}

  #check that the ll, ul and spread have been generated
  if(exists("lls")==FALSE){
    stop("Lower limit, upper limit and spread have not been created - check that input values make sense")}
  
    #this could be generalised for other types of whale calling depth bands 
  
  if (dist=="beta"){
  
  shape1<-as.numeric(distparams[4])  
  shape2<-as.numeric(distparams[5])
  
  }
 
  #NB: The maximum vertical bearing is an incidence angle i.e., vertical is 0 rad (0 deg), horizontal is pi/2 rad (90 deg).
  #If maxbearing is less than pi/2, then will need to calculate max_xy, as this will be different from maxdetrad
  

    max_xys<-maxdetrads*sin(maxbearings) 


  #####################################################################################  
  #Step 3: Start the simulation
  
  for (n in 1:iterations){
 
 results.all.obs<-multisensor_optimize_1(type, data, dist, maxdetrads, maxdepths, maxbearings, max_xys,linelengths, popn, iterations, lls, uls, spreads,det.func.params,totpops,shape1,shape2)
 
 
     for (s in 1:dim(params)[1]){
     results.obs[s,n]<-length(results.all.obs[[s]]) 
       if (type=="fixed") {
    
      results.cyl.mdr[s,n]<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[1]
      results.cyl.maxxy[s,n]<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[2]
      results.cone[s,n]<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[3]
      results.sph[s,n]<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[4]
      } else 
   
      if (type=="towed") {
     
          results.rect.mdr<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[1]
          results.rect.maxxy<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[2]
          results.triangle<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[3]
          results.semicircle<-as.numeric(results.all.obs[s+length(maxdetrads)][[1]])[4]
      } 
    }
  
    ############################################################################
    #STEP 4: Likelihood optimisation 
   
maxLHs<-optimize(multisensor_optimize_2,c(det.func.params[2],det.func.params[3]),type, data, dist, maxdetrads, maxdepths, maxbearings, max_xys,linelengths, popn, iterations, lls, uls, spreads,det.func.params,results.all.obs,shape1,shape2)

    results.sig[n]<-maxLHs$minimum
        
    #########################################################################
    #
    #STEP 5 - Estimate the average probability of detection.  This is the denominator of the LH function using the maximum LH value.  
    
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
        results.pdet[s,n]<-denom3pdet
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
      results.pdet[s,n]<-denom3pdet
    }
    
    }
   print(n)
  }    
    
    #########################################################################
    #STEP 6 - work up to estimating density
    #
    #corrected number of animals in the cylinder or rectangle.
    n.corrected<-results.obs/results.pdet
    if (type=="fixed"){
      results.est.cyl.max_xy<-n.corrected}
    if (type=="towed"){
      results.est.rect.max_xy<-n.corrected}
    


  #this is abundance, want density in the horizontal plane i.e., per km2
  if (type=="fixed"){
    for (s in 1:dim(params)[1]){
    results.dens[s,]<-results.est.cyl.max_xy[s,]/(pi*((max_xys[s]/1000)^2))}}
  if (type=="towed"){
    for (s in 1:dim(params)[1]){
    results.dens[s,]<-results.est.rect.max_xy[s,]/(2*(max_xys[s]/1000)*(linelengths[s]/1000))}}
  
  #percentage bias of true density
for (s in 1:dim(params)[1]){
  meanpercentbias[s]<-((mean(results.dens[s,])-truedensities.km2[s])/truedensities.km2[s])*100}
  #percentage bias of true sigma
  meanpercentbias_sigma<-((mean(results.sig)-det.func.params[1])/det.func.params[1])*100


  if (type=="fixed") {
    
    results.list<-list(
    "results.cyl.mdr"=results.cyl.mdr,
    "results.cyl.maxxy"=results.cyl.maxxy,
    "results.cone"=results.cone,
    "results.sph"=results.sph,
    "results.sig"=results.sig,
    "results.pdet"=results.pdet,
    "results.est.cyl.max_xy"=results.est.cyl.max_xy,
    "results.dens"=results.dens,
    "results.obs"=results.obs, 
    "meanpercentbias"=meanpercentbias, 
    "meanpercentbias_sigma"=meanpercentbias_sigma) } else
          
      if (type=="towed") {
        
        results.list<-list(
        "results.rect.mdr"=results.rect.mdr,
        "results.rect.maxxy"=results.rect.maxxy,
        "results.triangle"=results.triangle,
        "results.semicircle"=results.semicircle,
        "results.sig"=results.sig,
        "results.pdet"=results.pdet,
        "results.est.rect.max_xy"=results.est.rect.max_xy,
        "results.dens"=results.dens,
        "results.obs"=results.obs,
        "meanpercentbias"=meanpercentbias,
        "meanpercentbias_sigma"=meanpercentbias_sigma) } 
  
  #return(results.list)
  ##############################################################
  #STEP 7: Print useful output
  
  #Print information about the simulation set up
  print1<-print(paste("The survey is a", type, "survey"))
  cat(print1,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
print1aa<-print(paste("The number of separate monitoring lines or points with different parameters are",dim(params)[1]))
cat(print1aa,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  if (data=="range"){
    datatype<-"ranges."}
  if (data=="bearing.inc"){
    datatype<-"incident angles i.e., vertical is zero."}
  if (data=="bearing.elev"){
    datatype<-"elevation angles i.e., horizontal is zero."}
  print1a<-print(paste("The animal population metric set in this simulation was",popn[1]))
  print2a<-print("True densities are:")
  print2b<-print(round(truedensities.km2,2))
  print2c<-print("animals per km2") 
  print3<-print(paste("The collected data are", datatype))
  print4<-print(paste("The depth distribution is assumed to be",dist,"."))
  print5a<-print("The maximum detection radii were assumed to be:")
  print5b<-print(t(params[,1]))
  print5c<-print("m.")
  print6a<-print("The maximum depths were assumed to be")
  print6b<-print(params[,2])
  print6c<-print("m.")
  print7a<-print("The maximum vertical detection bearings were assumed to be ")
  print7b<-print(round(params[,3],2))
  print7c<-print("rad.")
  cat(c(print1a,print2a,print2b,print2c,print3,print4,print5a,print5b,print5c,print6a,print6b,print6c,print7a,print7b,print7c),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  if (type=="towed"){
    print8a<-print("The survey line lengths were assumed to be")
    print8b<-print(params[,4])
    print8c<-print("m.")
    cat(c(print8a,print8b,print8c),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)}
    print9a<-print("The number of simulated animals were")
    print9b<-print(totpops)
    print10<-print(paste("The simulation was run", iterations, "times."))
    cat(c(print9a,print9b,print10),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
print11a<-print(paste("The lower limit,upper limit and spread of the depth distribution parameters were set to be", distparams[1],",", distparams[2], "and", distparams[3],"."))
  cat(print11a,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  if (dist=="beta"){
    print11<-print(paste("The shape parameters of the beta distribution were", shape1, "and", shape2,"."))
    cat(print11,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)}
    print12<-print(paste("The sigma value used in the half normal detection function was",det.func.params[1]," and the range of values used in the optimisation were",det.func.params[2], "and",det.func.params[3]))
  
  #Now print some results
  print13<-print(paste("The mean estimate of sigma was", round(mean(results.sig),2)))
  print14a<-print("The mean estimates of the detection probabilities were")
  print14b<-print(t(round(apply(results.pdet,1,mean),2)))
  print15a<-print("The mean densities of animals were")
  print15b<-print(t(round(apply(results.dens,1,mean),0)))
  print15c<-print("animals/km2")
  print16a<-print("The range of estimated densities was between")
  print16b<-print(paste(round(min(results.dens),0), "and" ,round(max(results.dens),0),"animals/km2."))
  print17a<-print("The mean percentage bias in the density estimate was")
  print17b<-print(round(meanpercentbias,0))
  print17c<-print("%.")
  print18a<-print("The mean percentage bias in the sigma estimate was")
  print18b<-print(round(meanpercentbias_sigma,2))
  print18c<-print("%.")
  print19a<-print("The mean number of detections per line/point were")
  print19b<-print(t(round(apply(results.obs,1,mean),0)))
  cat(c(print12,print13,print14a,print14b,print15a,print15b,print15c,print16a,print16b,print17a,print17b,print17c,print18a,print18b,print18c,print19a,print19b),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  
  #plot histogram of sigma estimates
  
  jpeg(paste('Sim_plot_',label,'_hist_sig.jpg'))
  hist(results.sig)
  dev.off()

  save(results.list, file = paste(label,".RData"))

  #plot the most useful plot - the set up of the simulation - one plot for each instrument

  for (s in 1:dim(params)[1]){
  
  if (type=="towed"){
  x_totpop<-results.all.obs[[((s-1)*5)+1+(2*length(maxdetrads))]]
  z_totpop<-results.all.obs[[((s-1)*5)+2+(2*length(maxdetrads))]]
  vertr_totpop_edit<-results.all.obs[[((s-1)*5)+3+(2*length(maxdetrads))]]
  inctheta_totpop_edit<-results.all.obs[[((s-1)*5)+4+(2*length(maxdetrads))]]
  det_totpop<-results.all.obs[[((s-1)*5)+5+(2*length(maxdetrads))]]}
  
  if (type=="fixed"){
  x_totpop_edit<-results.all.obs[[((s-1)*5)+1+(2*length(maxdetrads))]]
  z_totpop_edit<-results.all.obs[[((s-1)*5)+2+(2*length(maxdetrads))]]
  vertr_totpop_edit<-results.all.obs[[((s-1)*5)+3+(2*length(maxdetrads))]]
  inctheta_totpop_edit<-results.all.obs[[((s-1)*5)+4+(2*length(maxdetrads))]]
  det_totpop<-results.all.obs[[((s-1)*5)+5+(2*length(maxdetrads))]]}
 
  if (type=="fixed"){
  
    jpeg(paste('Sim_plot_',label,s,'.jpg'))
    
    plot(x_totpop_edit,z_totpop_edit,xlim=c(-maxdetrads[s],maxdetrads[s]),ylim=c(0,maxdepths[s]),ylab=("Elevation (m)"),xlab=("Horizontal range (m)"))
    points(x_totpop_edit[!vertr_totpop_edit>maxdetrads[s]],z_totpop_edit[!vertr_totpop_edit>maxdetrads[s]],col="red")
    points(x_totpop_edit[!inctheta_totpop_edit>maxbearings[s]],z_totpop_edit[!inctheta_totpop_edit>maxbearings[s]],col="green")
    points(x_totpop_edit[!inctheta_totpop_edit>maxbearings[s] & !vertr_totpop_edit>maxdetrads[s]],z_totpop_edit[!inctheta_totpop_edit>maxbearings[s] & !vertr_totpop_edit>maxdetrads[s]],col="purple") 
    points(x_totpop_edit[det_totpop==1],z_totpop_edit[det_totpop==1],pch=16) 
    abline(v=c(-max_xys[s],max_xys[s]))
    
    dev.off()
  }
  
  if (type=="towed"){
    
    jpeg(paste('Sim_plot_',label,s,'.jpg'))
    
    plot(x_totpop,z_totpop,xlim=c(-maxdetrads[s],maxdetrads[s]),ylim=c(maxdepths[s],0),ylab=("Depth (m)"),xlab=("Horizontal range (m)"))
    points(x_totpop[!vertr_totpop_edit>maxdetrads[s]],z_totpop[!vertr_totpop_edit>maxdetrads[s]],col="red")
    points(x_totpop[!inctheta_totpop_edit>maxbearings[s]],z_totpop[!inctheta_totpop_edit>maxbearings[s]],col="green")
    points(x_totpop[!inctheta_totpop_edit>maxbearings[s] & !vertr_totpop_edit>maxdetrads[s]],z_totpop[!inctheta_totpop_edit>maxbearings[s] & !vertr_totpop_edit>maxdetrads[s]],col="purple") 
    points(x_totpop[det_totpop==1],z_totpop[det_totpop==1],pch=16)
    abline(v=c(-max_xys[s],max_xys[s]))
    
    dev.off()
  }
   }
}

