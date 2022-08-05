depth.simulation.vs2<-function(type, data, dist, params, popn, iterations, distparams,det.func.params,label){

  ############################################################################################
  #Step 1: Create results vectors

  if (type=="fixed") {
    
    results.cyl.mdr<-numeric()
    results.cyl.maxxy<-numeric()
    results.cone<-numeric()
    results.sph<-numeric()
    results.sig<-numeric()
    results.pdet<-numeric()
    results.est.cyl.max_xy<-numeric()
    results.dens<-numeric()
    results.obs<-numeric()} else 
      
  if (type=="towed") {
        
    results.rect.mdr<-numeric()
    results.rect.maxxy<-numeric()
    results.triangle<-numeric()
    results.semicircle<-numeric()
    results.sig<-numeric()
    results.pdet<-numeric()
    results.est.rect.max_xy<-numeric()
    results.dens<-numeric()
    results.obs<-numeric()} 

  ###################################################################################################
  #Step 2: Enter the parameters of the survey area and simulation 
  
  #Create a 3D survey area

  #Maximum radial distance of the detection process - this is the maximum slant detection range for either a towed or fixed transect:
  maxdetrad<-params[1]
  #Maximum depth (elevation in a fixed scenario)
  #Note that the calculations below work in exactly the same way for towed and fixed scenarios. 
  #The only difference is that depth in the towed scenario is actual depth i.e., 2000m = the seafloor
  #In contrast, in the fixed scenario, depth is actually elevation i.e., 2000m = the sea surface.
  #The metric is actually the distance between the instrument and a boundary, either the sea surface or sea floor.
  maxdepth<-params[2]
  #Maximum vertical bearing (folding observations over vertical axis)

  if (data=="bearing.inc" | data=="range") { 
    maxbearing<-params[3]}
  if (data=="bearing.elev") { 
    maxbearing<-(pi/2)-params[3]}
  
  #Line length - for towed scenarios only
  if (type=="towed") {
    linelength<-params[4] }
  
  #Density/number of animals in the population
  if (popn[1]=="abundance"){
  totpop<-as.numeric(popn[2])
  #Calculate the true density
  if (type=="fixed"){
    truedensity.km2<-totpop/((2*maxdetrad/1000)*(2*maxdetrad/1000))} 
  if  (type=="towed"){
    truedensity.km2<-totpop/((2*maxdetrad/1000)*(linelength/1000))}
  }
  if (popn[1]=="density"){
    truedensity.km2<-as.numeric(popn[2])
    #Calculate the true abundance 
    if (type=="fixed"){
      totpop<-as.numeric(popn[2])*((2*maxdetrad/1000)*(2*maxdetrad/1000))} 
    if  (type=="towed"){
      totpop<-as.numeric(popn[2])*((2*maxdetrad/1000)*(linelength/1000))} 
    }
  
  #Distribution parameters
  #for either uniform or beta - this allows for distributions that do not range between 0 and maxdepth
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
  #If maxbearing is less than pi/2, then will need to calculate max_xy (the maximum horizontal radius), as this will be different from maxdetrad
  
  if (maxbearing==pi/2){
    max_xy<-maxdetrad} else { 
      max_xy<-maxdetrad*sin(maxbearing)} 

  #####################################################################################  
  #Step 3: Start the simulation
  
  for (n in 1:iterations){
    
    #Need to assign each animal a location in the (x,z,y) dimension i.e., a radial distance with an azimuth and a depth
    
    x_totpop<-mat.or.vec(totpop, 1)
    z_totpop<-mat.or.vec(totpop, 1)
    y_totpop<-mat.or.vec(totpop, 1)
    
    for (i in 1:totpop){
      x_totpop[i]<-runif(1, min=-maxdetrad, max=maxdetrad)
      if  (type=="towed"){
        y_totpop[i]<-runif(1, min=0, max=linelength)} else 
          if (type=="fixed"){ 
            y_totpop[i]<-runif(1, min=-maxdetrad, max=maxdetrad)}
      
      if (dist == "beta"){ 
        z_totpop[i]<-(rbeta(1, shape1, shape2)*spread)+ll} else       
          if (dist == "uniform"){
            z_totpop[i]<-runif(1, min=ll, max=ul) } else {
              stop("Error - depth distribution not recognised")}
    }
    
    #USEFUL PLOTS - some plots are provided at the end of the function but the code for the plots below have been kept
    
    #PLOT 1 - look at depth vs horizontal range
    #plot(x_totpop,z_totpop, ylim = c(0,maxdepth),xaxt='n', main = "Depth against horizontal range")
    #axis(3)
    
    #PLOT 2 - look at a birdseye view of the point transect 
    #plot(x_totpop,y_totpop, main = "Birdseye view of simulated study area")
    
    #PLOT 3 - 3D plot
    # 3D Scatterplot package
    # libraries loaded once earlier
    # s3d<-scatterplot3d(x_totpop,y_totpop,z_totpop,cex.symbols=0.5,main="3D Scatterplot",zlim=c(0,maxdepth))
    # #add the location of the point - fixed only
    # s3d$points3d(0,0,0,col="red",pch=16)
    
    #Calculate the vertical radial distances and bearings for all animals - FIXED TRANSECTS ONLY
    
    if (type=="fixed"){
      #first need the horizontal radial distance for each animal
      #this is the relationship between x and y dimensions (use pythag)
      horizr_totpop<-sqrt((x_totpop^2)+(y_totpop^2))
      #NB:  don't want to use horizr values that are greater than the maxdetrad, as these are not in the population that I want to consider - don't want to pollute the detection process so need to remove them.
      horizr_totpop_edit<-horizr_totpop[!horizr_totpop>maxdetrad]
      #edit other vectors too
      x_totpop_edit<-x_totpop[!horizr_totpop>maxdetrad]
      y_totpop_edit<-y_totpop[!horizr_totpop>maxdetrad]
      z_totpop_edit<-z_totpop[!horizr_totpop>maxdetrad]
      
      #number animals in cylinder
      true.n.cylinder<-length(x_totpop_edit)
      #the line below is for when the true number changes for each iteration
      results.cyl.mdr[n]<-true.n.cylinder
      
      #now use the edited horizontal radial distance and the edited depth to estimate the vertical radial distance
      vertr_totpop_edit<-sqrt((horizr_totpop_edit^2)+(z_totpop_edit^2))
      #Note that these values include animals outside the detection area - any animal with a vertical radial distance greater than maxdetrad will have pdet of zero. 
      
      #check maximum vertical radial distance
      #Useful checks - left in for now
      #max_vertr<-sqrt((maxdepth^2)+(maxdetrad^2))
      #range(vertr_totpop_edit)
      
      #number animals in half sphere (which may be capped depending on depth)
      true.n.halfsphere<-length(x_totpop_edit[!vertr_totpop_edit>maxdetrad])
      #the line below is for when the true number changes for each iteration
      results.sph[n]<-true.n.halfsphere
      
      #now estimate the vertical bearing i.e., the elevation angle (we ignore the horizontal bearing - the detection function is symmetrical around the vertical axis)
      verttheta_totpop_edit<-atan(z_totpop_edit/horizr_totpop_edit)
      
      #now make these incident angles
      inctheta_totpop_edit<-(pi/2)-verttheta_totpop_edit
      #Note that these values include animals outside the detection area - any animal with a incidence angle greater than (pi/4) will have pdet of zero. 
      
      #number animals in cone - NB: this will be the same as the results.cyl.mdr/true cylinder if there is no bearing restriction
      true.n.cone<-length(x_totpop_edit[inctheta_totpop_edit<maxbearing])
      #the line below is for when the true number changes for each iteration
      results.cone[n]<-true.n.cone
      
      #number in cylinder with radius of max_xy:
      x_totpop_maxxy<-x_totpop[!horizr_totpop>max_xy]
      y_totpop_maxxy<-y_totpop[!horizr_totpop>max_xy]
      z_totpop_maxxy<-z_totpop[!horizr_totpop>max_xy]
      true.n.max_xy.cylinder<-length(x_totpop_maxxy)
      results.cyl.maxxy[n]<-true.n.max_xy.cylinder
      #Note that these values include animals outside the detection area - any animal with a vertical bearing greater than maxbearing will have pdet of zero. I.e., these will be detections outside the cone but in the cylinder.
      
    }
    
    #PLOT 4 - check largest horizontal ranges have been removed
    # plot(x_totpop,y_totpop,pch=16,col="red", main="Population to be detected")
    #points(x_totpop_edit,y_totpop_edit,pch=16)
    
    #PLOT 5 - 3D image of cylinder
    #s3dcylinder<-scatterplot3d(x_totpop_edit,y_totpop_edit,z_totpop_edit,cex.symbols=0.5,main="3D Scatterplot of study population",xlim=c(range(x_totpop)),ylim=c(range(y_totpop)),zlim=c(0,maxdepth))
    #s3dcylinder$points3d(x_totpop_edit[!vertr_totpop_edit>maxdetrad],y_totpop_edit[!vertr_totpop_edit>maxdetrad],z_totpop_edit[!vertr_totpop_edit>maxdetrad],col="red")
    
    #PLOT 6 - 2D image of hemisphere - may be capped depending on depth 
    #plot(x_totpop_edit,z_totpop_edit,xlim=c(-maxdetrad,maxdetrad),ylim=c(0,maxdepth))
    #points(x_totpop_edit[!vertr_totpop_edit>maxdetrad],z_totpop_edit[!vertr_totpop_edit>maxdetrad],col="red")
    
    #PLOT 7 - Incidence angle
    #s3dcylinder<-scatterplot3d(x_totpop_edit,y_totpop_edit,z_totpop_edit,cex.symbols=0.5,main="3D Scatterplot",xlim=c(range(x_totpop)),ylim=c(range(y_totpop)),zlim=c(range(z_totpop)))
    #s3dcylinder$points3d(x_totpop_edit[!inctheta_totpop_edit>maxbearing],y_totpop_edit[!inctheta_totpop_edit>maxbearing],z_totpop_edit[!inctheta_totpop_edit>maxbearing],col="red")
    
    #PLOT 8 - Incidence angle
    #2D plot
    #plot(x_totpop_edit[!inctheta_totpop_edit>maxbearing],z_totpop_edit[!inctheta_totpop_edit>maxbearing],ylim=c(0,maxdepth))
    
    #PLOT 9 - all relevant observations
    #detailed plot
    #plot(x_totpop_edit,z_totpop_edit,xlim=c(-maxdetrad,maxdetrad),ylim=c(0,maxdepth))
    #points(x_totpop_edit[!vertr_totpop_edit>maxdetrad],z_totpop_edit[!vertr_totpop_edit>maxdetrad],col="red")
    #points(x_totpop_edit[!inctheta_totpop_edit>maxbearing],z_totpop_edit[!inctheta_totpop_edit>maxbearing],col="green")
    #points(x_totpop_edit[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],z_totpop_edit[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],col="purple")  
    #abline(v=c(-max_xy,max_xy))
    
    #Calculate the vertical radial distances and bearings for all animals - TOWED TRANSECTS ONLY
    
    if (type=="towed"){
      
      #number animals in rectangle
      true.n.rectangle<-length(x_totpop)
      #the line below is for when the true number changes for each iteration
      results.rect.mdr[n]<-true.n.rectangle
      
      #now use the horizontal radial distance and the depths to estimate the vertical radial distance.
      #NB: the horizontal distances did not have to be edited first, but still called "edited" to be consistent with the fixed case.
      vertr_totpop_edit<-sqrt((x_totpop^2)+(z_totpop^2))
      #Note that these values include animals outside the detection area - any animal with a vertical radial distance greater than maxdetrad will have pdet of zero. 
      
      #check maximum vertical radial distance
      #Useful checks - left in for now
      #max_vertr<-sqrt((maxdepth^2)+(maxdetrad^2))
      #range(vertr_totpop_edit)
      
      #number animals in semicircle (which may be capped depending on depth)
      true.n.semicircle<-length(x_totpop[!vertr_totpop_edit>maxdetrad])
      #the line below is for when the true number changes for each iteration
      results.semicircle[n]<-true.n.semicircle
      
      #now estimate the vertical bearing i.e., the elevation angle 
      verttheta_totpop_edit<-atan(z_totpop/abs(x_totpop))
      
      #now make these incident angles
      inctheta_totpop_edit<-(pi/2)-verttheta_totpop_edit
      #Note that these values include animals outside the detection area - any animal with a incidence angle greater than (pi/4) will have pdet of zero. 
      
      #number animals in triangle - NB: this will be the same as the results.crect.mdr/true rectangle if there is no bearing restriction
      true.n.triangle<-length(x_totpop[inctheta_totpop_edit<maxbearing])
      #the line below is for when the true number changes for each iteration
      results.triangle[n]<-true.n.triangle
      
      #number in rectangle with radius of max_xy:
      x_totpop_maxxy<-x_totpop[!x_totpop>max_xy]
      y_totpop_maxxy<-y_totpop[!x_totpop>max_xy]
      z_totpop_maxxy<-z_totpop[!x_totpop>max_xy]
      true.n.max_xy.rectangle<-length(x_totpop_maxxy)
      results.rect.maxxy[n]<-true.n.max_xy.rectangle
      #Note that these values include animals outside the detection area - any animal with a vertical bearing greater than maxbearing will have pdet of zero. I.e., these will be detections outside the triangle but in the rectangle.
      
    } 
    
    #PLOT 10 - 2D image of semicircle - may be capped depending on depth 
    #plot(x_totpop,z_totpop,xlim=c(-maxdetrad,maxdetrad),ylim=c(0,maxdepth))
    #points(x_totpop[!vertr_totpop_edit>maxdetrad],z_totpop[!vertr_totpop_edit>maxdetrad],col="red")
    
    #PLOT 11 - Incidence angle - this will be the same if there is no bearing restriction
    #s3dcylinder<-scatterplot3d(x_totpop,y_totpop,z_totpop,cex.symbols=0.5,main="3D Scatterplot",xlim=c(-maxdetrad,maxdetrad),ylim=c(0,linelength),zlim=c(0,maxdepth))
    #s3dcylinder$points3d(x_totpop[!inctheta_totpop_edit>maxbearing],y_totpop[!inctheta_totpop_edit>maxbearing],z_totpop[!inctheta_totpop_edit>maxbearing],col="red")
    
    #PLOT 12 - Incidence angle
    #2D plot
    #plot(x_totpop[!inctheta_totpop_edit>maxbearing],z_totpop[!inctheta_totpop_edit>maxbearing],ylim=c(0,maxdepth))
    
    #PLOT 13 - all relevant observations
    #detailed plot
    #plot(x_totpop,z_totpop,xlim=c(-maxdetrad,maxdetrad),ylim=c(0,maxdepth))
    #points(x_totpop[!vertr_totpop_edit>maxdetrad],z_totpop[!vertr_totpop_edit>maxdetrad],col="red")
    #points(x_totpop[!inctheta_totpop_edit>maxbearing],z_totpop[!inctheta_totpop_edit>maxbearing],col="green")
    #points(x_totpop[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],z_totpop[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],col="purple") 
    #abline(v=c(-max_xy,max_xy)) 
    
    ###########################################################################
    #STEP 4 - Apply a detection function
    #A half normal detection function is being applied
    #500 - 2000 seem good values for sigma
    
    #For each radial distance, evaluate the detection function pdf and determine whether the animal was detected or not
    
    #create a place for results
    pdet_totpop<-mat.or.vec(length(vertr_totpop_edit), 1)
    det_totpop<-mat.or.vec(length(vertr_totpop_edit), 1)
    
    #apply detection function
    pdet_totpop<-grtheta(vertr_totpop_edit,det.func.params[1])
    
    for (i in 1:length(vertr_totpop_edit)){
      if (inctheta_totpop_edit[i]>maxbearing | vertr_totpop_edit[i]>maxdetrad) {
        det_totpop[i]==0} else
        {det_totpop[i]<-rbinom(1,1,pdet_totpop[i])}
    }
    
    #look at 3D plot - fixed instrument
    # scatterplot3d(x_totpop_edit[det_totpop==1],y_totpop_edit[det_totpop==1], z_totpop_edit[det_totpop==1],main="3D Scatterplot",xlim = c(-maxdetrad,maxdetrad),ylim = c(-maxdetrad,maxdetrad))
    
    # #look at depth vs x dimension
    # plot(x_totpop_edit[det_totpop==1],z_totpop_edit[det_totpop==1], xlim = c(-maxdetrad,maxdetrad),ylim = c(0,maxdepth))
    
    # #look at a birdseye view of the point transect
    # plot(x_totpop_edit[det_totpop==1],y_totpop_edit[det_totpop==1],xlim = c(-maxdetrad,maxdetrad), ylim = c(-maxdetrad,maxdetrad))
    
    if (data=="range"){
      #make the vertical radial distances of these detections the data for the likelihood
      vertr_obs<-vertr_totpop_edit[det_totpop==1] 
      results.obs[n]<-length(vertr_obs)} else {

        #note that in either bearing case, the code is written using the incidence angle
        #if in an analysis, a user specified that they had bearing.elev data, then these would simply be converted to incidence angles.
        #In either case, it would be assumed that the tilt of the instrument was known.
      
        verttheta_obs<-inctheta_totpop_edit[det_totpop==1]
  
        results.obs[n]<-length(verttheta_obs) }


    ############################################################################
    #STEP 5: Likelihood optimisation - using functions which need to be read in from "Depth.LH.functions.vs1.r"
    
    if (dist=="beta" & data=="range"){
      maxLH<-optimize(model.log.lik.fct.beta.r_obs,c(det.func.params[2],det.func.params[3]),shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,vertr_obs)}
    
    if (dist=="beta" && (data=="bearing.inc" | data=="bearing.elev")){
      maxLH<-optimize(model.log.lik.fct.beta.theta_obs,c(det.func.params[2],det.func.params[3]),shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,verttheta_obs)}
    
    if (dist=="uniform" & data=="range"){
      maxLH<-optimize(model.log.lik.fct.uniform.r_obs,c(det.func.params[2],det.func.params[3]),max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,vertr_obs)}
    
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="towed"){
      stop("A towed scenario with a uniform depth distribution and bearing data cannot be analysed")}
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="fixed"){
      maxLH<-optimize(model.log.lik.fct.uniform.theta_obs,c(det.func.params[2],det.func.params[3]),max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,verttheta_obs)}
    
    results.sig[n]<-maxLH$minimum
    
    
    #########################################################################
    #
    #STEP 6 - Estimate the average probability of detection.  This is the denominator of the LH function using the maximum LH value.  
    
    sigmaest<-as.numeric(maxLH[1])
    
    #Set up the integration for the denominator
    #############################################################################
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
      results.pdet[n]<-denom3pdet
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
      
      #Sum the results of denom2pdet to get the final integral.
      denom3pdet<-sum(denom2pdet)
      results.pdet[n]<-denom3pdet
    } 
    
    #########################################################################
    #STEP 8 - work up to estimating density
    #
    #corrected number of animals in the cylinder or rectangle.
    n.corrected<-results.obs[n]/results.pdet[n]
    if (type=="fixed"){
      results.est.cyl.max_xy[n]<-n.corrected}
    if (type=="towed"){
      results.est.rect.max_xy[n]<-n.corrected}
    print(n)
  } 
  
  #this is abundance, want density in the horizontal plane i.e., per km2
  if (type=="fixed"){
    results.dens<-results.est.cyl.max_xy/(pi*((max_xy/1000)^2))}
  if (type=="towed"){
    results.dens<-results.est.rect.max_xy/(2*(max_xy/1000)*(linelength/1000))}
  
  #percentage bias of true density
  meanpercentbias<-((mean(results.dens)-truedensity.km2)/truedensity.km2)*100
  #percentage bias of true sigma
  meanpercentbias_sigma<-((mean(results.sig)-det.func.params[1])/det.func.params[1])*100
  ############################################################
  #STEP 9 - Return useful output
  
  if (type=="fixed") {
    
    results.list<-list("results.cyl.mdr"=results.cyl.mdr,
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
        
        results.list<-list("results.rect.mdr"=results.rect.mdr,
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
  if (data=="range"){
    datatype<-"ranges."}
  if (data=="bearing.inc"){
    datatype<-"incident angles i.e., vertical is zero."}
  if (data=="bearing.elev"){
    datatype<-"elevation angles i.e., horizontal is zero."}
  print1a<-print(paste("The animal population metric set in this simulation was",popn[1]))
  print2<-print(paste("True density is",round(truedensity.km2,2),"animals per km2")) 
  print3<-print(paste("The collected data are", datatype))
  print4<-print(paste("The depth distribution is assumed to be",dist,"."))
  print5<-print(paste("The maximum detection radius was assumed to be", params[1],"m."))
  print6<-print(paste("The maximum depth was assumed to be", params[2],"m."))
  print7<-print(paste("The maximum vertical detection bearing was assumed to be ", round(params[3],2),"rad."))
  cat(c(print1a,print2,print3,print4,print5,print6,print7),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  if (type=="towed"){
    print8<-print(paste("The survey line length was assumed to be", params[4],"m."))
    cat(print8,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)}
    print9<-print(paste("The number of simulated animals was",totpop,"."))
    print10<-print(paste("The simulation was run", iterations, "times."))
    cat(c(print9,print10),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  print11a<-print(paste("The lower limit,upper limit and spread of the depth distribution parameters were", ll, ul, "and", spread,"."))
  cat(print11a,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  if (dist=="beta"){
    print11<-print(paste("The shape parameters of the beta distribution were", shape1, "and", shape2,"."))
    cat(print11,file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)}
    print12<-print(paste("The sigma value used in the half normal detection function was",det.func.params[1]," and the range of values used in the optimisation were",det.func.params[2], "and",det.func.params[3]))

  #Now print some results
  print13<-print(paste("The mean estimate of sigma was", round(mean(results.sig),2)))
  print14<-print(paste("The mean estimate of the detection probability was", round(mean(results.pdet),2)))
  print15<-print(paste("The mean density of animals was", round(mean(results.dens),2),"animals/km2"))
  print16<-print(paste("The range of estimated densities was between", round(min(results.dens)), "and" ,round(max(results.dens)),"animals/km2."))
  print17<-print(paste("The mean percentage bias in the density estimate was", round(meanpercentbias,2), "%."))
  print18<-print(paste("The mean percentage bias in the sigma estimate was", round(meanpercentbias_sigma,2), "%."))
  print19<-print(paste("The mean number of detections in the simulation was", round(mean(results.obs),2)))
  cat(c(print12,print13,print14,print15,print16,print17,print18,print19),file=paste("Sim_output_",label,".txt"),sep="\n",append=TRUE)
  
  #plot the most useful plot - the set up of the simulation
  if (type=="fixed"){
  
    jpeg(paste('Sim_plot_',label,'.jpg'))
    
    plot(x_totpop_edit,z_totpop_edit,xlim=c(-maxdetrad,maxdetrad),ylim=c(0,maxdepth),ylab=("Elevation (m)"),xlab=("Horizontal range (m)"))
    points(x_totpop_edit[!vertr_totpop_edit>maxdetrad],z_totpop_edit[!vertr_totpop_edit>maxdetrad],col="red")
    points(x_totpop_edit[!inctheta_totpop_edit>maxbearing],z_totpop_edit[!inctheta_totpop_edit>maxbearing],col="green")
    points(x_totpop_edit[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],z_totpop_edit[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],col="purple") 
    points(x_totpop_edit[det_totpop==1],z_totpop_edit[det_totpop==1],pch=16) 
    abline(v=c(-max_xy,max_xy))
    
    dev.off()
  }
  
  if (type=="towed"){
    
    jpeg(paste('Sim_plot_',label,'.jpg'))
    
    plot(x_totpop,z_totpop,xlim=c(-maxdetrad,maxdetrad),ylim=c(maxdepth,0),ylab=("Depth (m)"),xlab=("Horizontal range (m)"))
    points(x_totpop[!vertr_totpop_edit>maxdetrad],z_totpop[!vertr_totpop_edit>maxdetrad],col="red")
    points(x_totpop[!inctheta_totpop_edit>maxbearing],z_totpop[!inctheta_totpop_edit>maxbearing],col="green")
    points(x_totpop[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],z_totpop[!inctheta_totpop_edit>maxbearing & !vertr_totpop_edit>maxdetrad],col="purple") 
    points(x_totpop[det_totpop==1],z_totpop[det_totpop==1],pch=16)
    abline(v=c(-max_xy,max_xy))
    
    dev.off()
  }
  
  jpeg(paste('Sim_plot_',label,'_sig.vs.dens.jpg'))
  plot(results.sig,results.dens)
  dev.off()

  save(results.list, file = paste(label,".RData"))



}

