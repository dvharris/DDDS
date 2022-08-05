multisensor_optimize_1<-function(type, data, dist, maxdetrads, maxdepths, maxbearings, max_xys,linelengths, popn, iterations, lls, uls, spreads,det.func.params,totpops,shape1,shape2){
  
  #create a results.post vector that will return the individual LH values and also create the results vectors that will be combined for complete results matrices for all iterations 
  
  results.all.obs<-list() 
  #going to use results.all.obs as a temporary store for other results.  These will be stored below the observations and can be transferred once the results.all.obs is returned.
  #Also for plots need to save the following:
#x_totpop_edit (or just x_totpop for towed situations)
#z_totpop_edit (or just z_totpop for towed situations)
#vertr_totpop_edit
#inctheta_totpop_edit
#det_totpop

#So for each instrument, will have 5 extra lists at end of results.all.obs (After other summary results too).  So for 5 instruments:
#Lists 1 - 5: values
#6 - 10: summary results for each instrument
#11 - 15: above results for first instrument (1 iteration of the simulation)
#etc for the other instruments.
  
  for (s in 1:length(maxdetrads)){
    
    #Need to assign each animal a location in the (x,z,y) dimension i.e., a radial distance with an azimuth and a depth
    
    x_totpop<-mat.or.vec(totpops[s], 1)
    z_totpop<-mat.or.vec(totpops[s], 1)
    y_totpop<-mat.or.vec(totpops[s], 1)
    
    for (i in 1:totpops[s]){
      x_totpop[i]<-runif(1, min=-maxdetrads[s], max=maxdetrads[s])
      if  (type=="towed"){
        y_totpop[i]<-runif(1, min=0, max=linelengths[s])} else 
          if (type=="fixed"){ 
            y_totpop[i]<-runif(1, min=-maxdetrads[s], max=maxdetrads[s])}
      
      #Note change to both beta and uniform depth distributions
      
      if (dist == "beta"){ 
        z_totpop[i]<-(rbeta(1, shape1, shape2)*spreads[s])+lls[s]} else       
          if (dist == "uniform"){
            z_totpop[i]<-runif(1, min=lls[s], max=uls[s]) } else {
              stop("Error - depth distribution not recognised")}
    }
    
    #Calculate the vertical radial distances and bearings for all animals - FIXED TRANSECTS ONLY
    
    if (type=="fixed"){
      #first need the horizontal radial distance for each animal
      #this is the relationship between x and y dimensions (use pythag)
      horizr_totpop<-sqrt((x_totpop^2)+(y_totpop^2))
      #NB:  don't want to use horizr values that are greater than the maxdetrad, as these are not in the population that I want to consider - don't want to pollute the detection process so need to remove them.
      horizr_totpop_edit<-horizr_totpop[!horizr_totpop>maxdetrads[s]]
      #edit other vectors too
      x_totpop_edit<-x_totpop[!horizr_totpop>maxdetrads[s]]
      y_totpop_edit<-y_totpop[!horizr_totpop>maxdetrads[s]]
      z_totpop_edit<-z_totpop[!horizr_totpop>maxdetrads[s]]
      
      #number animals in cylinder
      true.n.cylinder<-length(x_totpop_edit)
      #the line below is for when the true number changes for each iteration
      #results.all.mdr[s]<-true.n.cylinder
      results.all.obs[[s+length(maxdetrads)]]<-true.n.cylinder
      
      #now use the edited horizontal radial distance and the edited depth to estimate the vertical radial distance
      vertr_totpop_edit<-sqrt((horizr_totpop_edit^2)+(z_totpop_edit^2))
      #Note that these values include animals outside the detection area - any animal with a vertical radial distance greater than maxdetrad will have pdet of zero. 
      
      #check maximum vertical radial distance
      #Useful checks - left in for now
      #max_vertr<-sqrt((maxdepth^2)+(maxdetrad^2))
      #range(vertr_totpop_edit)
      
      #number animals in half sphere (which may be capped depending on depth)
      true.n.halfsphere<-length(x_totpop_edit[!vertr_totpop_edit>maxdetrads[s]])
      #the line below is for when the true number changes for each iteration
      #results.sph[s]<-true.n.halfsphere
      results.all.obs[[s+length(maxdetrads)]][2]<-true.n.halfsphere
      
      #now estimate the vertical bearing i.e., the elevation angle (we ignore the horizontal bearing - the detection function is symmetrical around the vertical axis)
      verttheta_totpop_edit<-atan(z_totpop_edit/horizr_totpop_edit)
      
      #now make these incident angles
      inctheta_totpop_edit<-(pi/2)-verttheta_totpop_edit
      #Note that these values include animals outside the detection area - any animal with a incidence angle greater than (pi/4) will have pdet of zero. 
      
      #number animals in cone - NB: this will be the same as the results.cyl.mdr/true cylinder if there is no bearing restriction
      
      true.n.cone<-length(x_totpop_edit[inctheta_totpop_edit<maxbearings[s]])
      #the line below is for when the true number changes for each iteration
      #results.cone[s]<-true.n.cone
      results.all.obs[[s+length(maxdetrads)]][3]<-true.n.cone
      
      #number in cylinder with radius of max_xy:
      x_totpop_maxxy<-x_totpop[!horizr_totpop>max_xys[s]]
      y_totpop_maxxy<-y_totpop[!horizr_totpop>max_xys[s]]
      z_totpop_maxxy<-z_totpop[!horizr_totpop>max_xys[s]]
      true.n.max_xy.cylinder<-length(x_totpop_maxxy)
      #results.cyl.maxxy[s]<-true.n.max_xy.cylinder
      results.all.obs[[s+length(maxdetrads)]][4]<-true.n.max_xy.cylinder
      
      #Note that these values include animals outside the detection area - any animal with a vertical bearing greater than maxbearing will have pdet of zero. I.e., these will be detections outside the cone but in the cylinder.
      
    }
    
    #Calculate the vertical radial distances and bearings for all animals - TOWED TRANSECTS ONLY
    
    if (type=="towed"){
      
      #number animals in rectangle
      true.n.rectangle<-length(x_totpop)
      #the line below is for when the true number changes for each iteration
      #results.rect.mdr[s]<-true.n.rectangle
      results.all.obs[[s+length(maxdetrads)]]<-true.n.rectangle
      
      #now use the horizontal radial distance and the depths to estimate the vertical radial distance.
      #NB: the horizontal distances did not have to be edited first, but still called "edited" to be consistent with the fixed case.
      vertr_totpop_edit<-sqrt((x_totpop^2)+(z_totpop^2))
      #Note that these values include animals outside the detection area - any animal with a vertical radial distance greater than maxdetrad will have pdet of zero. 
      
      #check maximum vertical radial distance
      #Useful checks - left in for now
      #max_vertr<-sqrt((maxdepth^2)+(maxdetrad^2))
      #range(vertr_totpop_edit)
      
      #number animals in semicircle (which may be capped depending on depth)
      true.n.semicircle<-length(x_totpop[!vertr_totpop_edit>maxdetrads[s]])
      #the line below is for when the true number changes for each iteration
      #results.semicircle[s]<-true.n.semicircle
      results.all.obs[[s+length(maxdetrads)]][2]<-true.n.semicircle
      
      #now estimate the vertical bearing i.e., the elevation angle 
      verttheta_totpop_edit<-atan(z_totpop/abs(x_totpop))
      
      #now make these incident angles
      inctheta_totpop_edit<-(pi/2)-verttheta_totpop_edit
      #Note that these values include animals outside the detection area - any animal with a incidence angle greater than (pi/4) will have pdet of zero. 
      
      #number animals in triangle - NB: this will be the same as the results.crect.mdr/true rectangle if there is no bearing restriction
      true.n.triangle<-length(x_totpop[inctheta_totpop_edit<maxbearings[s]])
      #the line below is for when the true number changes for each iteration
      #results.triangle[s]<-true.n.triangle
      results.all.obs[[s+length(maxdetrads)]][3]<-true.n.triangle
      
      #number in rectangle with radius of max_xy:
      x_totpop_maxxy<-x_totpop[!x_totpop>max_xys[s]]
      y_totpop_maxxy<-y_totpop[!x_totpop>max_xys[s]]
      z_totpop_maxxy<-z_totpop[!x_totpop>max_xys[s]]
      true.n.max_xy.rectangle<-length(x_totpop_maxxy)
      #results.rect.maxxy[s]<-true.n.max_xy.rectangle
      results.all.obs[[s+length(maxdetrads)]][4]<-true.n.max_xy.rectangle
      
      #Note that these values include animals outside the detection area - any animal with a vertical bearing greater than maxbearing will have pdet of zero. I.e., these will be detections outside the triangle but in the rectangle.
      
    } 
    
    ###########################################################################
    #STEP 2 - Apply a detection function
    #A half normal detection function is being applied
    #500 - 2000 seem good values for sigma
    
    #For each radial distance, evaluate the detection function pdf and determine whether the animal was detected or not
    
    #create a place for results
    pdet_totpop<-mat.or.vec(length(vertr_totpop_edit), 1)
    det_totpop<-mat.or.vec(length(vertr_totpop_edit), 1)
    
    #apply detection function
    pdet_totpop<-grtheta(vertr_totpop_edit,det.func.params[1])
    
    for (i in 1:length(vertr_totpop_edit)){
      if (inctheta_totpop_edit[i]>maxbearings[s] | vertr_totpop_edit[i]>maxdetrads[s]) {
        det_totpop[i]==0} else
        {det_totpop[i]<-rbinom(1,1,pdet_totpop[i])}
    }
    
    
    if (data=="range"){
      #make the vertical radial distances of these detections the data for the likelihood
      vertr_obs<-vertr_totpop_edit[det_totpop==1]
      results.all.obs[[s]]<-vertr_obs} else {
        #note that in either bearing case, the code is written using the incidence angle
        #if in an analysis, a user specified that they had bearing.elev data, then these would simply be converted to incidence angles.
        #In either case, it would be assumed that the tilt of the instrument was known.
        #verttheta_obs<-verttheta_totpop_edit[det_totpop==1]
        verttheta_obs<-inctheta_totpop_edit[det_totpop==1]
        #incidence angle should be used - edit made Oct14 by Danielle
        
        results.all.obs[[s]]<-verttheta_obs}
  
  #here want to save the coordinate infomration from the last iteration of the simulation as an illustration of the simulation.
  #this is done for each instrument
  #this information is saved for every iteration but is only used in the last iteration to produce the plots
  if (type=="towed"){
  results.all.obs[[((s-1)*5)+1+(2*length(maxdetrads))]]<-x_totpop
  results.all.obs[[((s-1)*5)+2+(2*length(maxdetrads))]]<-z_totpop
  results.all.obs[[((s-1)*5)+3+(2*length(maxdetrads))]]<-vertr_totpop_edit
  results.all.obs[[((s-1)*5)+4+(2*length(maxdetrads))]]<-inctheta_totpop_edit
  results.all.obs[[((s-1)*5)+5+(2*length(maxdetrads))]]<-det_totpop}
  
  if (type=="fixed"){
  results.all.obs[[((s-1)*5)+1+(2*length(maxdetrads))]]<-x_totpop_edit
  results.all.obs[[((s-1)*5)+2+(2*length(maxdetrads))]]<-z_totpop_edit
  results.all.obs[[((s-1)*5)+3+(2*length(maxdetrads))]]<-vertr_totpop_edit
  results.all.obs[[((s-1)*5)+4+(2*length(maxdetrads))]]<-inctheta_totpop_edit
  results.all.obs[[((s-1)*5)+5+(2*length(maxdetrads))]]<-det_totpop}
  
  }
  return(results.all.obs)
} 

multisensor_optimize_2<-function(p,type, data, dist, maxdetrads, maxdepths, maxbearings, max_xys,linelengths, popn, iterations, lls, uls, spreads,det.func.params,results.all.obs,shape1,shape2){   
 #NB det.func.params not required for multisensor_optimize_2 - will need to remove from all calls of this function
 
 
  results.post<-array(0,length(maxdepths))
  
  for (s in 1:length(maxdepths)){
    ############################################################################
    #STEP 3: Likelihood optimisation - using functions which need to be read in from "Functions_Aug14.r"
    
    if (dist=="beta" & data=="range"){
vertr_obs<-results.all.obs[[s]]
      results.post[s]<-model.log.lik.fct.beta.r_obs(p,shape1,shape2,lls[s],uls[s],max_xys[s],maxdetrads[s],maxbearings[s],maxdepths[s],dist,type,vertr_obs)}
  
    if (dist=="beta" && (data=="bearing.inc" | data=="bearing.elev")){
      verttheta_obs<-results.all.obs[[s]]
      results.post[s]<-model.log.lik.fct.beta.theta_obs(p,shape1,shape2,lls[s],uls[s],max_xys[s],maxdetrads[s],maxbearings[s],maxdepths[s],dist,type,verttheta_obs)}
    
    if (dist=="uniform" & data=="range"){
      vertr_obs<-results.all.obs[[s]]
      results.post[s]<-model.log.lik.fct.uniform.r_obs(p,max_xys[s],maxdetrads[s],maxbearings[s],maxdepths[s],lls[s],uls[s],spreads[s],dist,type,vertr_obs)}
    
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="towed"){
      stop("A towed scenario with a uniform depth distribution and bearing data cannot be analysed")}
    if (dist=="uniform" && (data=="bearing.inc" | data=="bearing.elev") && type=="fixed"){
      verttheta_obs<-results.all.obs[[s]]
      results.post[s]<-model.log.lik.fct.uniform.theta_obs(p,max_xys[s],maxdetrads[s],maxbearings[s],maxdepths[s],lls[s],uls[s],spreads[s],dist,type,verttheta_obs)}
  }
 
  return(sum(results.post))
}
