
#Adjusted from vs1 by Danielle Harris in Dec14.  The adjustments have decreased the run time of the calculations.
#Adjusted from vs2 by Danielle Harris in Mar16.  Warning messages have been introduced if there are problems with the entry parameters relating to the maxdetrad, maxbearing, ll and ul. 
############################################################################
#(1)Loglikelihood function for beta depth distribution and range data
model.log.lik.fct.beta.r_obs<-function(p,shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,vertr_obs){
    #denote that the parameter that we want to optimise is sigma
    sigma<-p

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
denom1<-array(0,dim(meanrthetagrid)[2])
    for (i in 1:dim(meanrthetagrid)[2]){ 
        #evaluate the functions        
        if (type=="fixed"){
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else  
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
    if (type=="towed"){
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
        if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_towed(meanrthetagrid[1,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
      }    
#Multiply the results by the length of each range and theta cell.
denom2<-denom1*rrange[2]*thetarange[2]
    
#Sum the results of denom2 to get the final integral.
denom3<-sum(denom2)
}

#############################################################################################
    #create a results place for each numerator and likelihood value.  NB: the denominator value will be constant for all iterations 
    num<-array(NA,length(vertr_obs))
    loglik<-array(NA,length(vertr_obs))
    
    #calculate the numerator
    for (m in 1:length(vertr_obs)){     
      num1<-array(0,length(meant.vec))
      num2<-c()
       for (q in 1:length(meant.vec)){
      #evaluate the functions
      if (type=="fixed"){
        if (cos(meant.vec[q])*vertr_obs[m] > maxdepth) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] > ul) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] < ll) {num1[q]<-0} else 
        num1[q]<-grtheta(vertr_obs[m],sigma)*rxdist_fixed(vertr_obs[m],meant.vec[q],max_xy)*betascalepdf_rtheta(vertr_obs[m],meant.vec[q],shape1,shape2,ll,ul)}
      if (type=="towed"){
        if (cos(meant.vec[q])*vertr_obs[m] > maxdepth) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] > ul) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(vertr_obs[m],sigma)*rxdist_towed(vertr_obs[m],max_xy)*betascalepdf_rtheta(vertr_obs[m],meant.vec[q],shape1,shape2,ll,ul)}
      }
       #Multiply the results by the length of each range cell.
      num2<-num1*thetarange[2]
      #Sum the results of num2 to get the final integral.
      num3<-sum(num2)
      num[m]<-num3 
      loglik[m]<-log(num3)-log(denom3) 
      #print(m)
    }

    post<-sum(loglik)
    
    #if num1 is completely zero then all candidate whale positions for a given observation are not valid (for a given sigma value) and the whole loglik answer for that observation and sigma value will be -Inf.  Therefore, post will also be -Inf.  One solution would be to alter the distribution range of the whales.  Print a warning message if this occurs.
    if(post != "NaN" && post==-Inf){
      warning("Under certain tested values for sigma and parameters of maxdepth, maxdetrad, maxbearing, ll and ul, certain observed values are generating no possible positions for the calling animal.  In particular, check that the maxdetrad and maxbearing parameters are compatible.")}
    
    return(-post)
  } 
 
 ############################################################################
#(2) Loglikelihood function for beta depth distribution and bearing data 
 model.log.lik.fct.beta.theta_obs<-function(p,shape1,shape2,ll,ul,max_xy,maxdetrad,maxbearing,maxdepth,dist,type,verttheta_obs){
    #denote that the parameter that we want to optimise is sigma
    sigma<-p

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
  denom1<-array(0,dim(meanrthetagrid)[2])
      for (i in 1:dim(meanrthetagrid)[2]){  
        #evaluate the functions
        if (type=="fixed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
        if (type=="towed"){
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
          if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_towed(meanrthetagrid[1,i],max_xy)*betascalepdf_rtheta(meanrthetagrid[1,i],meanrthetagrid[2,i],shape1,shape2,ll,ul)}}
      }
       
#Multiply the results by the length of each range and theta cell.
denom2<-denom1*rrange[2]*thetarange[2]
    
#Sum the results of denom2 to get the final integral.
denom3<-sum(denom2)
}

#############################################################################################
#create a results place for each numerator and likelihood value.  NB: the denominator value will be constant for all iterations 
    num<-array(NA,length(verttheta_obs))
    loglik<-array(NA,length(verttheta_obs))
    
    #calculate the numerator
    
    for (m in 1:length(verttheta_obs)){     
      num1<-array(0,length(meanr.vec))
      num2<-c()
       for (q in 1:length(meanr.vec)){
      #evaluate the functions
      if (type=="fixed"){
        if (cos(verttheta_obs[m])*meanr.vec[q] > maxdepth) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] > ul) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] < ll) {num1[q]<-0} else 
         num1[q]<-grtheta(meanr.vec[q],sigma)*rxdist_fixed(meanr.vec[q],verttheta_obs[m],max_xy)*betascalepdf_rtheta(meanr.vec[q],verttheta_obs[m],shape1,shape2,ll,ul)}
      if (type=="towed"){
        if (cos(verttheta_obs[m])*meanr.vec[q] > maxdepth) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] > ul) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(meanr.vec[q],sigma)*rxdist_towed(meanr.vec[q],max_xy)*betascalepdf_rtheta(meanr.vec[q],verttheta_obs[m],shape1,shape2,ll,ul)}
       }

       #Multiply the results by the length of each range cell.
      num2<-num1*rrange[2]
      #Sum the results of num2 to get the final integral.
      num3<-sum(num2)
      num[m]<-num3 
      loglik[m]<-log(num3)-log(denom3) 
      #print(m)
    }

    post<-sum(loglik)
    
    #if num1 is completely zero then all candidate whale positions for a given observation are not valid (for a given sigma value) and the whole loglik answer for that observation and sigma value will be -Inf.  Therefore, post will also be -Inf.  One solution would be to alter the distribution range of the whales.  Print a warning message if this occurs.
    if(post != "NaN" && post==-Inf){
      warning("Under certain tested values for sigma and parameters of maxdepth, maxdetrad, maxbearing, ll and ul, certain observed values are generating no possible positions for the calling animal.  In particular, check that the maxdetrad and maxbearing parameters are compatible.")}
    
    return(-post)
  } 
  
############################################################################
#(3) Loglikelihood function for uniform depth distribution and range data
model.log.lik.fct.uniform.r_obs<-function(p,max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,vertr_obs){
    #denote that the parameter that we want to optimise is sigma
    sigma<-p

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
#calculate the denominator for either scenario with a uniform distribution

if (dist=="uniform"){   
  denom1<-array(0,dim(meanrthetagrid)[2])
  for (i in 1:dim(meanrthetagrid)[2]){  
    #evaluate the functions
    if (type=="fixed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*uniformpdf_rtheta(spread)}}
    if (type=="towed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_towed(meanrthetagrid[1,i],max_xy)*uniformpdf_rtheta(spread)}}
      }
    
    
#Multiply the results by the length of each range and theta cell.
denom2<-denom1*rrange[2]*thetarange[2]
    
#Sum the results of denom2 to get the final integral.
denom3<-sum(denom2)
}
#############################################################################################

#create a results place for each numerator and likelihood value.  NB: the denominator value will be constant for all iterations 
    num<-array(NA,length(vertr_obs))
    loglik<-array(NA,length(vertr_obs))
    
    #calculate the numerator
    
    for (m in 1:length(vertr_obs)){     
      num1<-array(0,length(meant.vec))
      num2<-c()
       for (q in 1:length(meant.vec)){ 
      #evaluate the functions
      if (type=="fixed"){
        if (cos(meant.vec[q])*vertr_obs[m] > maxdepth) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] > ul) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(vertr_obs[m],sigma)*rxdist_fixed(vertr_obs[m],meant.vec[q],max_xy)*uniformpdf_rtheta(spread)}
      if (type=="towed"){
        if (cos(meant.vec[q])*vertr_obs[m] > maxdepth) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] > ul) {num1[q]<-0} else 
        if (cos(meant.vec[q])*vertr_obs[m] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(vertr_obs[m],sigma)*rxdist_towed(vertr_obs[m],max_xy)*uniformpdf_rtheta(spread)}
      }
       #Multiply the results by the length of each range cell.
      num2<-num1*thetarange[2]
      #Sum the results of num2 to get the final integral.
      num3<-sum(num2)
      num[m]<-num3 
      loglik[m]<-log(num3)-log(denom3) 
      #print(m)
    }

    post<-sum(loglik)
    
    #if num1 is completely zero then all candidate whale positions for a given observation are not valid (for a given sigma value) and the whole loglik answer for that observation and sigma value will be -Inf.  Therefore, post will also be -Inf.  One solution would be to alter the distribution range of the whales.  Print a warning message if this occurs.
    if(post != "NaN" && post==-Inf){
      warning("Under certain tested values for sigma and parameters of maxdepth, maxdetrad, maxbearing, ll and ul, certain observed values are generating no possible positions for the calling animal.  In particular, check that the maxdetrad and maxbearing parameters are compatible.")}
    
    return(-post)
  } 
 
############################################################################
#(4) Loglikelihood function for uniform depth distribution and bearing data 
 model.log.lik.fct.uniform.theta_obs<-function(p,max_xy,maxdetrad,maxbearing,maxdepth,ll,ul,spread,dist,type,verttheta_obs){
    #denote that the parameter that we want to optimise is sigma
    sigma<-p

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
#calculate the denominator for either scenario with a uniform distribution

if (dist=="uniform"){   
  denom1<-array(0,dim(meanrthetagrid)[2])
  for (i in 1:dim(meanrthetagrid)[2]){  
    #evaluate the functions
    if (type=="fixed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_fixed(meanrthetagrid[1,i],meanrthetagrid[2,i],max_xy)*uniformpdf_rtheta(spread)}}
    if (type=="towed"){
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > maxdepth) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] > ul) {denom1[i]<-0} else 
      if (cos(meanrthetagrid[2,i])*meanrthetagrid[1,i] < ll) {denom1[i]<-0} else
        {denom1[i]<-grtheta(meanrthetagrid[1,i],sigma)*rxdist_towed(meanrthetagrid[1,i],max_xy)*uniformpdf_rtheta(spread)}}
      }
  
    
#Multiply the results by the length of each range and theta cell.
denom2<-denom1*rrange[2]*thetarange[2]
    
#Sum the results of denom2 to get the final integral.
denom3<-sum(denom2)
}
#############################################################################################

#create a results place for each numerator and likelihood value.  NB: the denominator value will be constant for all iterations 
    num<-array(NA,length(verttheta_obs))
    loglik<-array(NA,length(verttheta_obs))
    
    #calculate the numerator
    
    for (m in 1:length(verttheta_obs)){     
      num1<-array(0,length(meanr.vec))
      num2<-c()
       for (q in 1:length(meanr.vec)){
      #evaluate the functions
      if (type=="fixed"){
        if (cos(verttheta_obs[m])*meanr.vec[q] > maxdepth) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] > ul) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(meanr.vec[q],sigma)*rxdist_fixed(meanr.vec[q],verttheta_obs[m],max_xy)*uniformpdf_rtheta(spread)}
      if (type=="towed"){
        if (cos(verttheta_obs[m])*meanr.vec[q] > maxdepth) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] > ul) {num1[q]<-0} else 
        if (cos(verttheta_obs[m])*meanr.vec[q] < ll) {num1[q]<-0} else
         num1[q]<-grtheta(meanr.vec[q],sigma)*rxdist_towed(meanr.vec[q],max_xy)*uniformpdf_rtheta(spread)}
      }
       #Multiply the results by the length of each range cell.
      num2<-num1*rrange[2]
      #Sum the results of num2 to get the final integral.
      num3<-sum(num2)
      num[m]<-num3 
      loglik[m]<-log(num3)-log(denom3) 
      #print(m)
    }

    post<-sum(loglik)
    
    #if num1 is completely zero then all candidate whale positions for a given observation are not valid (for a given sigma value) and the whole loglik answer for that observation and sigma value will be -Inf.  Therefore, post will also be -Inf.  One solution would be to alter the distribution range of the whales.  Print a warning message if this occurs.
    if(post != "NaN" && post==-Inf){
      warning("Under certain tested values for sigma and parameters of maxdepth, maxdetrad, maxbearing, ll and ul, certain observed values are generating no possible positions for the calling animal.  In particular, check that the maxdetrad and maxbearing parameters are compatible.")}
    
    return(-post)
  }   