##############################################################################
#Code the detection function and depth distribution and the r divided by the x distribution (the r comes from the change of variables)

#Detection function
grtheta=function(r, sigma){
  #Inputs:
  #       r - radial distances
  #       sigma - shape parameter
  #-----------------------------------------------------------------
  exp((-r^2)/(2*sigma^2))}

betascalepdf_rtheta<-function(r,theta,shape1,shape2,ll,ul){
  ((((r*cos(theta))-ll)^(shape1-1))*((ul-(r*cos(theta)))^(shape2-1)))/((beta(shape1,shape2))*((ul-ll)^(shape1+shape2-1)))
}

uniformpdf_rtheta<-function(spread){
  1/spread}

rxdist_fixed<-function(r,theta,max_xy){  #r multiplied by the x distribution 
  #Inputs:
  #     r - radial distance
  #     verttheta - vertical bearing
  r*(2*(r*sin(theta)))/max_xy^2}

rxdist_towed<-function(r,max_xy){  #r multiplied by the x distribution 
  #Inputs:
  #     r - radial distance
  r/max_xy}

