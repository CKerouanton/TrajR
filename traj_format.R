########################################### GTRAJ ##############################################################################
# Largely inspired by adehabitatLT (Callenge, 2006), but compute and assign trajectory metrics between i-1 and i, on i.
# Metrics added : cumuled distance, cumuled timestamp.
# Angles are used with a cos transformation, in order to visualize turning backs in the trajectory.
# Those new metrics are useful in order to use pre-processing functions on trajectories (outliers detection, Douglas-Peucker), and in order to use ST-DBSCAN for stop detection


# Ref = Calenge, C. (2006). The package “adehabitat” for the R software: a tool for the analysis of space and habitat use by animals. Ecological modelling, 197(3), 516-519.

# Inputs :
# points = spatio-temporal traj with points id within the trajectory (idp), points id within the whole dataset of trajectories (idpall), trajectory id (idt)
# xy = matrix with longitude and latitude
# date = list with trajectory dates (character or POSIXct)
# angle = logical. Computation of points angles, but works with adehabitatLT

# Outputs :
# traj = new trajectory with computed metrics (x, y, date, dt = duration, dist = distance, r2n = distance to first point, sp = speed, angle, tcum = cumulated timestamp, dcum = cumulated timestamp)
#################################################################################################################################

gtraj = function(points, xy, date, angle = TRUE){
  
  traj = data.frame(idp = points$idp, idt = as.character(points$idt), idpall = points$idpall)
  
  traj$idt = as.character(traj$idt)
  
  idt1 = traj$idt[1]
  
  traj$x = xy[,1]
  traj$y = xy[,2]
  traj$date = as.POSIXct(date)
  
  traj$dt = c(0, unlist(lapply(c(2:nrow(traj)), function(x){as.numeric(difftime(traj$date[x], traj$date[x-1], units = "secs"))})))
  
  traj = traj[is.na(traj$x)!= TRUE,]
  
  trajb = spatial_traj93(traj, traj$x, traj$y)
  
  traj$dist = c(0, unlist(lapply(c(2:nrow(traj)), function(i){gDistance(trajb[i,], trajb[i-1,])})))
  
  traj$r2n = c(0, unlist(lapply(c(2:nrow(traj)), function(i){gDistance(trajb[i,], trajb[1,])})))
  
  rm(trajb)
  
  traj$sp = c(0, traj$dist[2:nrow(traj)]/traj$dt[2:nrow(traj)])
  
  if (angle){
  traj$angle = ltraj2spdf(as.ltraj(xy, date, idt1))$rel.angle
  traj$angle = cos(traj$angle)
  } else{}
  
  traj$tcum = cumsum(traj$dt)
  
  traj$dcum = cumsum(traj$dist)
  
  return(traj)
}





################################## GTRAJ2 #########################################################################
# Function to compute new metrics, relative to terrain context
# It works with trajectories computed with GTRAJ, and it needs the function realdistR, available in my repositery "realdistR"
#
# Inputs :
# traj = trajectory as a result of GTRAJ
# mnt = altitudes, as a spatial grid data frame
# seg = resolution for altitude extraction, advice : seg > mnt resolution
# deniv = logical. Compute deniveles from a point i-1 to i.
#
# Output :
# trajectory with metrics : 
# distr = real distance, using mnt
# deniv = denivele between two points
# spr = speed using real distance
# spd = denivele speed
##############################################################################################################################

gtraj2 = function(traj, mnt, seg, deniv = TRUE){
  traj$distr = 0
  traj$deniv = 0
  traj$spr = 0
  if (deniv){traj$spd = 0}
  for (i in 2:nrow(traj)){
  dd = realdistR(traj[i,], traj[i-1,], reso = seg, alti_grid)
  traj$distr[i] = dd$distr
  traj$spr[i] = traj$distr[i]/traj$dt[i]
  
  if (deniv){
  traj$deniv[i] = dd$deniv
  traj$spd[i] = traj$deniv[i]/traj$dt[i]
  }
  }
  return(traj)
}



