#################### GET_TRAJ #######################################
# Function for extracting a trajectory in a dataset of trajectories #
#####################################################################

# Inputs : n = trajectory, data = trajectories dataset, id_traj = trajectories id
# Output : one trajectory

get_traj = function(n, data, id_traj){
  li = c(levels(as.factor(as.character(id_traj))))
  traj = subset(data, id_traj == li[n])
  return(traj)
}
