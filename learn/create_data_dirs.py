import os

# Dirs
dir_names = ["data","data/learn_params","data/learn_params/ixn_params","data/learn_params/moments"]

for dir_name in dir_names:
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

dir_name = "data/sample_traj"
if not os.path.exists(dir_name):
	os.makedirs(dir_name)

for i in range(0,500):
	dir_name = "data/sample_traj/%04d"%i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

	dir_name = "data/sample_traj/%04d/lattices"%i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

	dir_name = "data/sample_traj/%04d/moments"%i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

	dir_name = "data/sample_traj/%04d/covs"%i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)