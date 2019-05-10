import os

# Dirs
dir_names = ["data","data/learn_params","data/learn_params/ixn_params","data/learn_params/moments"]

for dir_name in dir_names:
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)