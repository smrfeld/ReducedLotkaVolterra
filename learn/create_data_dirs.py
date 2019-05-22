import os

parent = "data"
if not os.path.exists(parent):
	os.makedirs(parent)

dir_name = parent + "/diagnose_sampling"
if not os.path.exists(dir_name):
	os.makedirs(dir_name)

dir_name = parent + "/learn_centered"
if not os.path.exists(dir_name):
	os.makedirs(dir_name)

subdir_name = dir_name + "/moments"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/ixn_params"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/ixn_params_lp"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/diff_eq_rhs"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/sample_traj"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/sample_traj_lp"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)