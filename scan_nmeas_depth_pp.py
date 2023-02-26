import numpy as np
from tqdm import tqdm

file_name = "data/230226/230226_nd2_"

job_arr_start = 0
job_arr_end = 99
job_arr_step = 1
job_arr_l = np.arange(job_arr_start,job_arr_end+1,job_arr_step)

for job_id in tqdm(job_arr_l):
    job_name = file_name + str(job_id)
    if job_id == job_arr_start:
        n_meas_l = np.loadtxt(job_name+"_scanx.csv",skiprows=1,delimiter=",")
        depth_l = np.loadtxt(job_name+"_scany.csv",skiprows=1,delimiter=",")
        data_raw =  np.loadtxt(job_name+"_data.csv",skiprows=1,delimiter=",")
        n_meas_length = len(n_meas_l)
        depth_length = len(depth_l)
        cmi_ave_l = data_raw[:,0].reshape((depth_length,n_meas_length))
        cmi_std_l = data_raw[:,1].reshape((depth_length,n_meas_length))

    else:
        data_raw =  np.loadtxt(job_name+"_data.csv",skiprows=1,delimiter=",")
        cmi_ave_l = cmi_ave_l + data_raw[:,0].reshape((depth_length,n_meas_length))
        cmi_std_l = cmi_std_l + data_raw[:,1].reshape((depth_length,n_meas_length))

cmi_ave_l = cmi_ave_l / len(job_arr_l)
cmi_std_l = cmi_std_l / len(job_arr_l)

np.savez(file_name+"pp.npz",n_meas_l=n_meas_l,depth_l=depth_l,cmi_ave_l=cmi_ave_l,cmi_std_l=cmi_std_l)