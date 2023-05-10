import numpy as np
from tqdm import tqdm
import os.path

file_name = "data/230313/230313_nd4_"

job_arr_start = 0
job_arr_end = 499
job_arr_step = 1
job_arr_l = np.arange(job_arr_start,job_arr_end+1,job_arr_step)

cnt = 0

for job_id in tqdm(job_arr_l):
    job_name = file_name + str(job_id)
    if os.path.isfile(job_name+"_scanx.csv") == False:
        print(job_id)
        continue

    if cnt == 0:
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

    cnt = cnt + 1

cmi_ave_l = cmi_ave_l / cnt
cmi_std_l = cmi_std_l / cnt

np.savez(file_name+"pp.npz",n_meas_l=n_meas_l,depth_l=depth_l,cmi_ave_l=cmi_ave_l,cmi_std_l=cmi_std_l)