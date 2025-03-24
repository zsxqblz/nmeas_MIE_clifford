import numpy as np
from tqdm import tqdm
import os.path

file_name = "data/250322/250322_n304_"

job_arr_start = 0
job_arr_end = 99
job_arr_step = 1
job_arr_l = np.arange(job_arr_start,job_arr_end+1,job_arr_step)

# dist_ave_l = []
# dist_std_l = []
cnt = 0
for job_id in tqdm(job_arr_l):
    job_name = file_name + str(job_id)
    if os.path.isfile(job_name) == False:
        print(job_name)
        continue

    if cnt == 0:
        data_raw =  np.loadtxt(job_name,skiprows=1,delimiter=",")
        dy_l = np.array(data_raw[:,0])
        dist_ave_l = np.array(data_raw[:,1])
        dist_std_l = np.array(data_raw[:,2])

    else:
        data_raw =  np.loadtxt(job_name,skiprows=1,delimiter=",")
        dist_ave_l = dist_ave_l + data_raw[:,1]
        dist_std_l = dist_std_l + data_raw[:,2]

    cnt = cnt + 1

dist_ave_l = dist_ave_l / cnt
dist_std_l = dist_std_l / cnt

np.savez(file_name+"pp.npz",dy_l=dy_l,dist_ave_l=dist_ave_l,dist_std_l=dist_std_l)