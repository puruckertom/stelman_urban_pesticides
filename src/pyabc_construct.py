#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import jn_setup
# from simulation_procedure import model, mode
from tools.paths import *
from simulation_procedure import make_model
import pandas as pd
import hydroeval as he
import numpy as np
import pyabc
import uuid

mode = "test"
model = make_model(mode = mode, swmm_cleanup = 'none', vvwm_cleanup = 'none')
# In[2]:


"""
Priors. Get the values from that csv. Just for SWMM at first.
"""
swmm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])

'''
Link up with the vvwm priors and make one big list with 36 params
'''
vvwm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges_vvwm.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])

param_ranges = pd.concat([swmm_ranges, vvwm_ranges], axis = 0)

if mode == "debug":
    param_ranges = param_ranges.loc[['NImperv','kd']]

priors = param_ranges.to_dict("index")

# borrowed from Jeff: <https://github.com/JeffreyMinucci/bee_neonic_abc/blob/master/pyabc_run.ipynb>
prior = pyabc.Distribution(**{key: pyabc.RV("uniform", loc = v['Min'], scale = v['Range']) for key, v in priors.items()})


# In[3]:


# priors


# In[4]:


# prior


# # Make the .new object
# ### 1. Import observed data

# In[5]:


# import it again to make inspect it
# specifically for TEST mode!
if mode == 'debug':
    with open(os.path.join(main_path, 'master_debug','debug_obs_data.txt'),'r') as read_file:
        obs_dict = eval(read_file.read())
elif mode == 'test':
    with open(os.path.join(main_path, 'master_test','test_obs_data.txt'),'r') as read_file:
        obs_dict = eval(read_file.read())
elif mode == 'run':
    with open(os.path.join(main_path, 'master','obs_data.txt'),'r') as read_file:
        obs_dict = eval(read_file.read())
# obs_dict


# ### 2. Initialize dask client for dask distributed sampler

# In[6]:


# from dask.distributed import Client#, LocalCluster
# cluster = LocalCluster()#n_workers=(90/2), threads_per_worker = 2)  # Set for 96 vCPU compute instance
# client = Client(cluster)#,timeout=400)

# make it simpler
# if __name__ == "__main__":
# client = Client()

# sampler = pyabc.sampler.DaskDistributedSampler(dask_client = client)

# See if this takes all those errors out
sampler = pyabc.sampler.SingleCoreSampler()
# make the process more transparent
sampler.sample_factory.record_rejected = True
sampler.show_progress = True

# ### 3. Set up a sqlite db directory

# In[7]:


# Initialize a new ABC inference run
dbid = uuid.uuid4().hex[0:8]
print(dbid)
database_dir = os.path.join(temp_path, 'results_db')  
if not os.path.exists(database_dir):
    os.mkdir(database_dir)
db_path = ("sqlite:///" + os.path.join(database_dir, "test_pyabc_" + dbid + ".db"))


# ### 4. Defining a Distance function
# 
# We need to refactor the NSE distance function using the pyabc.Distance class.
# We will need the hydroeval library and the pyabc.SimpleFunctionDistance to do this

# In[8]:


# make a file to hold onto these NSEs for our own record
with open(os.path.join(temp_path, "NSEs_" + dbid + ".txt"), "w") as nse_file:
    nse_file.write("NSEs\n")


# In[9]:


def nse(x, x_0):
    nse = he.evaluator(he.nse, simulation_s = np.array(list(x.values())), evaluation = np.array(list(x_0.values())))[0]
    print("nse ", nse)
    # make record
    with open(os.path.join(temp_path, "NSEs_" + dbid + ".txt"),"a") as nse_file:
        nse_file.write(str(nse)+"\n")
    return nse
    
NSE = pyabc.SimpleFunctionDistance(fun=nse)

# the best answer is 1
# make one that measures distance from 1
NSED = pyabc.SimpleFunctionDistance(fun = lambda x, x_0: 1 - nse(x, x_0))


# ### 5. Define ABCSMC object

# In[10]:


abc = pyabc.ABCSMC(model, prior, population_size = pyabc.ConstantPopulationSize(40), sampler = sampler, distance_function = NSED)
                #    # might fix the dask problem too
                #    population_size = pyabc.ConstantPopulationSize(4), # just to shorten the run
                #    sampler = sampler,
                #    distance_function = NSED)


# In[11]:


# abc


# ### 6. Initialize a new abc run

# In[12]:


abc.new(db_path, obs_dict)


# In[13]:


# add 1 to generations
history = abc.run(max_nr_populations=2, minimum_epsilon=0.2)

