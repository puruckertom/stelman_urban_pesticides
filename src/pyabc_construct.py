#!/usr/bin/env python
# coding: utf-8

# In[1]:


from simulation_procedure import model
from tools.paths import *
import pandas as pd, pyabc, hydroeval


# In[2]:


"""
Priors. Get the values from that csv. Just for SWMM at first.
"""
swmm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges.csv"), index_col=0,
                           usecols = ["Parameter","Min", "Range"])

'''
Link up with the vvwm priors and make one big list with 36 params
'''
vvwm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges_vvwm.csv"), index_col=0,
                           usecols = ["Parameter","Min", "Range"])

param_ranges = pd.concat([swmm_ranges, vvwm_ranges], axis = 0)

priors = param_ranges.to_dict("index")

# borrowed from Jeff: <https://github.com/JeffreyMinucci/bee_neonic_abc/blob/master/pyabc_run.ipynb>
prior = pyabc.Distribution(**{key: pyabc.RV("uniform", loc = v['Min'], scale = v['Range'])
                        for key, v in priors.items()})


# In[3]:


prior


# # Make the .new object
# ### 1. Import observed data

# In[4]:


# import it again to make inspect it
# specifically for TEST mode!
with open(os.path.join(main_path, 'master_test','test_obs_data.txt'),'r') as read_file:
    obs_dict = eval(read_file.read())
# obs_dict


# ### 2. Initialize dask client for dask distributed sampler

# In[5]:


from dask.distributed import Client, LocalCluster
cluster = LocalCluster()#n_workers=(90/2), threads_per_worker = 2)  # Set for 96 vCPU compute instance
client = Client(cluster)#,timeout=400)

sampler = pyabc.sampler.DaskDistributedSampler(dask_client = client)


# ### 3. Define ABCSMC object

# In[7]:


abc = pyabc.ABCSMC(model, prior, #distance_function = hydroeval.nse, 
                   population_size = pyabc.AdaptivePopulationSize(40, max_population_size = 40), # just to shorten the run
                   sampler = sampler,
                   distance_function = pyabc.PNormDistance(p=2))


# In[8]:


abc


# ### 4. Set up a sqlite db directory

# In[10]:


# Initialize a new ABC inference run
database_dir = os.path.join(temp_path, 'results_db')  
if not os.path.exists(database_dir):
    os.mkdir(database_dir)
db_path = ("sqlite:///" +
           os.path.join(database_dir, "test_pyabc.db"))


# ### 5. Initialize a new abc run

# In[11]:


abc.new(db_path, obs_dict)


# In[12]:


history = abc.run(max_nr_populations=3, minimum_epsilon=0.2)

