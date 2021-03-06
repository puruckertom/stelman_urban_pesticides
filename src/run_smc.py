import sys, argparse
from simulation_procedure import make_model
from tools.paths import *
import pandas as pd
import hydroeval as he
import numpy as np
import pyabc
import uuid

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode", help = """
    Specify if you want to use 'debug', 'test', or 'run' mode:\n
    Debug mode is designed to test code quickly to expose bugs in runs. 
    It only uses 2 of the 36 paramters, 1 of the 7 outfalls, 3 of the 106 
    summary statistics, and a 103-day-subset of the 3287-day timeframe of the data.\n
    Test mode is designed to impersonate a minimally time-consuming, but complete, 
    run of the simulation procedure by using a small sample of the data.\n
    Run mode is what you want to run when the code is fully baked, 
    and you're ready to get results based on the whole dataset.""", required=True, type = str)
    parser.add_argument("-c", "--cleanup", nargs = 2, metavar = ("SWMMcleanup", "VVWMcleanup"), help = """
    Specify SWMM and VVWM cleanup levels, each as 'full', 'some', or 'none':\n
    * SWMMcleanup:\n 
    - If full, each swmm-related file will be deleted after fulfilling its use to save memory.\n 
    - If some, executable and binary swmm-related files will be deleted after use, but human-readable files won't.\n 
    - If none, only the run-specific copy of the master executable (dll) file will be deleted after use.\n
    * VVWMcleanup:\n 
    - If full, each vvwm-related file will be deleted after fulfilling its use to save memory.\n 
    - If some, all vvwm-related files but one, the human-readable daily output csv file, will be deleted after use.\n 
    - If none, only the run-specific copies of the master executable and weather files will be deleted after use.
    """, required=True, type = str, choices = {"full","some","none"})
    parser.add_argument("-g", "--generations", help = "Number of generations to run SMC for.", required=True, type = int)
    parser.add_argument("-n", "--popsize", help = "Number of particles per generation.", required=True, type = int)
    parser.add_argument("-s", "--SWMMparameters", help = "Number of SWMM parameters to use.", required = False, type = int)
    parser.add_argument("-v", "--VVWMparameters", help = "Number of VVWM parameters to use.", required = False, type = int)
    parser.add_argument("-p", "--parallel", help = """Should a dask distributed (parallel) sampler 
    be used instead of a single core (non-parallel) sampler?""", action = 'store_true')#, required=False, type = bool)
    args = parser.parse_args()
    # for simulation part
    # make them lowercase and give them shorter names to go by
    mode, swmm_cleanup, vvwm_cleanup = args.mode.lower(), args.cleanup[0].lower(), args.cleanup[1].lower()
    # make sure user provided all legal values
    assert mode in ["debug", "test", "run"], 'Acceptable values of --mode include \'debug\', \'test\', and \'run\', not \'' + mode + '\'.'
    assert swmm_cleanup in ["full", "some", "none"], 'Acceptable values of --swmmcleanup include \'full\', \'some\', and \'none\', not \'' + swmm_cleanup + '\'.'
    assert vvwm_cleanup in ["full", "some", "none"], 'Acceptable values of --vvwmcleanup include \'full\', \'some\', and \'none\', not \'' + vvwm_cleanup + '\'.'
    # for smc part
    ngen, npop = min(max(1,args.generations), 15), min(max(args.popsize,2), 1000)
    # optionals
    debug_params = []
    if args.SWMMparameters:
        assert 1 <= args.SWMMparameters <= 16, 'SWMMparameters must be an integer in [1,16]'
        debug_params = [i for i in range(args.SWMMparameters)]
    if args.VVWMparameters:
        assert 1 <= args.VVWMparameters <= 18, 'VVWMparameters must be an integer in [1,18]'
        debug_params = debug_params + [i+16 for i in range(args.VVWMparameters)]
    # make the model using the user-input mode and cleanup levels
    model = make_model(mode = mode, swmm_cleanup = swmm_cleanup, vvwm_cleanup = vvwm_cleanup)#, debug_params = debug_params)
    #
    # pyabc_construct stage starts here!
    #
    # Priors. Get the values from that csv. SWMM at first. Then VVWM. Then link them together.
    swmm_ranges = pd.read_csv(os.path.join(master_path, "swmm_param_priors.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])
    vvwm_ranges = pd.read_csv(os.path.join(master_path, "vvwm_param_priors.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])
    param_ranges = pd.concat([swmm_ranges, vvwm_ranges], axis = 0)
    # take just a tiny subset if in debug mode
    if mode == "debug" and debug_params:
        param_ranges = param_ranges.iloc[debug_params]
    # make the dataframe into a dictionary
    priors = param_ranges.to_dict("index")
    # make the dictionary into a pyabc distribution object
    prior = pyabc.Distribution(**{key: pyabc.RV("uniform", loc = v['Min'], scale = v['Range']) for key, v in priors.items()})
    # import the dictionary of observed data (benchmarks for comparison) depending on the mode
    if mode == 'debug':
        with open(os.path.join(main_path, 'master_debug','debug_obs_data.txt'),'r') as read_file:
            obs_dict = eval(read_file.read())
    elif mode == 'test':
        with open(os.path.join(main_path, 'master_test','test_obs_data.txt'),'r') as read_file:
            obs_dict = eval(read_file.read())
    elif mode == 'run':
        with open(os.path.join(main_path, 'master','obs_data.txt'),'r') as read_file:
            obs_dict = eval(read_file.read())
    # use the single core sampler for now because dask is being fussy
    # Parallel
    if args.parallel:
        from dask.distributed import Client
        client = Client()
        sampler = pyabc.sampler.DaskDistributedSampler(dask_client = client)
    else:
        sampler = pyabc.sampler.SingleCoreSampler()
    # make the process more transparent
    sampler.sample_factory.record_rejected = True
    sampler.show_progress = True
    # Set up a sqlite db directory with a unique random identifier
    dbid = uuid.uuid4().hex[0:8]
    print("Database ID: " + dbid)
    # compose path to run-specific database
    database_dir = os.path.join(temp_path, 'results_db')  
    if not os.path.exists(database_dir):
        os.mkdir(database_dir)
    db_path = ("sqlite:///" + os.path.join(database_dir, "test_pyabc_" + dbid + ".db"))
    # make a file to hold onto these NSEs for our own record
    with open(os.path.join(temp_path, "NSEs_" + dbid + ".txt"), "w") as nse_file:
        nse_file.write("NSEs\n")
    # define NSE Distance function: Calculate NSE with hydroeval library and the subtract 1 from it to get NSED
    def nse(x, x_0):
        nse = he.evaluator(he.nse, simulation_s = np.array(list(x.values())), evaluation = np.array(list(x_0.values())))[0]
        print("nse ", nse)
        # make record
        with open(os.path.join(temp_path, "NSEs_" + dbid + ".txt"),"a") as nse_file:
            nse_file.write(str(nse)+"\n")
        return nse
    # NSEs live on [-infinity, 1] and the best NSE is 1
    # Distances live on [0, infinity and the best distance is 0
    # Therefore, let's make a formula that's informed by NSE, but has range and properties of a distance function
    NSED = pyabc.SimpleFunctionDistance(fun = lambda x, x_0: 1 - nse(x, x_0))
    # make the pyabc Seq Monte Carlo object
    abc = pyabc.ABCSMC(model, prior, population_size = pyabc.ConstantPopulationSize(npop), sampler = sampler, distance_function = NSED)
    # initialize a new run
    abc.new(db_path, obs_dict)
    # run it!
    history = abc.run(max_nr_populations=ngen, minimum_epsilon=0.2)