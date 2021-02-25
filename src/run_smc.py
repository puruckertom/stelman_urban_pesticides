import sys, argparse
from simulation_procedure import make_model
from tools.paths import *
import pandas as pd
import hydroeval as he
import numpy as np
import pyabc
import uuid
"""
import os, pandas as pd, numpy as np, pickle
import pytest_shutil, shutil, regex as re, dask, uuid
from pyswmm import Simulation
from pyswmm.swmm5 import SWMMException
import swmmtoolbox.swmmtoolbox as swmmtoolbox
import time
from pyswmm.lib import DLL_SELECTION
import subprocess

from tools.paths import *

from tools.functions import *

dll_path = DLL_SELECTION()
dll_bn = os.path.basename(dll_path)

# sub_list_area has the areas of each subcatchment in order
# used by 05
with open(os.path.join(master_path,'sub_list_area.txt'),'r') as read_file:
    sub_list_area = eval(read_file.read())

# sub_ids is the dictionary detailing which subcatchments make up which outfalls
# used by 05
with open(os.path.join(master_path,'outfall_partition.txt'),'r') as read_file:
    sub_ids = eval(read_file.read())

# from simulation_procedure import model
# from tools.paths import *
# import pandas as pd, pyabc, hydroeval as he, numpy as np

# # for debug mode 
# from pyDOE2 import lhs
# from scipy.stats import uniform
# swmm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges.csv"), index_col=0,
#                         usecols = ["Parameter","Min", "Range"])
# vvwm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges_vvwm.csv"), index_col=0,
#                         usecols = ["Parameter","Min", "Range"])
# param_ranges = pd.concat([swmm_ranges, vvwm_ranges], axis = 0)
# del swmm_ranges, vvwm_ranges

# Import Observed Data
obs_data = pd.read_csv(obs_path, usecols=["Sample_date", "Site_code"],
                       parse_dates=["Sample_date"])


def make_model(mode, swmm_cleanup, vvwm_cleanup):
    # activate test mode
    inp_path = set_inp_path(mode)
    print(inp_path)

    if mode == "debug":
        outfalls = ['outfall_31_28']
    else:
        outfalls = ['outfall_31_26', 'outfall_31_28', 'outfall_31_29', 'outfall_31_35',
                    'outfall_31_36', 'outfall_31_38', 'outfall_31_42']


    def model(params1):
    
        # import test parameter input
        if mode == "debug":

            # get them into the objects we need
            swmm_keys = list(params1.keys())[:1]
            vvwm_keys = list(params1.keys())[1:]
            swmm_params = {key: params1[key] for key in swmm_keys}
            vvwm_params = {key: params1[key] for key in vvwm_keys}
            # lhs1 = lhs(n=36, samples=1)[0]
            # for i in range(0,36):
            #     lhs1[i] = param_ranges["Min"][i] + (lhs1[i])*(param_ranges["Range"][i])
            # params = {}
            # for key,value in list(zip(param_ranges.index, lhs1)):
            #     params[key] = value
            # # now the arguments being passed in
            # for key in list(params1.keys()):
            #     params[key] = params1[key]
            # params1 = params
            # # error prevention
            # if params1['MaxRate'] < params1['MinRate']:
            #     params1['MaxRate'], params1['MinRate'] = params1['MinRate'], params1['MaxRate']
            # if params1['FC'] < params1['WP']:
            #     params1['FC'], params1['WP'] = params1['WP'], params1['FC']
                
            # # get them into the objects we need
            # swmm_keys = list(params1.keys())[:17]

            # # wp and fc in this object should be switched. FOR NOW:
            # swmm_keys = swmm_keys[:10] + ['WP', 'FC'] + swmm_keys[12:]

            # vvwm_keys = list(params1.keys())[17:]
            # swmm_params = {key: params1[key] for key in swmm_keys}
            # vvwm_params = {key: params1[key] for key in vvwm_keys}
            
        if mode == "test":

            
            # error prevention
            if params1['MaxRate'] < params1['MinRate']:
                params1['MaxRate'], params1['MinRate'] = params1['MinRate'], params1['MaxRate']
            if params1['FC'] < params1['WP']:
                params1['FC'], params1['WP'] = params1['WP'], params1['FC']
                
            # get them into the objects we need
            swmm_keys = list(params1.keys())[:17]

            # wp and fc in this object should be switched. FOR NOW:
            swmm_keys = swmm_keys[:10] + ['WP', 'FC'] + swmm_keys[12:]

            vvwm_keys = list(params1.keys())[17:]
            swmm_params = {key: params1[key] for key in swmm_keys}
            vvwm_params = {key: params1[key] for key in vvwm_keys}

        # spin up simulation id with this cool too for generating random ids
        sid = uuid.uuid4().hex[0:8]
        # set up logger
        loginfo, logerror = log_prefixer(sid)

        # ## Step 1: make simulation-specific SWMM items

        # make paths to (soon-to-be) directory & files for this simulation using sid
        # directory
        sdir_path = os.path.join(temp_path, sid)
        # input file
        sinp_path = os.path.join(sdir_path, sid) + ".inp"
        # output file
        sout_path = os.path.join(sdir_path, sid) + ".out"
        # report file
        srpt_path = os.path.join(sdir_path, sid) + ".rpt"
        # dynamic link library file
        ## basename
        sdll_bn = dll_bn[:dll_bn.rindex(".")] + '-' + sid + dll_bn[dll_bn.rindex("."):]
        ## full path
        sdll_path = os.path.join(sdir_path, sdll_bn) 

        # this ifelse stuff is needless inside this fxn, but it's harmless, so I'll leave it for now
        # make the directory
        if not os.path.exists(sdir_path):
            os.mkdir(sdir_path)
            print("Folder ", sid, " created", "\n")
        else:
            print("Folder ", sid, "already exists")

        # write input file
        with open(inp_path, "r") as read_file, open(sinp_path, "w") as write_file:
            filelines = read_file.readlines()

            # first we need to correct some absolute paths, because they are currently only set to work on the author's computer
            filelines = replace_infile_abspaths(filelines = filelines)
            
            # if mode == "debug" or mode == "test":
            if mode == 'debug':
                filelines[172:(172 + 113)] = editted_lines(swmm_dict = swmm_params, Num = 113, row_0 = 172, parameter = "NImperv", Col = 1, flines = filelines)
            
            #elif mode == "test":
            else:
                # 113 = number of subcatchments
                #for c, par in enumerate(["NImperv", "NPerv", "SImperv", "SPerv", "PctZero"]):
                for c, par in enumerate(swmm_keys[0:5]):
                    filelines[172:(172 + 113)] = editted_lines(swmm_dict = swmm_params, Num = 113, row_0 = 172, parameter = par, Col = c+1, flines = filelines)
                #for c, par in enumerate(["MaxRate", "MinRate", "Decay", "DryTime"]):
                for c, par in enumerate(swmm_keys[5:9]):
                    filelines[289:(289 + 113)] = editted_lines(swmm_dict = swmm_params, Num = 113, row_0 = 289, parameter = par, Col = c+1, flines = filelines)

                # 1 = number of aquifers
                #for c, par in enumerate(["Por", "WP", "FC", "Ksat"]):
                for c, par in enumerate(swmm_keys[9:13]):
                    filelines[406:(406 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 406, parameter = par, Col = c+1, flines = filelines)

                # 195 = number of conduits
                # filelines[734:(734 + 195)] = editted_lines(swmm_dict = swmm_params, Num = 195, row_0 = 734, parameter = "Rough", Col = 4, flines = filelines)

                # 1 = number of pollutants
                filelines[1125:(1125 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1125, parameter = "Kdecay", Col = 5, flines = filelines)
                filelines[1371:(1371 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1371, parameter = "BCoeff2", Col = 4, flines = filelines)
                filelines[1377:(1377 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1377, parameter = "WCoeff2", Col = 4, flines = filelines)
            
            write_file.writelines(filelines)

        # copy a dll file into sdll_path
        shutil.copyfile(dll_path, sdll_path)

        # #### Execute simulation
        # In test mode, should take about 2-3 minutes. 
        # In run mode, should take a while.

        # This is not needed in practice
        # # delete pre-existing .out, if present, in order to run swmm agreeably
        # if os.path.exists(sout_path):
        #     loginfo("Deleting current copy of <" + sout_path + "> so new copy can be created.")
        #     #print("Deleting current copy of <NPlesantCreek.out> so new copy can be created.")
        #     os.remove(sout_path)

        # load the model {no interaction, write (binary) results to sout_path, use the specified dll}

        sim = Simulation(inputfile=sinp_path, reportfile=srpt_path, outputfile=sout_path, swmm_lib_path=sdll_path)
        # if no errors were thrown, we procede with the simulation

        # simulate the loaded model
        loginfo("Executing SWMM simmulation with no interaction.")# Input from <" + sinp_path + ">. Will store output in <" + sout_path + ">.")
        with sim as s:
            for step in s:
                pass

        # extract swmm outputs with swmmtoolbox and delete expensive binary files
        lab1, lab2 = 'subcatchment,,Runoff_rate', 'subcatchment,,Bifenthrin'
        runf = swmmtoolbox.extract(sout_path, lab1)
        bif = swmmtoolbox.extract(sout_path, lab2)

        if swmm_cleanup == 'full':
            loginfo("Deleting swmm temp files to free up memory.")#<" + sdir_path + "> contents to free up memory.")
            os.system("rm " + sdir_path + "/*")
        elif swmm_cleanup == 'some':
            loginfo("Deleting large swmm temp files to free up memory.")#<" + sdir_path + "> contents to free up memory.")
            os.system("rm " + sout_path + " " + sdll_path)
        elif swmm_cleanup == 'none':
            loginfo("Deleting swmm dll file to free up memory.")#<" + sdir_path + "> contents to free up memory.")
            os.system("rm " + sdll_path)
    
        # compute daily averages
        runf = runf.resample('D').mean()
        bif = bif.resample('D').mean()

        # conversion for vvwm for runoff and bifenthrin
        runf = runf.mul(86400).mul(0.01).div(sub_list_area)
        bif = bif.mul(runf.values)

        for o in outfalls:
            # set pathways
            outfall_dir = os.path.join(sdir_path, o)
            
            # this ifelse stuff is needless inside this fxn, but it's harmless, so I'll leave it for now
            # make the directory
            if not os.path.exists(outfall_dir):
                os.mkdir(outfall_dir)
                print("Folder ", sid, o, " created", "\n")
            else:
                print("Folder ", sid, o, "already exists")
            
            # create .zts file 
            # has 1 at each subcatchment of outfall[o], and 0 at the rest
            weights = np.array([(1 if x in sub_ids[o] else 0) for x in range(113)])
            # we want the daily totals the runoff and bifenthrin within each outfall
            # use dot product and the weights vector to evaluate this on the df
            runf_sum = runf @ weights
            bif_sum = bif @ weights
            
            # we want to combine tables and add filler and date-part columns
            vvwm_df = pd.DataFrame({"year": runf_sum.index.year,
                                "month": runf_sum.index.month,
                                "day": runf_sum.index.day,
                                "runf_sum": runf_sum,
                                "B": 0,
                                "bif_sum": bif_sum,
                                "MEp":0
                            })
            
            # read out into comma-delimited .txt file
            vvwm_df.to_csv(os.path.join(outfall_dir, "output.zts"), 
                        header=False, index=False, sep=',')

            # for this to work, we need to write 3 blank lines to the beginning
            # read in the file we just wrote
            with open(os.path.join(outfall_dir, "output.zts"), "r") as read_file:
                filelines = read_file.readlines()
            # write this back into it with 3 lines added to the beginning
            with open(os.path.join(outfall_dir, "output.zts"), "w") as write_file:
                # write blanks to dummy file
                write_file.write('\n\n\n')
                # read lines from original and append to dummy file
                write_file.writelines(filelines) 

        # create a blank df to populate
        output_df = pd.DataFrame()

        for o in outfalls:
            # set pathways
            outfall_dir = os.path.join(sdir_path, o)
            outfall_file = os.path.join(outfall_dir, "vvwmTransfer.txt")
            
            with open(vvwmTransfer_path,"r") as read_file:
                filelines = read_file.readlines()
            
            # if mode == "debug" or mode == "test":
            if mode == 'debug':
                filelines[4] = str(vvwm_params[vvwm_keys[0]]) + "\n"
                
            #elif mode == "test":
            else:

                # why is it 1:8 and not 0:8? Is it a typo?
                for c, param in enumerate(list(vvwm_keys)[0:8]): # changed from 1:8 2/19/21
                    filelines[c+4] = str(vvwm_params[param]) + "\n"

                filelines[17] = str(vvwm_params[vvwm_keys[8]]) + "\n"

                for c, param in enumerate(list(vvwm_keys)[9:15]):
                    filelines[c+40] = str(vvwm_params[param]) + "\n"

                for c, param in enumerate(list(vvwm_keys)[15:19]):
                    filelines[c+47] = str(vvwm_params[param]) + "\n"

            # update pathways
            filelines[0] = os.path.join(outfall_dir, "output") + '\n'
            filelines[29] = os.path.join(outfall_dir, "vvwm_wet.dvf") + '\n'
            filelines[68] = os.path.join(outfall_dir, "output_NPlesant_Custom_parent_daily.csv") + '\n'
            filelines[69] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg1_daily.csv") + '\n'
            filelines[70] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg2_daily.csv") + '\n'
            filelines[71] = os.path.join(outfall_dir, "output_NPlesant_Custom_parent_analysis.txt") + '\n'
            filelines[72] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg1_analysis.txt") + '\n'
            filelines[73] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg2_analysis.txt") + '\n'
            filelines[74] = os.path.join(outfall_dir, "output_NPlesant_Custom_parent_deem.rdf") + '\n'
            filelines[75] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg1_deem.rdf") + '\n'
            filelines[76] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg2_deem.rdf") + '\n'
            filelines[77] = os.path.join(outfall_dir, "output_NPlesant_Custom_parent_calendex.rdf") + '\n'
            filelines[78] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg1_calendex.rdf") + '\n'
            filelines[79] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg2_calendex.rdf") + '\n'
            filelines[80] = os.path.join(outfall_dir, "output_NPlesant_Custom_parent_esa.txt") + '\n'
            filelines[81] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg1_esa.txt") + '\n'
            filelines[82] = os.path.join(outfall_dir, "output_NPlesant_Custom_deg2_esa.txt") + '\n'
                
            with open(outfall_file, "w") as write_file:
                # write out file
                write_file.writelines(filelines)
            
            # copy weather file into new file location
            if mode == 'debug':
                old_wet_path = os.path.join(main_path, "master_debug", "vvwm_wet.dvf")
            elif mode == 'test':
                old_wet_path = os.path.join(main_path, "master_test", "vvwm_wet.dvf")
            elif mode =='run':
                old_wet_path = os.path.join(weather_path, "vvwm_wet.dvf")
            new_wet_path = os.path.join(outfall_dir, "vvwm_wet.dvf")
            shutil.copyfile(old_wet_path, new_wet_path)

            # copy exe into new file location
            if sys.platform == "linux" or sys.platform == "linux2":
                exe_bn = "vvwm"
            if sys.platform == "win32":
                exe_bn = "VVWM.exe"
            old_exe_path = os.path.join(exe_path, exe_bn)
            new_exe_path = os.path.join(outfall_dir, exe_bn)
            shutil.copyfile(old_exe_path, new_exe_path)
            
            # run vvwm.exe (vvwm.exe [...]/outfall_31_xx/vvwmTransfer.txt)
            command = new_exe_path + ' ' + outfall_file
            subprocess.call(command)
            
            # read in produced data from the output of the vvwm run we just completed
            output = pd.read_csv(filelines[68][:-1], 
                                usecols = [1], skiprows=5, names = ["davg_bif_conc"])*1000000
            if mode == 'debug':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=103, freq='D')
            elif mode == 'test':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=778, freq='D')
            elif mode == 'run':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=3287, freq='D')
            output['Site_code'] = o[-5:]
            output = output.merge(obs_data, how = "inner", on = ['Sample_date','Site_code'])
            output.set_index([np.datetime_as_string(output.Sample_date, unit = 'D'),'Site_code'], inplace = True)
            output_df = output_df.append(output[['davg_bif_conc',]], ignore_index = False)

        # sort by date and site
        output_df = output_df.set_index([["_".join([a,b[3:]]) for a,b in output_df.index]]).sort_index()
        output_dict = output_df.to_dict()['davg_bif_conc']

        # vvwm cleanup
        if vvwm_cleanup == 'none' or vvwm_cleanup == 'some' or vvwm_cleanup == 'full':
            exe_ = os.path.join(sdir_path, "outfall_31_??", exe_bn)
            wet_ = os.path.join(sdir_path, "outfall_31_??", "vvwm_wet.dvf")
            loginfo("Deleting vvwm exe and weather files to free up memory.")
            os.system("rm " + exe_ + " " + wet_)
        if vvwm_cleanup == 'some' or vvwm_cleanup == 'full':
            zts_ = os.path.join(sdir_path, "outfall_31_??", "output.zts")
            analysis_ = os.path.join(sdir_path, "outfall_31_??", "output_NPlesant_Custom_parent_analysis.txt")
            transfer_ = os.path.join(sdir_path, "outfall_31_??", "vvwmTransfer.txt")
            loginfo("Deleting internal vvwm files to free up memory.")
            os.system("rm " + zts_ + " " + analysis_ + " " + transfer_)
        if vvwm_cleanup == 'full':
            daily_ = os.path.join(sdir_path, "outfall_31_??", "output_NPlesant_Custom_parent_daily.csv")
            loginfo("Deleting vvwm results file to free up memory.")
            os.system("rm" + daily_)
        if vvwm_cleanup == 'full' and swmm_cleanup == 'full':
            loginfo("Deleting temp folder.")
            os.system("rm -r " + sdir_path + "/")
        
        # os.system("rm -r " + sdir_path + "/")
        return(output_dict)

    return model
"""

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
    parser.add_argument("-s", "--swmmcleanup", help = """
    Specify SWMM cleanup level as 'full', 'some', or 'none':\n
    If full, each swmm-related file will be deleted after fulfilling its use to save memory.\n
    If some, executable and binary swmm-related files will be deleted after use, but human-readable files won't.\n
    If none, only the run-specific copy of the master executable (dll) file will be deleted after use.
    """, required=True, type = str)
    parser.add_argument("-v", "--vvwmcleanup", help = """
    Specify vvwm cleanup level as 'full', 'some', or 'none':\n
    If full, each vvwm-related file will be deleted after fulfilling its use to save memory.\n
    If some, all vvwm-related files but one, the human-readable daily output csv file, will be deleted after use.\n
    If none, only the run-specific copies of the master executable and weather files will be deleted after use.
    """, required=True, type = str)
    args = parser.parse_args()
    # make them lowercase and give them shorter names to go by
    mode, swmm_cleanup, vvwm_cleanup = args.mode.lower(), args.swmmcleanup.lower(), args.vvwmcleanup.lower()
    # make sure user provided all legal values
    assert mode in ["debug", "test", "run"], 'Acceptable values of --mode include \'debug\', \'test\', and \'run\', not \'' + mode + '\'.'
    assert swmm_cleanup in ["full", "some", "none"], 'Acceptable values of --swmmcleanup include \'full\', \'some\', and \'none\', not \'' + swmm_cleanup + '\'.'
    assert vvwm_cleanup in ["full", "some", "none"], 'Acceptable values of --vvwmcleanup include \'full\', \'some\', and \'none\', not \'' + vvwm_cleanup + '\'.'
    '''if mode not in ["debug", "test", "run"]:
        raise ValueError('Acceptable values of --mode include \'debug\', \'test\', and \'run\', not \'' + mode + '\'.')
    if swmm_cleanup not in ["full", "some", "none"]:
        raise ValueError('Acceptable values of --swmmcleanup include \'full\', \'some\', and \'none\', not \'' + swmm_cleanup + '\'.')
    if vvwm_cleanup not in ["full", "some", "none"]:
        raise ValueError('Acceptable values of --vvwmcleanup include \'full\', \'some\', and \'none\', not \'' + vvwm_cleanup + '\'.')'''
    # make the model using the user-input mode and cleanup levels
    model = make_model(mode = mode, swmm_cleanup = swmm_cleanup, vvwm_cleanup = vvwm_cleanup)
    #
    # pyabc_construct stage starts here!
    #
    # Priors. Get the values from that csv. SWMM at first. Then VVWM. Then link them together.
    swmm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])
    vvwm_ranges = pd.read_csv(os.path.join(master_path, "lhs_param_ranges_vvwm.csv"), index_col=0, usecols = ["Parameter","Min", "Range"])
    param_ranges = pd.concat([swmm_ranges, vvwm_ranges], axis = 0)
    # take just a tiny subset if in debug mode
    if mode == "debug":
        param_ranges = param_ranges.loc[['NImperv','kd']]
    # make the dataframe into a dictionary
    priors = param_ranges.to_dict("index")
    # make the dictionary into a pyabc distribution object
    # borrowed from Jeff: <https://github.com/JeffreyMinucci/bee_neonic_abc/blob/master/pyabc_run.ipynb>
    prior = pyabc.Distribution(**{key: pyabc.RV("uniform", loc = v['Min'], scale = v['Range']) for key, v in priors.items()})
    # import the dictionary of observed data (benchmarks for summary statistic comparison) depending on the mode
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
    sampler = pyabc.sampler.SingleCoreSampler()
    # make the process more transparent
    sampler.sample_factory.record_rejected = True
    sampler.show_progress = True
    # Set up a sqlite db directory with a unique random identifier
    dbid = uuid.uuid4().hex[0:8]
    print("Database ID: " + dbid)
    database_dir = os.path.join(temp_path, 'results_db')  
    if not os.path.exists(database_dir):
        os.mkdir(database_dir)
    db_path = ("sqlite:///" + os.path.join(database_dir, "test_pyabc_" + dbid + ".db"))
    # make a file to hold onto these NSEs for our own record
    with open(os.path.join(temp_path, "NSEs_" + dbid + ".txt"), "w") as nse_file:
        nse_file.write("NSEs\n")
    # make NSE and NSE distance calculating functions
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
    abc = pyabc.ABCSMC(model, prior, population_size = pyabc.ConstantPopulationSize(40), sampler = sampler, distance_function = NSED)
    # initialize a new run
    abc.new(db_path, obs_dict)
    # run it!
    history = abc.run(max_nr_populations=2, minimum_epsilon=0.2)