#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import required modules
import os, pandas as pd, numpy as np
import pytest_shutil, shutil, regex as re, uuid
from pyswmm import Simulation
import swmmtoolbox.swmmtoolbox as swmmtoolbox
from pyswmm.lib import DLL_SELECTION
import subprocess
# import objects from tools scripts
from tools.paths import *
from tools.functions import *

# save the path to the default dynamic link library used for running SWMM
dll_path = DLL_SELECTION()
dll_bn = os.path.basename(dll_path)


# ## Import objects from files

# In[2]:


# Import a list of the areas of each subcatchment in order
with open(os.path.join(master_path,'sub_list_area.txt'),'r') as read_file:
    sub_list_area = eval(read_file.read())

# Import the dictionary detailing which subcatchments make up which outfalls
with open(os.path.join(master_path,'outfall_partition.txt'),'r') as read_file:
    sub_ids = eval(read_file.read())

# Import (a dictionary of) a test set of inputs for the params1 argument of the model() function
with open(os.path.join(main_path,'master_test','test_params.txt'),'r') as read_file:
    params = eval(read_file.read())

# Import observed data
obs_data = pd.read_csv(obs_path, usecols=["Sample_date", "Site_code"],
                       parse_dates=["Sample_date"])

# ## Define make_model() function:

# In[7]:


'''
Creates model() function and sets up environment according to specifics provided
 Inputs: mode <str> -Either 'debug', 'test', or 'run'-
     -Debug mode tests code quickly to expose bugs in runs-
     -Test mode impersonates a shortened, but complete, run of the simulation procedure, using a small sample of the data-
     -Run mode (when the code is fully baked) takes forever and gives results based on the whole dataset-
   swmm_cleanup <str> -Either 'full', 'some', or 'none'-
     -If full, each swmm-related file will be deleted after fulfilling its use to save memory-
     -If some, executable and binary swmm-related files will be deleted after use, but human-readable files won't-
     -If none, only the run-specific copy of the master executable (dll) file will be deleted after use-
   vvwm_cleanup <str> -Either 'full', 'some', or 'none'-
     -If full, each vvwm-related file will be deleted after fulfilling its use to save memory-
     -If some, all vvwm-related files but one, the human-readable daily output csv file, will be deleted after use-
     -If none, only the run-specific copies of the master executable and weather files will be deleted after use-
   debug_params <list of unique integers in [0,33]> -indices of parameters to use in model() instead of full set-
     -If empty (default), signals model() to use full set-
 Output: model <function> -Takes in dict of parameters, runs a simulation, and outputs dict of summary statistics-
'''
def make_model(mode, swmm_cleanup, vvwm_cleanup, debug_params = []):


    # ### Activate mode

    # In[8]:


    # get outfall names into a list: 1 outfall name if debug mode, 7 outfall names otherwise
    if mode == "debug":
        outfalls = ['outfall_31_28']
    else:
        outfalls = ['outfall_31_26', 'outfall_31_28', 'outfall_31_29', 'outfall_31_35',
                    'outfall_31_36', 'outfall_31_38', 'outfall_31_42']


    # In[9]:


    # set input path according to mode given
    inp_path = set_inp_path(mode)
    print(inp_path)


    # ## Define model() function:

    # In[12]:


    '''
    Simulates bifenthin application events and rain events over the course of time in an urban setting.
    Uses multistage model to predict periodic subcatchment contamination fate based on simulated events and given parameters.
    Returns bifenthrin concentration estimates specific to a set of date-site pairs for which field data exists.
    Inputs: params1 <dict> -Dictionary of parameters (of interest) used in SWMM and VVWM modeling procedures and their values-
    Output: output_dict <dict> -Dictionary summarizing model outcomes-
        -Keys: a date followed by a site identifier (ie. yyyy-mm-dd_xx, where xx is the 2-digit suffix of the outfall name)-
        -Values: a (date, site)-specific average of the simulation-estimated bifenthrin concentrations its key pertains to-
    '''
    def model(params1):


        # #### Set up input variables and simulation id

        # In[13]:


        if mode == "debug":
            # Fill in params1 with whichever of the 34 parameters are not provided as inputs with the defaults from test_params
            for key, value in params1.items():
                params[key] = value
            params1 = params

        # Error prevention: 
        # In real life, these situations would be impossible.
        # Since these parameters are simulated, they need to be watched out for and fixed.
        # If they are overlooked, the model will fail early during simulation.
        if params1['MaxRate'] < params1['MinRate']:
            params1['MaxRate'], params1['MinRate'] = params1['MinRate'], params1['MaxRate']
        if params1['FC'] < params1['WP']:
            params1['FC'], params1['WP'] = params1['WP'], params1['FC']

        # partition parameters into the ones for swmm and the ones for vvwm
        swmm_keys = list(params1.keys())[:16] # 3/2
        vvwm_keys = list(params1.keys())[16:] # 3/2
        swmm_params = {key: params1[key] for key in swmm_keys}
        vvwm_params = {key: params1[key] for key in vvwm_keys}


        # In[14]:


        # spin up simulation id with this cool too for generating random ids
        sid = uuid.uuid4().hex[0:8]
        # set up logging
        loginfo, logerror = log_prefixer(sid)


        # ## Step 1: make simulation-specific SWMM items

        # In[15]:


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
        ## dll basename
        sdll_bn = dll_bn[:dll_bn.rindex(".")] + '-' + sid + dll_bn[dll_bn.rindex("."):]
        ## dll full path
        sdll_path = os.path.join(sdir_path, sdll_bn) 


        # In[16]:


        # make the directory for which we just constructed a path
        if not os.path.exists(sdir_path):
            os.mkdir(sdir_path)
            print("Folder ", sid, " created", "\n")
        else:
            print("Folder ", sid, "already exists")

        # make swmm input file:
        with open(inp_path, "r") as read_file, open(sinp_path, "w") as write_file:
            filelines = read_file.readlines()

            # first we need to correct some absolute paths, because they are currently only set to work on the author's computer
            filelines = replace_infile_abspaths(filelines = filelines)
            
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

            # 1 = number of pollutants
            filelines[1125:(1125 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1125, parameter = "Kdecay", Col = 5, flines = filelines)
            filelines[1371:(1371 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1371, parameter = "BCoeff2", Col = 4, flines = filelines)
            filelines[1377:(1377 + 1)] = editted_lines(swmm_dict = swmm_params, Num = 1, row_0 = 1377, parameter = "WCoeff2", Col = 4, flines = filelines)
            
            # write the changes to the file
            write_file.writelines(filelines)

        # copy a dll file into sdll_path
        shutil.copyfile(dll_path, sdll_path)


        # #### Execute simulation
        # In test mode, should take about 2-3 minutes. 
        # In run mode, should take a while.

        # In[17]:


        # Error prevention
        # delete pre-existing .out, if present, in order to run swmm agreeably
        if os.path.exists(sout_path):
            loginfo("Deleting current copy of <" + sout_path + "> so new copy can be created.")
            #print("Deleting current copy of <NPlesantCreek.out> so new copy can be created.")
            os.remove(sout_path)

        # load the model {no interaction, write (binary) results to sout_path, use the specified dll}
        sim = Simulation(inputfile=sinp_path, reportfile=srpt_path, outputfile=sout_path, swmm_lib_path=sdll_path)
        # simulate the loaded model
        loginfo("Executing SWMM simmulation with no interaction. Input from <" + sinp_path + ">. Will store output in <" + sout_path + ">.")
        # use for-loop to avoid runtime costs due to excessive logging
        with sim as s:
            for step in s:
                pass


        # #### Get the info to a safe place and then clean up

        # In[18]:


        # extract swmm outputs with swmmtoolbox and delete expensive binary files
        lab1, lab2 = 'subcatchment,,Runoff_rate', 'subcatchment,,Bifenthrin'
        runf = swmmtoolbox.extract(sout_path, lab1)
        bif = swmmtoolbox.extract(sout_path, lab2)


        # In[19]:


        # clean up
        if swmm_cleanup == 'full':
            loginfo("Deleting swmm temp files to free up memory.")
            rm(os.path.join(sdir_path,"*"))
        elif swmm_cleanup == 'some':
            loginfo("Deleting large swmm temp files to free up memory.")
            rm(sout_path, sdll_path)
        elif swmm_cleanup == 'none':
            loginfo("Deleting swmm dll file to free up memory.")
            rm(sdll_path)


        # In[20]:


        # compute daily averages
        runf = runf.resample('D').mean()
        bif = bif.resample('D').mean()


        # In[21]:


        # unit conversion for vvwm: runoff and bifenthrin
        ## Runoff
        ### 1. multiply by 86400 to convert days into seconds
        ### 2. multiply by 0.01 to convert square meters to hectares
        ### 3. divide by area of subcatchment (hectares) to get per-hectare volumes
        runf = runf.mul(86400).mul(0.01).div(sub_list_area)
        ## Bifenthrin
        ### multiply by runoff volume to get ????
        bif = bif.mul(runf.values)


        # In[22]:


        for o in outfalls:
            # create outfall directory path
            outfall_dir = os.path.join(sdir_path, o)
            
            # make the directory for which we just constructed a path
            if not os.path.exists(outfall_dir):
                os.mkdir(outfall_dir)
                print("Folder ", sid, o, " created", "\n")
            else:
                print("Folder ", sid, o, "already exists")
            
            # create a vector to strain the outfall[o]-specific subcatchments out of the 113 subcatchments
            # weights vector has 1 at each index corresponding to a subcatchment of outfall[o], and 0 at the rest
            weights = np.array([(1 if x in sub_ids[o] else 0) for x in range(113)])
            # we want the daily totals of the runoff and bifenthrin within outfall[o]
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
            
            # make swmm output data into vvwm input file
            outfall_path = os.path.join(outfall_dir, "output.zts")
            # read out into comma-delimited .zts file
            vvwm_df.to_csv(outfall_path, header=False, index=False, sep=',')

            # for this to work, we need to write 3 blank lines to the beginning
            with open(outfall_path, "r") as read_file:
                filelines = read_file.readlines()
            # write this back into it with 3 lines added to the beginning
            with open(outfall_path, "w") as write_file:
                # write blanks to file
                write_file.write('\n\n\n')
                # append original lines to file after blank lines
                write_file.writelines(filelines)


        # In[23]:


        # create a blank df to populate
        output_df = pd.DataFrame()

        for o in outfalls:
            # create outfall directory path and outfall file path
            outfall_dir = os.path.join(sdir_path, o)
            outfall_file = os.path.join(outfall_dir, "vvwmTransfer.txt")
            
            # make vvwm setup file for outfall o:
            with open(vvwmTransfer_path,"r") as read_file:
                filelines = read_file.readlines()
            
            # update parameter values
            for c, param in enumerate(list(vvwm_keys)[0:6]): 
                filelines[c+4] = str(vvwm_params[param]) + "\n"
            filelines[11] = str(vvwm_params[vvwm_keys[6]]) + "\n"
            filelines[17] = str(vvwm_params[vvwm_keys[7]]) + "\n"
            for c, param in enumerate(list(vvwm_keys)[8:14]):
                filelines[c+40] = str(vvwm_params[param]) + "\n"
            for c, param in enumerate(list(vvwm_keys)[14:18]):
                filelines[c+47] = str(vvwm_params[param]) + "\n"

            # Update path to (swmm) output (aka, vvwm input) data file
            filelines[0] = os.path.join(outfall_dir, "output") + '\n'
            # Update path to weather file
            filelines[29] = os.path.join(outfall_dir, "vvwm_wet.dvf") + '\n'
            # Insert paths that don't and won't exist, but need to be included as a technicality or else the model won't run
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
            elif sys.platform == "win32":
                exe_bn = "VVWM.exe"
            old_exe_path = os.path.join(exe_path, exe_bn)
            new_exe_path = os.path.join(outfall_dir, exe_bn)
            shutil.copyfile(old_exe_path, new_exe_path)
            
            # run vvwm.exe (vvwm.exe [...]/outfall_31_xx/vvwmTransfer.txt)
            command = new_exe_path + ' ' + outfall_file
            subprocess.call(command)
            
            # simulated data vs field data units conversion:
            ## Simulation data conc: units = kg/m^3
            ## Observed data conc: units = micrograms/L
            # unit conversion: (meters^3 / liter) & (micrograms / kilogram)
            m3_per_L, mcrg_per_kg = 0.001, 1000000000
            
            # read in produced data from the output of the vvwm run we just completed
            output = pd.read_csv(filelines[68][:-1], 
                                usecols = [1], skiprows=5, names = ["davg_bif_conc"]) * m3_per_L * mcrg_per_kg
            ## Conversion process: Kg/m^3 * 1m^3/1000L * 1000000000 micrograms/kg = 1000000 * 1microgram/L

            # add a date column
            if mode == 'debug':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=103, freq='D')
            elif mode == 'test':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=778, freq='D')
            elif mode == 'run':
                output['Sample_date'] = pd.date_range(start='1/1/2009', periods=3287, freq='D')
            # add a location columns
            output['Site_code'] = o[-5:]
            # sift out only the rows with dates and sites that there is also field data available for
            output = output.merge(obs_data, how = "inner", on = ['Sample_date','Site_code'])
            # change index to a (date, site) pair
            output.set_index([np.datetime_as_string(output.Sample_date, unit = 'D'),'Site_code'], inplace = True)
            # add this df to the df that contains the dfs of all outfalls
            output_df = output_df.append(output[['davg_bif_conc',]], ignore_index = False)

        # label and sort by date and site
        output_df = output_df.set_index([["_".join([a,b[3:]]) for a,b in output_df.index]]).sort_index()


        # In[24]:


        # conver data frame to dictionary
        output_dict = output_df.to_dict()['davg_bif_conc']


        # In[25]:


        # cleanup
        if vvwm_cleanup == 'none' or vvwm_cleanup == 'some' or vvwm_cleanup == 'full':
            exe_ = os.path.join(sdir_path, "outfall_31_??", exe_bn)
            wet_ = os.path.join(sdir_path, "outfall_31_??", "vvwm_wet.dvf")
            loginfo("Deleting vvwm exe and weather files to free up memory.")
            rm(exe_, wet_)
        if vvwm_cleanup == 'some' or vvwm_cleanup == 'full':
            zts_ = os.path.join(sdir_path, "outfall_31_??", "output.zts")
            analysis_ = os.path.join(sdir_path, "outfall_31_??", "output_NPlesant_Custom_parent_analysis.txt")
            transfer_ = os.path.join(sdir_path, "outfall_31_??", "vvwmTransfer.txt")
            loginfo("Deleting internal vvwm files to free up memory.")
            rm(zts_, analysis_, transfer_)
        if vvwm_cleanup == 'full':
            ofdir_ = os.path.join(sdir_path, "outfall_31_??")
            loginfo("Deleting vvwm results files and folder to free up memory.")
            rm(ofdir_)
        if vvwm_cleanup == 'full' and swmm_cleanup == 'full':
            loginfo("Deleting temp folder.")
            rm(sdir_path)


        # In[26]:


        return(output_dict)
        # output_dict


        # In[27]:


    return(model)
