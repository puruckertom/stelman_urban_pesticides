#%%
import sys, os

# %%
main_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# %%
master_path = os.path.join(main_path,"master")

# %%
exe_path = os.path.join(master_path,"exe")

# %%
vvwmTransfer_path = os.path.join(master_path,"vvwmTransfer.txt")

# %%
partition_path = os.path.join(master_path,"outfall_components.txt")

# %%
weather_path = os.path.join(master_path,"weather")

# %%
temp_path = os.path.join(main_path,"temp")
if not os.path.exists(temp_path):
    os.mkdir(temp_path)
    print("Temp folder created", "\n")

# %%
obs_path = os.path.join(master_path,"SURF_water_placer_bifenthrin.csv")

# %%
def set_inp_path(mode:str):
    if mode.lower() == 'debug':
        return os.path.join(main_path, "master_debug", "NPlesantCreek.inp")
    elif mode.lower() == 'test':
        return os.path.join(main_path, "master_test", "NPlesantCreek.inp")
    elif mode.lower() == 'run':
        return os.path.join(master_path, "NPlesantCreek.inp")
    else:
        raise ValueError('Acceptable values of mode include \'debug\', \'test\', and \'run\', not \'' + mode + '\'.')

# %%
#print("To activate debug mode, run 'inp_path = set_inp_path('debug')'.\nTo activate test mode, run 'inp_path = set_inp_path('test')'.\nTo activate run mode, run 'inp_path = set_inp_path('run')'.")
