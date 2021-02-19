#%%
import sys, os

# %%
main_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# %%
master_path = os.path.join(main_path,"master")

# %%
exe_path = os.path.join(master_path,"exe")

# %%
#inp_path = os.path.join(master_path,"NPlesantCreek.inp")

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

# # %%
# def test_mode(on = True):
#     if not isinstance(on, bool):
#         raise TypeError("expected {0}, received {1}".format(bool, type(on)))
#     else:
#         global inp_path
#         inp_path = os.path.join(main_path, "master_test", "NPlesantCreek.inp") if on else os.path.join(master_path, "NPlesantCreek.inp")
#         return

# # %%
# def test_mode_on():
#     global inp_path
#     inp_path = os.path.join(main_path, "master_test", "NPlesantCreek.inp")

# # %%
# def test_mode_off():
#     global inp_path
#     inp_path = os.path.join(master_path, "NPlesantCreek.inp")

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
# print("To activate test mode, run 'test_mode()'.\nFrom test mode, run 'test_mode(False)' to deactivate.")
print("To activate debug mode, run 'inp_path = set_inp_path('debug')'.\nTo activate test mode, run 'inp_path = set_inp_path('test')'.\nTo activate run mode, run 'inp_path = set_inp_path('run')'.")
