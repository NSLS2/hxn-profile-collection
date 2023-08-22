import numpy as np
import os
import sys
import json

with open("/nsls2/data/hxn/shared/config/bluesky/profile_collection/startup/plot_elems.json", "r") as fp:
    xrf_elems = json.load(fp)
    fp.close()

roi_elems = xrf_elems["roi_elems"]
live_plot_elems  = xrf_elems["live_plot_elems"]
line_plot_elem = xrf_elems["line_plot_elem"]

update_elements = reload_bsui

'''
roi_elems = ['Cu','Ge','W_L','Ti','Si','Cl','Ga','S','Cr','Mn','Fe','Co','Ni','Zn','Pt_L','Au_L'] #limit = 16 
live_plot_elems = roi_elems[:4] # no limts but has to be in the above roi_elems
#live_plot_elems = ['Cu','Ge','Ti','W_L'] # no limts but has to be in the above roi_elems
line_plot_elem = roi_elems[0] #only one element
'''

def create_user_dir(users_name):
    
    month = datetime.now().month
    year = datetime.now().year
    
    month_to_quarter = {1:"Q1",2:"Q1",3:"Q1",4:"Q1",
                        5:"Q2",6:"Q2",7:"Q2",8:"Q2",
                        9:"Q3",10:"Q3",11:"Q3",12:"Q3"}
    year_quarter = f"{year}{month_to_quarter[month]}"
    dir_name = f"/data/users/{year_quarter}/{users_name}_{year_quarter}"

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        return dir_name
    else:
        print("user directory already exists")
        return dir_name

def create_user_symlink(src_dir = "/data/2023Q1/"):
    
    if os.path.exists(src_dir):
        dest = "/data/users/current_user"
        
        #recreate symlink
        if os.path.exists(dest):
            os.unlink(dest)
        os.symlink(src_dir,dest)
        print(f"Successfully created a symlink for {src_dir} as {dest}")
        return 

    else:
        print(f"{src_dir} does not exists; Failed to create a symlink")

        return
    

def setup_new_user(name = "Lastname",experimenters = "HX,NU,SER", sample = "Gold", sample_image = "/data/users/hxn_logo.png"):
    RE.md["PI"] = name
    RE.md["experimenters"] = experimenters
    RE.md["sample"] = sample
    RE.md["scan_name"] = sample+"_1"

    udir = create_user_dir(name)
    print(f"User directory is; {udir}")
    create_user_symlink(src_dir = udir)
    setup_pdf_function(sample_name = sample, experimenters = experimenters, img_to_add = sample_image)
    insertTitle()

from bluesky_queueserver_api import BPlan
from bluesky_queueserver_api.zmq import REManagerAPI
RM = REManagerAPI()
# RM.item_execute((BPlan("fly2d", ["fs", "zebra", "sclr1", "merlin1", "xspress3"], "dssx", -1, 1, 10, "dssy", -1, 1, 10, 0.1)))
# RM.item_add((BPlan("fly2d", ["fs", "zebra", "sclr1", "merlin1", "xspress3"], "dssx", -1, 1, 10, "dssy", -1, 1, 10, 0.3)))


#XF:03IDC-VA{VT:Chm-TCG:2}P-I