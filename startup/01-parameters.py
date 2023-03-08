import numpy as np
import os


roi_elems = ['Au_M','Si','Cl','P','Ti','Sr_L','La_L','Cr','Mn','Fe','Co','Ni','Cu','Zn','Pt_L','K'] #limit = 16 
live_plot_elems = ['Sr_L','Ti','Mn','La_L'] # no limts but has to be in the above roi_elems
line_plot_elem = ["Mn"] #only one element


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
    
    udir = create_user_dir(name)
    print(f"User directory is; {udir}")
    create_user_symlink(src_dir = udir)
    setup_pdf_function(sample_name = sample, experimenters = experimenters, img_to_add = sample_image)
    insertTitle()


'''
merlin_pos = {'diff_x': 8, 'diff_y1':-12.9, 'diff_y2':-12.9, 'diff_z': -50, 'diff_cz': -24.7}
cam11_pos = {'diff_x': 216.33, 'diff_y1': 19.177, 'diff_y2': 19.177, 'diff_z': -50, 'diff_cz': -24.7}
telescope_pos = {'diff_x': -342,'diff_z': -50, 'diff_cz': -24.7}
'''