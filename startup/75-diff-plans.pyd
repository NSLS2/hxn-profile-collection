#Last Update 09/11/2021 by AP

"""

ReadMe:



EXAMPLE OF USAGE:


"""


import numpy as np
import pandas as pd
import time,json
import scipy
from datetime import datetime
import scipy.constants as consts
from scipy.sparse.linalg import gmres, lgmres, LinearOperator

#from phantominator import shepp_logan
import tifffile as tf



def make_zp_diff_plan(save_as = "/data/users/current_user/zp_diff_params_template.json" ):

    zp_diff_scan = {   
                    "angle_info":{'th_motor':'zpsth',
                                "start":70, 
                                "end":72, 
                                "angle_step":0.05}, 

                    "fly2d_scan":{'det':'dets4',
                                "x_motor":'zpssx',
                                "x_start":-1,
                                "x_end":1,
                                "x_num":100, 
                                "y_motor":'zpssy',
                                "y_start":-1,
                                "y_end":1,
                                "y_num":100,  
                                "exposure":0.03},

                    "xalign":{"do_align":True,
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5,
                            "move_coarse":True,
                            "negative_flag":True},
                    
                    "yalign":{"do_align":True,
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5,
                            "move_coarse":True,
                            "negative_flag":True},

                    "align_2d_com":{"do_align":False,
                            "x_start":-2,
                            "x_end": 2,
                            "x_num": 100,
                            "y_start":-2,
                            "y_end": 2,
                            "y_num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "threshold": 0.5,
                            "move_x":True,
                            "move_y":True},
                    
                    "stop_iter":False,
                    "add_angles":[],
                    "remove_angles":[-90,-91],
                    "stop_pdf":False,
                    "pdf_elems":["Ni"],
                    "pause_scan":False,
                    "test":False,
                    "ic_threshold":0.9,
                    "scan_label":"HXN_diff_Scan",
                    "save_log_to":"/data/users/current_user/"

                }


    with open(save_as,"w") as fp:
            json.dump(zp_diff_scan,fp, indent=6)

    fp.close()
    
    print(f"{save_as} plan is created")


def make_mll_diff_plan(save_as = "/data/users/current_user/mll_diff_params_template.json" ):

    zp_diff_scan = {   
                    "angle_info":{'th_motor':'dsth',
                                "start":70, 
                                "end":72, 
                                "angle_step":0.05}, 

                    "fly2d_scan":{'det':'dets4',
                                "x_motor":'dssx',
                                "x_start":-1,
                                "x_end":1,
                                "x_num":100, 
                                "y_motor":'dssy',
                                "y_start":-1,
                                "y_end":1,
                                "y_num":100,  
                                "exposure":0.03},

                    "xalign":{"do_align":True,
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5,
                            "move_coarse":True,
                            "negative_flag":True},
                    
                    "yalign":{"do_align":True,
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5,
                            "move_coarse":True,
                            "negative_flag":True},

                    "align_2d_com":{"do_align":False,
                            "x_start":-2,
                            "x_end": 2,
                            "x_num": 100,
                            "y_start":-2,
                            "y_end": 2,
                            "y_num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "threshold": 0.5,
                            "move_x":True,
                            "move_y":True},
                    
                    "stop_iter":False,
                    "add_angles":[],
                    "remove_angles":[-90,-91],
                    "stop_pdf":False,
                    "pdf_elems":["Ni"],
                    "pause_scan":False,
                    "test":False,
                    "ic_threshold":0.9,
                    "scan_label":"HXN_diff_Scan",
                    "save_log_to":"/data/users/current_user/"

                }


    with open(save_as,"w") as fp:
            json.dump(zp_diff_scan,fp, indent=6)

    fp.close()
    
    print(f"{save_as} plan is created")




        
                        
def zp_diff_scan_to_loop(angle, diff_params, ic_init):

        #caput("XF:03IDC-ES{Merlin:2}HDF1:NDArrayPort","ROI1") #patch for merlin2 issuee
        
        #get parameters from json
        xalign = diff_params["xalign"]
        yalign = diff_params["yalign"]
        align_2d = diff_params["align_2d_com"]
        image_scan = diff_params["fly2d_scan"]
        dets = eval(image_scan["det"])
        x_motor = eval(image_scan["x_motor"])
        y_motor = eval(image_scan["y_motor"])
        th_motor = eval(diff_params["angle_info"]["th_motor"])
        elems_to_pdf = diff_params["pdf_elems"]


        yield from bps.mov(th_motor, angle)
        
        #look for beam dump and ic3 threshold, ignores for code tests using json
        if not diff_params["test"]:
       
            yield from check_for_beam_dump()

            while (sclr2_ch2.get() < (diff_params["ic_threshold"]*ic_init)):
                 yield from bps.movr(smary,-0.040)
                 yield from peak_the_flux()
                 yield from bps.movr(smary,0.040)
                 ic_0 = sclr2_ch2.get()


        if xalign["do_align"]:
            yield from align_scan(x_motor, 
                            xalign["start"], 
                            xalign["end"], 
                            xalign["num"], 
                            xalign["exposure"],
                            xalign["elem"],
                            align_with = xalign["center_with"],
                            threshold = xalign["threshold"],
                            move_coarse = xalign["move_coarse"],
                            neg_flag =xalign["negative_flag"]
                            )
                            
        if yalign["do_align"]:
            yield from align_scan(  y_motor, 
                                    yalign["start"], 
                                    yalign["end"], 
                                    yalign["num"], 
                                    yalign["exposure"],
                                    yalign["elem"],
                                    align_with = yalign["center_with"],
                                    threshold = yalign["threshold"],
                                    move_coarse = xalign["move_coarse"],
                                    neg_flag =xalign["negative_flag"]
                    )

        #2d alignemnt using center of mass if condition is true
        elif align_2d["do_align"]:

            yield from align_2d_com_scan(x_motor,
                                         align_2d['x_start'],
                                         align_2d['x_end'],
                                         align_2d["x_num"],
                                         y_motor,
                                         align_2d["y_start"], 
                                         align_2d["y_end"], 
                                         align_2d["y_num"], 
                                         align_2d["exposure"],
                                         align_2d["elem"],
                                         align_2d["threshold"],
                                         align_2d["move_x"],
                                         align_2d["move_y"],)

        else:
            pass
        

            yield from fly2d( 
                            eval(image_scan["det"]),
                            x_motor,
                            image_scan["x_start"],
                            image_scan["x_end"],
                            image_scan["x_num"],
                            y_motor,
                            image_scan["y_start"],
                            image_scan["y_end"],
                            image_scan["y_num"],
                            image_scan["exposure"]
                            )

        #save images to pdf if
        if not diff_params["stop_pdf"]:

            insert_xrf_map_to_pdf(-1,elements=elems_to_pdf, 
                                  title_ = ['energy','zpsth'], 
                                  note = diff_params["scan_label"])
            plt.close()


def run_zp_diff(path_to_json):
    beamDumpOccured = False
                    
    #open json file for angle info first
    with open(path_to_json,"r") as fp:
        diff_params = json.load(fp)
    fp.close()
    print("json file loaded")

    #create angle list for iteration
    angle_info = diff_params["angle_info"]
    print(angle_info)

    angles = np.linspace(angle_info["start"], 
                        angle_info["end"],
                        int(1+abs(angle_info["end"] - angle_info["start"])/angle_info["angle_step"])
                        )
                        
    print(f"total angles  = {len(angles)}")
                        
    angle_list = pd.DataFrame()

    angle_list["angles"] = angles

    #add real energy to the dataframe
    angle_list['E Readback'] = np.nan 
    
    #add scan_id to the dataframe
    angle_list['scan_id'] = np.nan 
    
    #recoed time
    angle_list['TimeStamp'] = pd.Timestamp.now()
    
    
    #record if peak beam happed before the scan   
    angle_list['Peak Flux'] = False 
    
    print(angle_list)
    yield from bps.sleep(1)#time to quit if anything wrong
    
    #get intal ic1 value
    ic_0 = sclr2_ch2.get()
    
    #opening fast shutter for initial ic3 reading
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
    yield from bps.sleep(5)

    logger.info("Reading IC3 value")
    
    #get the initial ic3 reading for peaking the beam
    ic_3_init =  sclr2_ch4.get()
    #open the json file to catch any updates 
    
     #Ic values are useful for calibration
    angle_list['IC3'] = ic_3_init
    angle_list['IC0'] = sclr2_ch2.get()
    angle_list['IC3_before_peak'] = ic_3_init
    angle_list["zpsth"] = np.nan
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 

    image_scan_i = diff_params["fly2d_scan"]

    tot_time_ = (image_scan_i["x_num"]*image_scan_i["y_num"]*image_scan_i["exposure"]*len(angle_list))
    tot_time = tot_time_/3600
    overhead = 1.5
    end_datetime = time.ctime(time.time()+tot_time_*overhead)
    check = input(f"This plan takes about {tot_time*overhead :.1f} hours,"
                    f"Projected to {end_datetime} continue (y/n)?")
    if check == "y":

        #loop with list of angles
        for n,angle in enumerate(tqdm.tqdm(angles,desc = 'Diff Scan')):
            print(f"{angle = }")
            yield from bps.sleep(1)

            #open the json file to catch any updates 
            with open(path_to_json,"r") as fp:
                diff_params = json.load(fp)
                fp.close()
                
                    #stop data collection if necessary.user input taken 
            if diff_params["stop_iter"]:
                save_page()
                break
        

            while diff_params["pause_scan"]:
                yield from bps.sleep(10) #check if this freezes the gui or not
                with open(path_to_json,"r") as fp:
                    diff_params = json.load(fp)
                    fp.close() 

                if not diff_params["pause_scan"]:   
                    break
                    
            if diff_params["remove_angles"]==None:
                diff_params["remove_angles"] = []
            

            if sclr2_ch2.get()<1000:
                beamDumpOccured = True
                yield from check_for_beam_dump()
                
                
            if beamDumpOccured:
                angle = angles[n-1]
                yield from bps.sleep(360) #time for uofb kick off
                yield from recover_from_beamdump()
                beamDumpOccured = False
                        
            #open fast shutter to check if ic3 reading is satistactory
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
            yield from bps.sleep(3)
                        
            #for df
            ic_3 = sclr2_ch4.get()
            ic_0 = sclr2_ch2.get()
                
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
            
            if not angle in np.array(diff_params["remove_angles"]):
                #diff scan at a single angle
                yield from zp_diff_scan_to_loop(angle, diff_params,ic_0)

            else:
                print(f"{angle} skipped")
                pass
            
            last_sid = int(caget('XF:03IDC-ES{Status}ScanID-I'))
            
            #Add more info to the dataframe
            angle_list['E Readback'].at[n] = e.position #add real energy to the dataframe
            angle_list["zpsth"].at[n] = zpsth.position
            angle_list['scan_id'].at[n] = int(last_sid) #add scan_id to the dataframe
            angle_list['TimeStamp'].at[n] = pd.Timestamp.now()
            angle_list['IC3'].at[n] = ic_3 #Ic values are useful for calibration
            angle_list['IC0'].at[n] = ic_0 #Ic values are useful for calibration
            #angle_list['Peak Flux'].at[n] = fluxPeaked # recoed if peakflux was excecuted
            #angle_list['IC3_before_peak'].at[n] = ic3 #ic3 right after e change, no peaking
            #fluxPeaked = False #reset
            
            # save the DF in the loop so quitting a scan won't affect
            filename = f"hxn_zp_diff_{diff_params.get('scan_label','')}_startID{int(angle_list['scan_id'][0])}.csv"
            angle_list.to_csv(os.path.join(diff_params["save_log_to"], filename), float_format= '%.5f')
            
            
            
        #TODO add angles to scan; need to be better
        #sort based on what current angle is
        if not diff_params["add_angles"]==None:
    
            added_angles = diff_params["add_angles"]
        
        for angle in tqdm.tqdm(added_angles,desc = 'Diff Scan; Additional Angles'):
            
            #open the json file to catch any updates 
            with open(path_to_json,"r") as fp:
                diff_params = json.load(fp)
                fp.close()
                
            if sclr2_ch2.get()<1000:
                beamDumpOccured = True
                yield from check_for_beam_dump()
                
                
            if beamDumpOccured:
                angle = angles[n-1]
                yield from bps.sleep(360) #time for uofb kick off
                yield from recover_from_beamdump()
                beamDumpOccured = False    

            #stop data collection if necessary.user input taken 
            if diff_params["stop_iter"]:
                save_page()
                break

            while diff_params["pause_scan"]:
                yield from bps.sleep(10) #check if this freezes the gui or not
                with open(path_to_json,"r") as fp:
                    diff_params = json.load(fp)
                    fp.close() 

                if not diff_params["pause_scan"]:   
                    break
            
            if not angle in np.array(diff_params["remove_angles"]):
                yield from zp_diff_scan_to_loop(angle, diff_params,ic_0)

            else:
                print(f"{angle} skipped")
                pass
                
                
                        #Add more info to the dataframe
            angle_list['angle'] = angle
            angle_list['E Readback'].at[n] = e_pos #add real energy to the dataframe
            angle_list["zpsth"].at[n] = zpsth.position
            angle_list['scan_id'].at[n] = int(last_sid) #add scan_id to the dataframe
            angle_list['TimeStamp'].at[n] = pd.Timestamp.now()
            angle_list['IC3'].at[n] = ic_3 #Ic values are useful for calibration
            angle_list['IC0'].at[n] = ic_0 #Ic values are useful for calibration

            
            # save the DF in the loop so quitting a scan won't affect
            angle_list.to_csv(os.path.join(diff_params["save_log_to"], filename), float_format= '%.5f')

        else:
            pass

        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        save_page() #save the pdf

    else:
        return
