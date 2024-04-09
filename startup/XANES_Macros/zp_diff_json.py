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

                    "fly2d_scan":{'det':'dets1',
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
                    "pause_scan":False,
                    "test":False,
                    "ic_threshold":0.9,
                    "scan_label":"HXN_diff_Scan"

                }


    with open(save_as,"w") as fp:
            json.dump(zp_diff_scan,fp, indent=6)

    fp.close()
    
    print(f"{save_as} plan is created")


def piezos_to_zero():
    yield from bps.mov(zpssx,0,zpssy,0,zpssz,0)


def align_scan(mtr,start,end,num,exp,elem_, align_with="line_center", 
               threshold = 0.5,move_coarse = True, neg_flag = False):

    """
    scan to align samples to field of view using using fly1d scan 

    mtr--> scanning motor, dssx, dssy, dssz etc.
    start,end,num,exp --> flyscan paramters
    elem_ --> element to use for alignemnt
    align_with --> choose bettween "edge" or "line_center"
    threshold --> threshold for line centering
    
    """
    fly_to_coarse = {"zpssx":"smarx","zpssy":"smary","zpssz":"smarz",
                     "dssx":"dsx","dssy":"dsy","dssz":"dsz"}
    uni_conv = 1

    yield from fly1d(dets_fs,
                    mtr, 
                    start, 
                    end, 
                    num,
                    exp
                    )
    if align_with == "line_center":
        xc = return_line_center(-1,elem_,threshold, neg_flag = neg_flag)

    elif align_with == "edge":
        xc,_ = erf_fit(-1,elem_,linear_flag=False)

    else:
        xc = mtr.position
        
    if mtr.name.startswith('zp'):
        uni_conv = 1000

    if move_coarse:
        yield from bps.movr(eval(fly_to_coarse[mtr.name]),xc/uni_conv)
        yield from piezos_to_zero()
        
    else:
        yield from bps.mov(mtr,xc)
                        
                        
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
        th_motor = eval(diff_params["anle_info"]["th_motor"])
        elems_to_pdf = diff_params["pdf_elems"]


        yield from bps.mov(th_motor, angle)
        
        #look for beam dump and ic3 threshold, ignores for code tests using json
        if not diff_params["test"]:
       
            yield from check_for_beam_dump()

            while (sclr2_ch2.get() < (diff_params["ic_threshold"]*ic_init)):
                 yield from peak_the_flux()
                 ic_0 = sclr2_ch2.get()


        if xalign["do_align"]:
            yield from align_scan(x_motor, 
                            xalign["start"], 
                            xalign["end"], 
                            xalign["num"], 
                            xalign["exposure"],
                            xalign["elem"],
                            xalign["center_with"],
                            xalign["threshold"],
                            xalign["move_coarse"]
                            )
                            
        if yalign["do_align"]:
            yield from align_scan(  y_motor, 
                                    yalign["start"], 
                                    yalign["end"], 
                                    yalign["num"], 
                                    yalign["exposure"],
                                    yalign["elem"],
                                    yalign["center_with"],
                                    yalign["threshold"],
                                    xalign["move_coarse"]
                    )

        #2d alignemnt using center of mass if condition is true
        elif align_2d["do_align"]:

            yield from align_2d_com_scan(x_motor,
                                         x_start_real,
                                         x_end_real,
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
        


        #2d scan sequence, based on angle x or z are scanned
        yield from fly2d(dets_, 
                        x_motor,
                        x_start_real,
                        x_end_real,
                        x_num,
                        y_motor,
                        y_start, 
                        y_end, 
                        y_num, 
                        exp
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
    
    #add scan id to the dataframe
    angle_list['Scan ID'] = np.nan 
    
    #recoed time
    angle_list['TimeStamp'] = pd.Timestamp.now()
    
    
    #record if peak beam happed before the scan   
    angle_list['Peak Flux'] = False 
    
    print(angle_list.head())
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


    image_scan_i = diff_params["fly2d_scan"]

    tot_time_ = (image_scan_i["x_num"]*image_scan_i["y_num"]*image_scan_i["exposure"]*len(angle_list))
    tot_time = tot_time_/3600
    overhead = 1.5
    end_datetime = time.ctime(time.time()+tot_time_*overhead)
    check = input(f"This plan takes about {tot_time*overhead :.1f} hours,"
                    f"Projected to {end_datetime} continue (y/n)?")
    if check == "y":

        #loop with list of angles
        for n,angle in enumerate(tqdm.tqdm(angles,desc = 'ZP diff Scan')):
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
            angle_list['Scan ID'].at[n] = int(last_sid) #add scan id to the dataframe
            angle_list['TimeStamp'].at[n] = pd.Timestamp.now()
            angle_list['IC3'].at[n] = ic_3 #Ic values are useful for calibration
            angle_list['IC0'].at[n] = ic_0 #Ic values are useful for calibration
            #angle_list['Peak Flux'].at[n] = fluxPeaked # recoed if peakflux was excecuted
            #angle_list['IC3_before_peak'].at[n] = ic3 #ic3 right after e change, no peaking
            #fluxPeaked = False #reset
            
            # save the DF in the loop so quitting a scan won't affect
            filename = f"hxn_zp_diff_{diff_params.get('scan_label','')}_startID{int(angle_list['Scan ID'][0])}.csv"
            angle_list.to_csv(os.path.join(diff_params["save_log_to"], filename), float_format= '%.5f')
            
            
            
        #TODO add angles to scan; need to be better
        #sort based on what current angle is
        if not diff_params["add_angles"]==None:
    
            added_angles = diff_params["add_angles"]
        
        for angle in tqdm.tqdm(added_angles,desc = 'MLL diff Scan; Additional Angles'):
            
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
            angle_list['Scan ID'].at[n] = int(last_sid) #add scan id to the dataframe
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
