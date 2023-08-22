#Last Update 09/11/2021 by AP

"""

ReadMe:

The workflow for xanes experiment is define below. This macro aims to use one flow for XANES of any given element.
This macro is designed to work with the GUI inputs as well.
To add a new element add the paramer file in the format given below


EXAMPLE OF USAGE:


For XANES Scan: <zp_list_xanes2d(FeXANES,dets1,zpssx,-13,11,150,zpssy,-13,11,150,0.05,
                    xcen = 0, ycen = 0,doAlignScan = False, alignElem = 'Fe',
                    alignX = (-1,1,100,0.1),
                    alignY = (-1,1,100,0.1), pdfElem = ['Fe','Cr'],
                    saveLogFolder = '/data/users/2021Q3/Ajith_2021Q3')


For Foil Calibration: <zp_list_xanes2d(e_list,dets6,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,
                    xcen = 0, ycen = 0, doAlignScan = False, pdfLog = False,
                    foilCalibScan = True, peakBeam = False)

"""


import numpy as np
import pandas as pd
import time,json
from datetime import datetime
import scipy.constants as consts


#Paramer list from previous runs in the order of atomic number of the element

CrXANES = {'high_e':6.0, 'high_e_zpz1':10.48, 'zpz1_slope':-5.04,
          'energy':[(5.97,5.98,0.005),(5.981,6.03,0.001), (6.032,6.046,0.005)] }
          
MnXANES = {'high_e':6.6, 'high_e_zpz1':68.3165, 'zpz1_slope':-5.04,
          'energy':[(6.520,6.530,0.005),(6.531,6.580,0.001),(6.585,6.601,0.005)]}
               
FeXANES = {'high_e':7.2, 'high_e_zpz1':6.41, 'zpz1_slope':-5.04}

NiXANES = {'high_e':8.300, 'high_e_zpz1':0.98, 'zpz1_slope':-5.04,
          'energy':[(8.30,8.325,0.005),(8.326,8.360,0.001),(8.360,8.430,0.006)],}

CuXANES = {'high_e':9.06,  'high_e_zpz1':-4.905, 'zpz1_slope':-5.04,
          'energy':[(8.96,8.975,0.005),(8.976,9.003,0.001)],}

ZnXANES =  {'high_e':9.7, 'high_e_zpz1':50.87, 'zpz1_slope':-5.04,
          'energy':[(9.64,9.666,0.005),(9.6665,9.681,.0005),(9.682,9.701,0.002),(9.705,9.725,0.005)]}

HfXANES =  {'high_e':9.6, 'high_e_zpz1':-7.775, 'zpz1_slope':-5.04,
          'energy':[(9.500,9.540,0.005),(9.541,9.6,0.001)]}

LuL3XANES =  {'high_e':9.3, 'high_e_zpz1':-5.4246, 'zpz1_slope':-5.04,
          'energy':[(9.150,9.200,0.005),(9.201,9.350,0.001),(9.352,9.400,0.002)]}

As_MLL_XANES = {'high_e':11.94, 
                'low_e':11.84,
                'high_e_hmll_z':0,
                'high_e_sbz':0,
                'low_e_hmll_z':9,
                'low_e_sbz':-39,
                'energy':[(11.84,11.86,0.005),
                          (11.861,11.88,0.001),
                          (11.881,11.90,0.002),
                          (11.90,11.94,0.005)]
                          
                }

As_MLL_XANES_minE = {'high_e':11.94, 
                'low_e':11.84,
                'high_e_hmll_z':0,
                'high_e_sbz':0,
                'low_e_hmll_z':9,
                'low_e_sbz':-39,
                'energy':[11.84,11.869,11.870,
                          11.872,11.878,11.880,
                          11.905,11.94]
                          
                }

                                ######################################
                                ######### FUNCTIONS BELOW ############
                                ######################################


def cbpm_on(action = True):
    cbpm_x = "XF:03IDC-CT{FbPid:04}PID:on"
    cbpm_y = "XF:03IDC-CT{FbPid:03}PID:on"

    if action:
        caput(cbpm_x,1)
        caput(cbpm_y,1)

        time.sleep(2)
    else:
        caput(cbpm_x,0)
        caput(cbpm_y,0)
        time.sleep(2)

def piezos_to_zero():
    yield from bps.mov(zpssx,0,zpssy,0,zpssz,0)



def move_energy_and_angle(e,angle,zpz_):

    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")
    yield from mov_zpz1(zpz_)
    yield from bps.mov(zpsth, angle)
    yield from bps.sleep(2)
    #cbpm_on(True)


def create_energy_angle_df(filename, XANESParam = FeXANES):

    """

    Generates a pandas dataframe of optics motor positions. Function uses high E and low E values in the dictionary
    to generate motor positions for all the energy points, assuming linear relationship.

    input: Dictionary conating optics values at 2 positions (high E and low E), option to start from high E or low E

    return : Dataframe looks like below;

       energy     ZP focus
    0   7.175     65.6575
    1   7.170     65.6870
    2   7.165     65.7165
    3   7.160     65.7460
    4   7.155     65.7755

    """
    # empty dataframe
    e_list = pd.DataFrame()
    
    enegy_and_angle = np.loadtxt(filename)
    e_points = enegy_and_angle[:,1]/1000
    angles = enegy_and_angle[:,0]

    #add list of energy as first column to DF
    e_list['energy'] = e_points
    e_list['angle'] = angles
    

    #read the paramer dictionary and calculate ugap list
    high_e = XANESParam['high_e']

    #zone plate increament is very close to the theorticla value , same step as above for zp focus
    zpz1_ref, zpz1_slope = XANESParam['high_e_zpz1'],XANESParam['zpz1_slope']
    zpz1_list = zpz1_ref + (e_list['energy'] - high_e)*zpz1_slope
    e_list['ZP focus'] = zpz1_list

    #return the dataframe
    return e_list

def alignment_scan(mtr, start,end,num,exp,elem_, align_with="line_center", threshold = 0.5, move_coarse = True):

    """
    scan to align samples to field of view using using fly1d scan 

    mtr--> scanning motor, dssx, dssy, dssz etc.
    start,end,num,exp --> flyscan paramters
    elem_ --> element to use for alignemnt
    align_with --> choose bettween "edge" or "line_center"
    threshold --> threshold for line centering
    
    """
    
    fly_to_coarse = {"zpssx":"smarx","zpssy":"smary","zpssz":"smarz"}

    yield from fly1d(dets_fs,
                    mtr, 
                    start, 
                    end, 
                    num,
                    exp
                    )
    if align_with == "line_center":
        xc = return_line_center(-1,elem_,threshold)

    elif align_with == "edge":
        xc,_ = erf_fit(-1,elem_,linear_flag=False)

    else:
        raise KeyError(f"{align_with}  is not defined")
    print(f"{mtr.name} centered to {xc :.2f}")
    
    if move_coarse:
        yield from bps.movr(eval(fly_to_coarse[mtr.name]),xc/1000)
        
    else:
        yield from bps.mov(mtr,xc)

def zp_tomo_2d_scan(angle,dets_,x_start,x_end,x_num,y_start,y_end,y_num,exp):
    print("zp tomo 2d scan")
    
    x_scale_factor = 0.9542
    z_scale_factor = 1.0309

    if np.abs(angle) < 44.99:
                
        x_start_real = x_start / np.cos(angle * np.pi / 180.)/ x_scale_factor
        x_end_real = x_end / np.cos(angle * np.pi / 180.)/ x_scale_factor

        yield from fly2d(dets_, 
                        zpssx,
                        x_start_real,
                        x_end_real,
                        x_num,
                        zpssy,
                        y_start, 
                        y_end, 
                        y_num, 
                        exp
                        )

    else:

        x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))/ z_scale_factor
        x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))/ z_scale_factor
        print(x_start_real,x_end_real)

        yield from fly2d(dets_, 
                        zpssz,
                        x_start_real,
                        x_end_real,
                        x_num,
                        zpssy,
                        y_start, 
                        y_end, 
                        y_num, 
                        exp
                        )

def zp_spectro_tomo_scan(elemParam,path_to_json,pdfElem = ['Fe','Cr'],
                        pdfLog = True, peakBeam = True,
                        saveLogFolder = '/data/users/current_user'):

    """ 
    Function to run XANES Scan. 
    
    Arguments:
           1. elemParam: Dictionary -  containg low and high energy optics positions and other useful info 
           2. dets: list - detector system in use
           3. mot1, mot2: EpicsMotors- Motors used for 2D scanning (eg: zpssx, zpssy, etc)
           4. xs,xe,ys,ye: float - scan start and end positions in X&Y directions
           5. x_num,y_num: float - number of steps in X&Y directions
           6. accq_t: float - aquistion (dwell) time for flyscan
           7. highEStart: boolean - if True start the stack with high energies first (Descenting order)
           8. doAlignScan: boolean - if True registration scans will be performed before the 2D scan
           9. xcen, ycen; positions where alignemnt scan would be done. This number updates after each alignment scan
           10. Options for reginstration scans
           11. Options to save XRFs to pdf after each scan
           12. Options to do foil calibration scans
           13. Save important information in CSV format to selected forlder 
           14. The user can turn on and off alignemnt scans
    
    
    """   
    #load paramfile
    with open(path_to_json,"r") as fp:
        scan_params = json.load(fp)
        fp.close()
    # marker to track beam dump             
    beamDumpOccured = False
                    
    e_list = create_energy_angle_df(scan_params["energy_angle_list_file"], XANESParam = elemParam)

    #add real energy to the dataframe
    e_list['E Readback'] = np.nan 
    
    #add scan id to the dataframe
    e_list['Scan ID'] = np.nan 
    
    #recoed time
    e_list['TimeStamp'] = pd.Timestamp.now()
    
    #Ic values are useful for calibration
    e_list['IC3'] = sclr2_ch4.get() 
    e_list['IC0'] = sclr2_ch2.get()
    e_list['IC3_before_peak'] = sclr2_ch4.get()
    
    
    #record if peak beam happed before the scan   
    e_list['Peak Flux'] = False 
    
    print(e_list.head())
    yield from bps.sleep(1)#time to quit if anything wrong
    
    #get intal ic1 value
    ic_0 = sclr2_ch2.get()
    
    #opening fast shutter for initial ic3 reading
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
    yield from bps.sleep(2)

    logger.info("Reading IC3 value")
    
    #get the initial ic3 reading for peaking the beam
    ic_3_init =  sclr2_ch4.get()
    #open the json file to catch any updates 


    image_scan_i = scan_params["fly2d_scan"]

    tot_time_ = (image_scan_i["x_num"]*image_scan_i["y_num"]*image_scan_i["exposure"]*len(e_list))
    tot_time = tot_time_/3600
    overhead = 1.5
    end_datetime = time.ctime(time.time()+tot_time_*overhead)
    check = input(f"This plan takes about {tot_time*overhead :.1f} hours,"
                    f"Projected to {end_datetime} continue (y/n)?")
    if check == "y":

        for i in tqdm.tqdm(range(len(e_list)),desc = 'Energy-Angle Scan'):
        #for i in range (len(e_list)):

            #open the json file to catch any updates 
            with open(path_to_json,"r") as fp:
                scan_params = json.load(fp)
                fp.close()

            while scan_params["pause_scan"]:
                yield from bps.sleep(10) #check if this freezes the gui or not
                with open(path_to_json,"r") as fp:
                    scan_params = json.load(fp)
                    fp.close() 

                if not scan_params["pause_scan"]:   
                    break

            #stop data collection if necessary.user input taken 
            if scan_params["stop_iter"]:
                save_page()
                break

            #if beam dump occur turn the marker on
            if sclr2_ch2.get()<1000:
                beamDumpOccured = True

            #wait if beam dump occured beamdump
            yield from check_for_beam_dump(threshold=5000)
            
            if beamDumpOccured:
                #wait for about 3 minutes for all the feedbacks to kick in
                yield from bps.sleep(120)

                #redo the previous energy
                e_t,angle_t,zpz_t, *others = e_list.iloc[i-1]

                #turn off the beamdump marker
                beamDumpOccured = False
                
            else:
                #unwrap df row for energy change
                e_t,angle_t,zpz_t, *others = e_list.iloc[i]
            
            yield from move_energy_and_angle(e_t,angle_t,zpz_t)

            #open fast shutter to check if ic3 reading is satistactory
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
            yield from bps.sleep(5)
            
            #get ic3 value before peaking, e change
            ic3_ = sclr2_ch4.get()
            
            # if ic3 value is below the threshold, peak the beam
            if ic3_ < ic_3_init*0.8:
                
                if peakBeam: yield from peak_the_flux()
                fluxPeaked = True # for df record
            else:
                fluxPeaked = False
            
            #for df
            ic_3 = sclr2_ch4.get()
            ic_0 = sclr2_ch2.get()
            
            alignX = scan_params["xalign"]
            alignY = scan_params["yalign"]       

            if e_list['energy'][i]<0: # for special scans if no align elem available
                
                '''
                yield from fly1d(dets,zpssx,-1,1,100,0.1)
                xcen = return_line_center(-1,'Cl',0.7)
                yield from bps.mov(zpssx, xcen)
                yield from fly1d(dets,zpssy,-1,1 ,100,0.1)
                ycen = return_line_center(-1,'Cl',0.7)
                yield from bps.mov(zpssy, ycen)
                '''
                pass

            else:
                if np.abs(angle_t) < 44.99:
                    mtr = zpssx
                else:
                    mtr = zpssz
                
                if alignX["do_align"]:
                    yield from alignment_scan(  mtr, 
                                                alignX["start"],
                                                alignX["end"],
                                                alignX["num"],
                                                alignX["exposure"],
                                                alignX["elem"],
                                                align_with=alignX["center_with"], 
                                                threshold = alignX["threshold"])                

                if alignY["do_align"]:
                    yield from alignment_scan(  zpssy, 
                                                alignY["start"],
                                                alignY["end"],
                                                alignY["num"],
                                                alignY["exposure"],
                                                alignY["elem"],
                                                align_with=alignY["center_with"], 
                                                threshold = alignY["threshold"]
                                                ) 
            # alignment_scan(mtr, start,end,num,exp,elem_, align_with="line_center", threshold = 0.5):
                                                                       

            print(f'Current scan: {i+1}/{len(e_list)}')
            image_scan = scan_params["fly2d_scan"]

            yield from zp_tomo_2d_scan( angle_t,
                                        eval(image_scan["det"]),
                                        image_scan["x_start"],
                                        image_scan["x_end"],
                                        image_scan["x_num"],
                                        image_scan["y_start"],
                                        image_scan["y_end"],
                                        image_scan["y_num"],
                                        image_scan["exposure"]
                                        )
            yield from bps.sleep(1)

            #close fast shutter
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
            # get some scan details and add to the list of scan id and energy
            last_sid = int(caget('XF:03IDC-ES{Status}ScanID-I'))
            e_pos = e.position
            
            #Add more info to the dataframe
            e_list['E Readback'].at[i] = e_pos #add real energy to the dataframe
            e_list['Scan ID'].at[i] = int(last_sid) #add scan id to the dataframe
            e_list['TimeStamp'].at[i] = pd.Timestamp.now()
            e_list['IC3'].at[i] = ic_3 #Ic values are useful for calibration
            e_list['IC0'].at[i] = ic_0 #Ic values are useful for calibration
            e_list['Peak Flux'].at[i] = fluxPeaked # recoed if peakflux was excecuted
            e_list['IC3_before_peak'].at[i] = ic3_ #ic3 right after e change, no peaking
            fluxPeaked = False #reset
            
            if pdfLog:
                try:
                    insert_xrf_map_to_pdf(-1,pdfElem,title_=['energy', 'zpsth'])# plot data and add to pdf
                except:
                    pass
            # save the DF in the loop so quitting a scan won't affect
            filename = f"HXN_nanoXANES_StartID{int(e_list['Scan ID'][0])}_{len(e_list)}_e_points.csv"
            e_list.to_csv(os.path.join(saveLogFolder, filename), float_format= '%.5f')

        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if pdfLog: save_page() #save the pdf

    else:
        return
