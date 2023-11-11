#Last Update 09/11/2021 by AP

"""

ReadMe:



EXAMPLE OF USAGE:

<run_zp_xanes(CuXANES,"/data/current_user/xanes_2d_test.json")


"""


import numpy as np
import pandas as pd
import time,json
from datetime import datetime
import scipy.constants as consts



#Paramer list from previous runs in the order of atomic number of the element

CrXANES = {'high_e':6.0, 'high_e_zpz1':10.48, 'zpz1_slope':-5.04,
          'energy':[(5.97,5.98,0.005),(5.981,6.03,0.001), (6.032,6.046,0.005)] }
          
MnXANES = {'high_e':6.6, 'high_e_zpz1':9.31, 'zpz1_slope':-5.04,
          'energy':[(6.510,6.530,0.005),(6.531,6.570,0.001),(6.575,6.600,0.0075)]}
               
FeXANES = {'high_e':7.2, 'high_e_zpz1':6.41, 'zpz1_slope':-5.04,
          'energy':[(7.09,7.105,0.005),(7.106,7.141,0.001),(7.14,7.18,0.005)],}
          
CoXANES = {'high_e':7.8, 'high_e_zpz1':3.197, 'zpz1_slope':-5.04,          
            'energy':[(7.690,7.705,0.005),(7.706,7.760,0.001),(7.765,7.800,0.005)],}

#CoXANES = {'high_e':7.8, 'high_e_zpz1':3.2725, 'zpz1_slope':-5.04,
#          'energy':[(7.736,7.760,0.001),(7.765,7.800,0.005)],}


NiXANES = {'high_e':8.300, 'high_e_zpz1':0.98, 'zpz1_slope':-5.04,
          'energy':[(8.30,8.325,0.005),(8.326,8.360,0.001),(8.360,8.430,0.006)],}

CuXANES = {'high_e':9.05,  'high_e_zpz1':-2.735, 'zpz1_slope':-5.04,
          'energy':[(8.950,8.975,0.005),(8.976,9.005,0.001),(9.009,9.033,0.004)],}

ZnXANES =  {'high_e':9.7, 'high_e_zpz1':-6.25, 'zpz1_slope':-5.04,
          'energy':[(9.620,9.650,0.005),(9.651,9.700,.001),(9.705,9.725,0.005)]}
          
HfXANES =  {'high_e':9.6, 'high_e_zpz1':-7.775, 'zpz1_slope':-5.04,
          'energy':[(9.500,9.540,0.005),(9.541,9.6,0.001)]}

LuL3XANES =  {'high_e':9.3, 'high_e_zpz1':-5.4246, 'zpz1_slope':-5.04,
          'energy':[(9.150,9.200,0.005),(9.201,9.350,0.001),(9.352,9.400,0.002)]}

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



def move_energy(e,zpz_ ):

    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")
    yield from mov_zpz1(zpz_)
    yield from bps.sleep(2)

    #cbpm_on(True)




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
        yield from piezos_to_zero()
        yield from bps.movr(eval(fly_to_coarse[mtr.name]),xc/1000)
        
    else:
        yield from bps.mov(mtr,xc)

                        

def generateEPoints(ePointsGen = [(9.645,9.665,0.005),(9.666,9.7,0.0006),(9.705,9.725,0.005)],reversed = True):

    """

    Generates a list of energy values from the given list

    input: Tuples in the format (start energy, end energy, energy resolution),
    if reversed is true the list will be transposed

    return : list of energy points

    """

    e_points = []

    if isinstance(ePointsGen[0], tuple):

        for values in ePointsGen:
            #use np.arange to generate values and extend it to the e_points list
            e_points.extend(np.arange(values[0],values[1],values[2]))

    elif isinstance(ePointsGen, list):
        e_points = ePointsGen

    else:
        raise TypeError("Invalid energy format")

    if reversed:
        #retrun list in the reversted order
        return e_points[::-1]
    else:
        return e_points                 
                        
                        
def generateEList(XANESParam = CrXANES, highEStart = True):

    """

    Generates a pandas dataframe of optics motor positions. Function uses high E and low E values in the dictionary
    to generate motor positions for all the energy points, assuming linear relationship.

    input: Dictionary conating optics values at 2 positions (high E and low E), option to start from high E or low E

    return : Dataframe looks like below;

       energy    ugap  crl_theta  ZP focus
    0   7.175  7652.5       1.75   65.6575
    1   7.170  7648.0       1.30   65.6870
    2   7.165  7643.5       0.85   65.7165
    3   7.160  7639.0       0.40   65.7460
    4   7.155  7634.5      -0.05   65.7755

    """
    # empty dataframe
    e_list = pd.DataFrame()

    #add list of energy as first column to DF
    e_list['energy'] = generateEPoints (ePointsGen = XANESParam ['energy'], reversed = highEStart)

    #read the paramer dictionary and calculate ugap list
    high_e = XANESParam['high_e']

    #zone plate increament is very close to the theorticla value , same step as above for zp focus
    zpz1_ref, zpz1_slope = XANESParam['high_e_zpz1'],XANESParam['zpz1_slope']
    zpz1_list = zpz1_ref + (e_list['energy'] - high_e)*zpz1_slope
    e_list['ZP focus'] = zpz1_list

    #return the dataframe
    return e_list

def generateEList_MLL(XANESParam = CuXANES, highEStart = False):
    print("generating e_list")

    """

    Generates a pandas dataframe of optics motor positions. Function uses high E and low E values in the dictionary
    to generate motor positions for all the energy points, assuming linear relationship.

    input: Dictionary conating optics values at 2 positions (high E and low E), option to start from high E or low E

    return : Dataframe looks like below;

       energy    ugap  crl_theta  ZP focus
    0   7.175  7652.5       1.75   65.6575
    1   7.170  7648.0       1.30   65.6870
    2   7.165  7643.5       0.85   65.7165
    3   7.160  7639.0       0.40   65.7460
    4   7.155  7634.5      -0.05   65.7755

    """
    # empty dataframe
    e_list = pd.DataFrame()

    #add list of energy as first column to DF
    e_list['energy'] = generateEPoints (ePointsGen = XANESParam ['energy'], reversed = highEStart)

    #read the paramer dictionary and calculate ugap list
    high_e = XANESParam['high_e']
    low_e = XANESParam['low_e']

    #lens increament

    high_hmll = XANESParam['high_e_hmll_z']
    high_sbz = XANESParam['high_e_sbz']
    low_hmll = XANESParam['low_e_hmll_z']
    low_sbz = XANESParam['low_e_sbz']


    hmll_z_slope = (high_hmll-low_hmll)/(high_e-low_e)
    sbz_slope = (high_sbz-low_sbz)/(high_e-low_e)
    print(sbz_slope)


    hmll_list = high_hmll + (e_list['energy'] - high_e)*hmll_z_slope
    sbz_list = high_sbz + (e_list['energy'] - high_e)*sbz_slope
    
    e_list['hmll_hz'] = hmll_list
    e_list['sbz'] = sbz_list


    #return the dataframe
    return e_list                        

def run_zp_xanes(elemParam,path_to_parameter_file):

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
          "start_from_high_e":false,
      "stop_iter": false,
      "pause_scan": false,
      "pdf_log":true,
      "pdf_elems:["Cu"],
      "peak_flux":true,
      "save_log_to":"/data/users/current_user"
    
    """   
    #load paramfile
    with open(path_to_parameter_file,"r") as fp:
        scan_params = json.load(fp)
        fp.close()
    # marker to track beam dump             
    beamDumpOccured = False    
                    
    e_list = generateEList(elemParam, highEStart =  scan_params["start_from_high_e"])

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
    e_list["zpsth"] = zpsth.position
    
    
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
            with open(path_to_parameter_file,"r") as fp:
                scan_params = json.load(fp)
                fp.close()

            while scan_params["pause_scan"]:
                yield from bps.sleep(10) #check if this freezes the gui or not
                with open(path_to_parameter_file,"r") as fp:
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
                yield from recover_from_beamdump()

                #redo the previous energy
                e_t,angle_t,zpz_t, *others = e_list.iloc[i-1]

                #turn off the beamdump marker
                beamDumpOccured = False
                
            else:
                #unwrap df row for energy change
                e_t,zpz_t, *others = e_list.iloc[i]
            
            yield from move_energy(e_t,zpz_t)

            #open fast shutter to check if ic3 reading is satistactory
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
            yield from bps.sleep(5)
            
            #get ic3 value before peaking, e change
            ic3_ = sclr2_ch4.get()
            
            # if ic3 value is below the threshold, peak the beam
            #if ic3_ < ic_3_init*0.9:
            if ic3_ < ic_3_init*scan_params["flux_threshold"]:
                
                if scan_params["peak_flux"]: yield from peak_the_flux()
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
                
                if alignX["do_align"]:
                    yield from alignment_scan(  zpssx, 
                                                alignX["start"],
                                                alignX["end"],
                                                alignX["num"],
                                                alignX["exposure"],
                                                alignX["elem"],
                                                align_with=alignX["center_with"], 
                                                threshold = alignX["threshold"])                


            # alignment_scan(mtr, start,end,num,exp,elem_, align_with="line_center", threshold = 0.5):
            
                                                                       

            print(f'Current scan: {i+1}/{len(e_list)}')
            image_scan = scan_params["fly2d_scan"]

            yield from fly2d( 
                            eval(image_scan["det"]),
                            zpssx,
                            image_scan["x_start"],
                            image_scan["x_end"],
                            image_scan["x_num"],
                            zpssy,
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
            e_list["zpsth"].at[i] = zpsth.position
            e_list['Scan ID'].at[i] = int(last_sid) #add scan id to the dataframe
            e_list['TimeStamp'].at[i] = pd.Timestamp.now()
            e_list['IC3'].at[i] = ic_3 #Ic values are useful for calibration
            e_list['IC0'].at[i] = ic_0 #Ic values are useful for calibration
            e_list['Peak Flux'].at[i] = fluxPeaked # recoed if peakflux was excecuted
            e_list['IC3_before_peak'].at[i] = ic3_ #ic3 right after e change, no peaking
            fluxPeaked = False #reset
            
            if scan_params["pdf_log"]:
                try:
                    insert_xrf_map_to_pdf(-1,scan_params["pdf_elems"],title_=['energy', 'zpsth'])# plot data and add to pdf
                except:
                    pass
            # save the DF in the loop so quitting a scan won't affect
            filename = f"HXN_2d_xanes_startID{int(e_list['Scan ID'][0])}_{len(e_list)}_e_angle_points.csv"
            e_list.to_csv(os.path.join(scan_params["save_log_to"], filename), float_format= '%.5f')

        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if scan_params["pdf_log"]: save_page() #save the pdf

    else:
        return
