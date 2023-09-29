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
                    
<zp_list_xanes2d(CuXANES, dets_fs, zpssx,-3, 3, 100,zpssy,-2,2,100,0.025, highEStart=False)


"""


import numpy as np
import pandas as pd
import time
from datetime import datetime
import scipy.constants as consts






#Paramer list from previous runs in the order of atomic number of the element

CrXANES = {'high_e':6.0, 'high_e_zpz1':10.48, 'zpz1_slope':-5.04,
          'energy':[(5.97,5.98,0.005),(5.981,6.03,0.001), (6.032,6.046,0.005)] }
          
MnXANES = {'high_e':6.6, 'high_e_zpz1':9.31, 'zpz1_slope':-5.04,
          'energy':[(6.510,6.530,0.005),(6.531,6.570,0.001),(6.575,6.600,0.005)]}
               
FeXANES = {'high_e':7.2, 'high_e_zpz1':6.41, 'zpz1_slope':-5.04,
          'energy':[(7.09,7.105,0.005),(7.106,7.141,0.001),(7.14,7.20,0.005)],}
          
CoXANES = {'high_e':7.8, 'high_e_zpz1':3.197, 'zpz1_slope':-5.04,          
            'energy':[(7.690,7.705,0.005),(7.706,7.760,0.001),(7.765,7.800,0.005)],}

#CoXANES = {'high_e':7.8, 'high_e_zpz1':3.2725, 'zpz1_slope':-5.04,
#          'energy':[(7.736,7.760,0.001),(7.765,7.800,0.005)],}


NiXANES = {'high_e':8.300, 'high_e_zpz1':0.98, 'zpz1_slope':-5.04,
          'energy':[(8.30,8.325,0.005),(8.326,8.360,0.001),(8.360,8.430,0.006)],}

CuXANES = {'high_e':9.05,  'high_e_zpz1':-2.735, 'zpz1_slope':-5.04,
          'energy':[(8.960,8.975,0.005),(8.976,9.010,0.001),(9.0125,9.05,0.004)],}

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


def peak_the_flux():
    
    #yield from peak_xy_volt(1) #temp


    #cbpm_on(True)

    """ Scan the c-bpm set points to find IC3 maximum """

    print("IC3is below threshold; Peaking the beam.")
    yield from bps.sleep(2)
    yield from peak_bpm_y(-4,4,10)
    yield from bps.sleep(1)
    yield from peak_bpm_x(-10,10,6)
    yield from bps.sleep(1)
    yield from peak_bpm_y(-2,2,4)
    
    

def move_energy(e,zpz_ ):

    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")
    yield from mov_zpz1(zpz_)
    yield from bps.sleep(2)

    #cbpm_on(True)

def move_energy_mll(e, hmll_z = 0,sbz_=0):

    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")

    sbz_pos_ = sbz.position
    if abs(sbz_- sbz_pos_)>1:
        yield from bps.mov(sbz,np.round(sbz_,1))

    hmll_hz_pos = hmll.hz.position
    if abs(hmll_z-hmll_hz_pos)>1:
        yield from bps.mov(hmll.hz,np.round(hmll_z,1))

    else:
        pass



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

def generateEList_MLL(XANESParam = As_MLL_XANES, highEStart = False):
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

def zp_list_xanes2d(elemParam,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = False,
                	doAlignScan = True, alignX = (-15,15,100,0.1,'Cu',0.2,True),
                	alignY = (-15,15,100,0.1,'Cu',0.2, True), xy_offset = (0,0),
                	pdfElem = ['Cu'],doScan = True, moveOptics = True,pdfLog = True,
                	foilCalibScan = False, peakBeam = True,
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
    # marker to track beam dump             
    beamDumpOccured = False
                    
    e_list = generateEList(elemParam, highEStart =  highEStart)

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
    
    #remeber the start positions
    mot1_i = mot1.position
    mot2_i = mot2.position

    tot_time_ = (x_num*y_num*accq_t*len(e_list))
    tot_time = tot_time_/3600
    
    if doAlignScan:
        overhead = 1.5
    else:
        overhead = 1.25

    end_datetime = time.ctime(time.time()+tot_time_*overhead)

    check = input(f"This plan takes about {tot_time*overhead :.1f} hours,"
                    f"Projected to {end_datetime} continue (y/n)?")

    if check == "y":

        for i in tqdm.tqdm(range(len(e_list)),desc = 'Energy Scan'):
        #for i in range (len(e_list)):

            #if beam dump occur turn the marker on
            if sclr2_ch2.get()<1000:
                beamDumpOccured = True
                cbpm_on(False)

            #wait if beam dump occured beamdump
            yield from check_for_beam_dump(threshold=5000)
            
            if beamDumpOccured:
                #wait for about 3 minutes for all the feedbacks to kick in
                yield from bps.sleep(120)

                #redo the previous energy
                e_t, zpz_t, *others = e_list.iloc[i-1]

                #turn off the beamdump marker
                beamDumpOccured = False
                
            else:
                #unwrap df row for energy change
                e_t, zpz_t, *others = e_list.iloc[i]
            
            yield from move_energy(e_t,zpz_t)

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

            # move to particle location for alignemnt scan
            #if doAlignScan:
            
                #yield from bps.mov(zpssx, xcen)
                #yield from bps.mov(zpssy, ycen)
            
            #do the alignemnt scan on the xanes elem after it excited , 
            #otherwise skip or use another element

        

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

            elif doAlignScan:

                try:
                
                    if alignY[-1]:
                        yield from fly1d(dets_fs,zpssy,alignY[0],alignY[1],alignY[2],alignY[3])
                        ycen = return_line_center(-1,alignX[4],alignY[5])
                        #ycen,_ = erf_fit(-1,alignX[4])
                        yield from bps.movr(smary, (ycen)*0.001)
                        #yield from bps.mov(zpssy, ycen)
                        print(f"zpssy centered to {ycen}")
                        plt.close()
                        yield from piezos_to_zero()

                    if alignX[-1]:
                        yield from fly1d(dets_fs,zpssx,alignX[0],alignX[1],alignX[2],alignX[3])
                        xcen = return_line_center(-1,alignX[4],alignX[5])
                        yield from bps.movr(smarx, (xcen)*0.001)
                        #xcen,_ = erf_fit(-1,alignX[4])
                        #yield from bps.mov(zpssx, xcen)
                        print(f"zpssx centered to {xcen}")
                        plt.close()
                        yield from piezos_to_zero()




                except:
                    pass
                
                #yield from bps.movr(smarx,xy_offset[0]/1000)
                #yield from bps.movr(smary,xy_offset[1]/1000)

            print(f'Current scan: {i+1}/{len(e_list)}')

            # do the fly2d scan
            #cbpm_on(False)

            if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

                if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.002) 
                #dead_time = 0.001 for 0.015 dwell

            else:

                if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
            yield from bps.sleep(1)

            #cbpm_on(True)
            yield from piezos_to_zero()

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

        #go back to max energy point if scans done reverese
        max_e_id = e_list['energy'].idxmax()
        e_max, zpz_max, *others = e_list.iloc[max_e_id]
        
        if not np.isclose(e_list['energy'].max(), e.position):
        
            yield from move_energy(e_max,zpz_max)
            
            yield from peak_the_flux()

        
        else: pass
            
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if pdfLog: save_page() #save the pdf

    else:
        return


def mll_list_xanes2d(elemParam,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = False,
                    doAlignScan = True, alignX = (-2.5,2.5,100,0.05,'Fe',0.5,True),
                    alignY = (-2.5,2.5,100,0.05,'Fe',0.5, True), xy_offset = (0,0),
                    pdfElem = ['Fe','As'],pdfLog = True,peakBeam = True,
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
    # marker to track beam dump             
    beamDumpOccured = False
                    
    e_list = generateEList_MLL(elemParam, highEStart =  highEStart)

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
    
    #remeber the start positions
    mot1_i = mot1.position
    mot2_i = mot2.position

    #calculate approximate time
    tot_time_ = (x_num*y_num*accq_t*len(e_list))
    tot_time = tot_time_/3600
    
    if doAlignScan:
        overhead = 1.5
    else:
        overhead = 1.25

    end_datetime = time.ctime(time.time()+tot_time_*overhead)

    check = input(f"This plan takes about {tot_time*overhead :.1f} hours,"
                    f"Projected to {end_datetime} continue (y/n)?")

    if check == "y":

        for i in tqdm.tqdm(range(len(e_list)),desc = 'Energy Scan'):
        #for i in range (len(e_list)):

            #if beam dump occur turn the marker on
            if sclr2_ch2.get()<1000:
                beamDumpOccured = True
                cbpm_on(False)

            #wait if beam dump occured beamdump
            yield from check_for_beam_dump(threshold=5000)
            
            if beamDumpOccured:
                #wait for about 3 minutes for all the feedbacks to kick in
                yield from bps.sleep(120)

                #redo the previous energy
                e_t, hmll_hz_t,sbz_t, *others = e_list.iloc[i-1]

                #turn off the beamdump marker
                beamDumpOccured = False
                
            else:
                #unwrap df row for energy change
                e_t, hmll_hz_t,sbz_t, *others = e_list.iloc[i]
            

            yield from move_energy_mll(e_t,hmll_z = hmll_hz_t,sbz_ = sbz_t)

            
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

            elif doAlignScan:

                try:

                    if alignX[-1]:
                        yield from fly1d(dets_fs,mot1,alignX[0],alignX[1],alignX[2],alignX[3])
                        xcen = return_line_center(-1,alignX[4],alignX[5])
                        #xcen,_ = erf_fit(-1,alignX[4])
                        yield from bps.mov(mot1, xcen)
                        print(f"mot1 centered to {xcen}")
                        plt.close()


                    if alignY[-1]:
                        yield from fly1d(dets_fs,mot2,alignY[0],alignY[1],alignY[2],alignY[3])
                        ycen = return_line_center(-1,alignX[4],alignY[5])
                        #ycen,_ = erf_fit(-1,alignX[4])
                        yield from bps.mov(mot2, ycen)
                        print(f"mot2 centered to {ycen}")
                        plt.close()



                except:
                    pass

            print(f'Current scan: {i+1}/{len(e_list)}')

            # do the fly2d scan
            #cbpm_on(False)

            if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

                yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.002) 
                #dead_time = 0.001 for 0.015 dwell

            else:

                yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
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
                    insert_xrf_map_to_pdf(-1,pdfElem,title_=['energy', 'sbz', 'hz'])# plot data and add to pdf
                except:
                    pass
            # save the DF in the loop so quitting a scan won't affect
            filename = f"HXN_nanoXANES_StartID{int(e_list['Scan ID'][0])}_{len(e_list)}_e_points.csv"
            e_list.to_csv(os.path.join(saveLogFolder, filename), float_format= '%.5f')
        '''
        #go back to max energy point if scans done reverese
        max_e_id = e_list['energy'].idxmax()
        e_max, hmll_max,v_mll_max *others = e_list.iloc[max_e_id]
        
        if not np.isclose(e_list['energy'].max(), e.position):
        
            yield from move_energy_mll(e_max,hmll_max,v_mll_max)
            
            yield from peak_the_flux()

        
        else: pass
        '''   
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if pdfLog: save_page() #save the pdf

    else:
        return


def mll_list_xanes2d_no_input(elemParam,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = False,
                    doAlignScan = True, alignX = (-2.5,2.5,100,0.05,'Fe',0.5,True),
                    alignY = (-2.5,2.5,100,0.05,'Fe',0.5, True), xy_offset = (0,0),
                    pdfElem = ['Fe','As'],pdfLog = True,peakBeam = True,
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
    # marker to track beam dump             
    beamDumpOccured = False
                    
    e_list = generateEList_MLL(elemParam, highEStart =  highEStart)

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
    
    #remeber the start positions
    mot1_i = mot1.position
    mot2_i = mot2.position

    #calculate approximate time
    tot_time_ = (x_num*y_num*accq_t*len(e_list))
    tot_time = tot_time_/3600
    
    if doAlignScan:
        overhead = 1.5
    else:
        overhead = 1.25

    end_datetime = time.ctime(time.time()+tot_time_*overhead)



    for i in tqdm.tqdm(range(len(e_list)),desc = 'Energy Scan'):
    #for i in range (len(e_list)):

        #if beam dump occur turn the marker on
        if sclr2_ch2.get()<1000:
            beamDumpOccured = True
            cbpm_on(False)

        #wait if beam dump occured beamdump
        yield from check_for_beam_dump(threshold=5000)
        
        if beamDumpOccured:
            #wait for about 3 minutes for all the feedbacks to kick in
            yield from bps.sleep(120)

            #redo the previous energy
            e_t, hmll_hz_t,sbz_t, *others = e_list.iloc[i-1]

            #turn off the beamdump marker
            beamDumpOccured = False
            
        else:
            #unwrap df row for energy change
            e_t, hmll_hz_t,sbz_t, *others = e_list.iloc[i]
        

        yield from move_energy_mll(e_t,hmll_z = hmll_hz_t,sbz_ = sbz_t)

        
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

        elif doAlignScan:

            try:

                if alignX[-1]:
                    yield from fly1d(dets_fs,mot1,alignX[0],alignX[1],alignX[2],alignX[3])
                    xcen = return_line_center(-1,alignX[4],alignX[5])
                    #xcen,_ = erf_fit(-1,alignX[4])
                    yield from bps.mov(mot1, xcen)
                    print(f"mot1 centered to {xcen}")
                    plt.close()


                if alignY[-1]:
                    yield from fly1d(dets_fs,mot2,alignY[0],alignY[1],alignY[2],alignY[3])
                    ycen = return_line_center(-1,alignX[4],alignY[5])
                    #ycen,_ = erf_fit(-1,alignX[4])
                    yield from bps.mov(mot2, ycen)
                    print(f"mot2 centered to {ycen}")
                    plt.close()



            except:
                pass

        print(f'Current scan: {i+1}/{len(e_list)}')

        # do the fly2d scan
        #cbpm_on(False)

        if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.002) 
            #dead_time = 0.001 for 0.015 dwell

        else:

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
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
                insert_xrf_map_to_pdf(-1,pdfElem,title_=['energy', 'sbz', 'hz'])# plot data and add to pdf
            except:
                pass
        # save the DF in the loop so quitting a scan won't affect
        filename = f"HXN_nanoXANES_StartID{int(e_list['Scan ID'][0])}_{len(e_list)}_e_points.csv"
        e_list.to_csv(os.path.join(saveLogFolder, filename), float_format= '%.5f')
    '''
    #go back to max energy point if scans done reverese
    max_e_id = e_list['energy'].idxmax()
    e_max, hmll_max,v_mll_max *others = e_list.iloc[max_e_id]
    
    if not np.isclose(e_list['energy'].max(), e.position):
    
        yield from move_energy_mll(e_max,hmll_max,v_mll_max)
        
        yield from peak_the_flux()

    
    else: pass
    '''   
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
    if pdfLog: save_page() #save the pdf



def repeated_scan():
    yield from peak_the_flux()
    yield from mll_list_xanes2d_no_input(As_MLL_XANES,dets1,dssx,-0.9,0.7,40,dssy,-0.6,1,40,0.03, highEStart=False)
    yield from peak_the_flux()
    yield from mll_list_xanes2d_no_input(As_MLL_XANES,dets1,dssx,-0.9,0.7,40,dssy,-0.6,1,40,0.03, highEStart=True)
    yield from peak_the_flux()
    yield from mll_list_xanes2d_no_input(As_MLL_XANES,dets1,dssx,-0.9,0.7,40,dssy,-0.6,1,40,0.03, highEStart=False)
    yield from peak_the_flux()
    yield from mll_list_xanes2d_no_input(As_MLL_XANES,dets1,dssx,-0.9,0.7,40,dssy,-0.6,1,40,0.03, highEStart=True)






#<zp_list_xanes2d(LuL3XANES, dets_fs, zpssx,-15,15,100,zpssy, -15,15,100,0.02, highEStart=False,  alignX = (-10,10,10
   #...: 0,0.05,"Au_M",0.5, True), alignY = (-10,10,100,0.05,"Au_M",0.5,True), pdfElem=["Lu_L", "Au_M"], peakBeam=False,saveLogFolder="/GPFS/XF03ID1/users/2022
   #...: Q2/Tyson_2022Q2")
#<zp_list_xanes2d(LuL3XANES, dets_fs, zpssx,-15,15,100,zpssy, -15,15,100,0.02, highEStart=False,  alignX = (-10,10,100,0.05,"Au_M",0.5, True), alignY = (-10
   #...: ,10,100,0.05,"Au_M",0.5,True), pdfElem=["Lu_L", "Au_M"], peakBeam=False,saveLogFolder="/GPFS/XF03ID1/users/2022Q2/Tyson_2022Q2")

#<zp_list_xanes2d(FeXANES, dets1, zpssx,-2.0, 2.0, 80,zpssy,-3.5,0.5,80,0.030, highEStart=True)

def two_xanes():
    yield from mll_list_xanes2d_no_input(As_MLL_XANES,dets_fs,dssx,-2,2,125,dssy,-2,2,125,0.025, highEStart=True)
    yield from mll_list_xanes2d_no_input(As_MLL_XANES_minE,dets_fs,dssx,-2,2,200,dssy,-2,2,200,0.02, highEStart=False)
