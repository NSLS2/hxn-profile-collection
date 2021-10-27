#Last Update 09/11/2021 by AP

"""
ReadMe:

The workflow for xanes experiment is define below. This macro aims to use one flow for XANES of any given element. 
This macro is designed to work with the GUI inputs as well.
To add a new element add the paramer file in the format given below


EXAMPLE OF USAGE:

run -i XANES_Macros/hxn_nanoXANES.py 


For XANES Scan: <zp_list_xanes2d(ZnXANES,dets1,zpssx,-9,-1,80,zpssy,0,7,70,0.2,xcen = 0, ycen = 0,
alignX = (-2,2,100,0.1,'Fe',0.5),alignY = (-2,2,100,0.1, 'Fe',0.5), 
pdfElem = ['Fe','Zn'],saveLogFolder = '/data/users/2021Q3/Satish_2021Q3')


For Foil Calibration: <zp_list_xanes2d(ZnXANES,dets6,zpssx,-1,1,10,zpssy,-1,1,10,0.03,
highEStart = False, doAlignScan = False, pdfLog = False, foilCalibScan = True, peakBeam = False,
saveLogFolder ='/data/users/2021Q3/Satish_2021Q3') 
"""


import numpy as np
from datetime import datetime
import pandas as pd

#Paramer list from previous runs in the order of atomic number of the element

CrXANES = {'high_e':6.03, 'low_e':5.97, 
          'high_e_ugap':6660, 'low_e_ugap':6616,
          'high_e_crl':15, 'low_e_crl':15,'crl_comb':(8),
          'high_e_zpz1':71.395, 'zpz1_slope':-5.9,
          'energy':[(5.97,5.98,0.005),(5.981,6.03,0.001), (6.032,6.046,0.005)], 
          'mirrorCoating': 'Si'}
          
MnXANES = {'high_e':6.6, 'low_e':6.5, 
          'high_e_ugap':7142, 'low_e_ugap':7057,
          'high_e_crl':-12, 'low_e_crl':-12,'crl_comb':(8,6),
          'high_e_zpz1':68.3165, 'zpz1_slope':-5.9,
          'energy':[(6.520,6.530,0.005),(6.531,6.580,0.001),(6.585,6.601,0.005)],
          'mirrorCoating': 'Si'}
               
FeXANES = {'high_e':7.2, 'low_e':7.1, 
          'high_e_ugap':7678, 'low_e_ugap':7588,
          'high_e_crl':4, 'low_e_crl':-6,'crl_comb':(12),
          'high_e_zpz1':65.54, 'zpz1_slope':-5.9,
          'energy':[(7.085,7.11,0.005),(7.111,7.144,0.001)],
          'mirrorCoating': 'Si or Rh', 'zposaz':5000}

ZnXANES =  {'high_e':9.7, 'low_e':9.6, 
          'high_e_ugap':6485, 'low_e_ugap':6435,
          'high_e_crl':5., 'low_e_crl':2.,'crl_comb':(22+4),
          'high_e_zpz1':50.87, 'zpz1_slope':-5.9,
          'energy':[(9.64,9.666,0.005),(9.6665,9.681,.0005),(9.682,9.701,0.002),(9.705,9.725,0.005)]}


NiXANES = {'high_e':8.427, 'low_e':8.3, 
          'high_e_ugap':5834, 'low_e_ugap':5768,
          'high_e_crl':11.81, 'low_e_crl':8,'crl_comb':(12+4),
          'high_e_zpz1':58.32, 'zpz1_slope':-5.9,
          'energy':[(8.30,8.32,0.005),(8.321,8.38,0.001),(8.382,8.430,0.005)],
          'mirrorCoating': 'Cr','zposaz':3500}
          

    
                                ######################################
                                ######### FUNCTIONS BELOW ############
                                ######################################


def generateEPoints(ePointsGen = [(9.645,9.665,0.005),(9.666,9.7,0.0006),(9.705,9.725,0.005)],reversed = True):

    """ 
    
    Generates a list of energy values from the given list
    
    input: Tuples in the format (start energy, end energy, energy resolution), 
    if reversed is true the list will be transposed 
    
    return : list of energy points 
    
    """
                
    e_points = []
    
    for values in ePointsGen:
        #use np.arange to generate values and extend it to the e_points list
        e_points.extend(np.arange(values[0],values[1],values[2]))
    
    if reversed:
        #retruns list in the reversted order
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
    high_e, low_e = XANESParam['high_e'],XANESParam['low_e']
    high_e_ugap, low_e_ugap = XANESParam['high_e_ugap'],XANESParam['low_e_ugap']
    
    #slope = dUgap/dE
    ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
    
    #y = highvalues+dE*increment per eV
    ugap_list = high_e_ugap + (e_list['energy'] - high_e)*ugap_slope
    
    # add the list to DF
    e_list['ugap'] = ugap_list

    #same steps as ugap for CRL theta
    high_e_crl, low_e_crl =XANESParam['high_e_crl'],XANESParam['low_e_crl']
    crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
    crl_list = high_e_crl + (e_list['energy'] - high_e)*crl_slope
    e_list['crl_theta'] = crl_list
    
    #zone plate increament is very close to the theorticla value , same step as above for zp focus
    zpz1_ref, zpz1_slope = XANESParam['high_e_zpz1'],XANESParam['zpz1_slope']
    zpz1_list = zpz1_ref + (e_list['energy'] - high_e)*zpz1_slope
    e_list['ZP focus'] = zpz1_list
    
    #return the dataframe
    return e_list 

def peak_the_flux():

    """ Scan the c-bpm set points to find IC3 maximum """

    print("IC3is below threshold; Peaking the beam.")
    yield from bps.sleep(2)
    yield from peak_bpm_y(-5,5,10)
    yield from bps.sleep(1)
    yield from peak_bpm_x(-15,15,6)
    yield from bps.sleep(1)
    yield from peak_bpm_y(-2,2,4)
    
def move_energy(e_,ugap_,zpz_,crl_th_, ignoreCRL= False, ignoreZPZ = False):
    
    
    """ Function to change energy knowing ugap, 
    crl_theta (optional) and zone plate focus (optional)"""

    print(f'{15*"*"} \n, Target Energy = {e_:.4f} \n,UGap = {ugap_:.4f} \n,CRL Theta = {crl_th_:.4f} \n,ZPZ = {zpz_:4f}, \n {15*"*"} ')
            
    #tuning the scanning pv on to dispable c bpms
    caput('XF:03IDC-ES{Status}ScanRunning-I', 1)  

    if abs(e.position-e_)<0.2:

        #move mono e
        yield from bps.mov(e,e_)
        yield from bps.sleep(1)
        
        #move ugap
        yield from bps.mov(ugap, ugap_)
        yield from bps.sleep(2)
        
        #move crl_theta, if true
        if not ignoreZPZ: yield from mov_zpz1(zpz_)
        yield from bps.sleep(1)
        
        #move zpz1, if true
        if not ignoreCRL: yield from bps.mov(crl.p,crl_th_)
        yield from bps.sleep(1)
        
        #scanning pv off to dispable c bpms
        caput('XF:03IDC-ES{Status}ScanRunning-I', 0)
        yield from bps.sleep(5)

    else: 
        raise ValueError("The requested energy change is more than 200 eV. This must be a mistake")
    

def zp_list_xanes2d(elemParam,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = True,
                    doAlignScan = True, xcen = 0, ycen = 0,
                    alignX = (-2,2,100,0.1,'Fe',0.7, True),alignY = (-2,2,100,0.1,'Fe',0.7, True), 
                    pdfElem = ('Fe','Cr'),doScan = True, moveOptics = True,  
                    pdfLog = True, foilCalibScan = False, peakBeam = True,
                    saveLogFolder = '/home/xf03id/Downloads'):
                    
                    
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
           9. 
           3. Options for reginstration scans
           4. Options to save XRFs to pdf after each scan
           5. Options to do foil calibration scans
           6. Save important information in CSV format to selected forlder 
           7. The user can turn on and off alignemnt scans
    
    
    """   
    # marker to track beam dump             
    beamDumpOccured = False
                    

    e_list = generateEList(XANESParam = elemParam, highEStart =  highEStart)

    #add real energy to the dataframe
    e_list['E Readback'] = np.nan 
    
    #add scan id to the dataframe
    e_list['Scan ID'] = np.nan 
    
    #recoed time
    e_list['TimeStamp'] = pd.Timestamp.now()
    
    #Ic values are useful for calibration
    e_list['IC3'] = sclr2_ch4.get() 
    e_list['IC0'] = sclr2_ch2.get()
    e_list['IC3_before_peak'] = sclr2_ch2.get()
    
    
    #record if peak beam happed before the scan   
    e_list['Peak Flux'] = False 
    
    print(e_list.head())
    yield from bps.sleep(10)#time to quit if anything wrong
    
    #get intal ic1 value
    ic_0 = sclr2_ch2.get()
    
    #opening fast shutter for initial ic3 reading
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
    yield from bps.sleep(2)
    
    #get the initial ic3 reading for peaking the beam
    ic_3_init =  sclr2_ch4.get()
     
    #close fast shutter after initial ic3 reading
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
    
    #remeber the start positions
    zpssx_i = zpssx.position
    zpssy_i = zpssy.position


    for i in range (len(e_list)):

        #if beam dump occur turn the marker on
        if sclr2_ch2.get()<0.1*ic_0:
            beamDumpOccured = True

        #wait if beam dump occured beamdump
        yield from check_for_beam_dump(threshold=0.1*ic_0)
        
        if beamDumpOccured:
            #wait for about 3 minutes for all the feedbacks to kick in
            yield from bps.sleep(200)

            #redo the previous energy
            e_t, ugap_t, crl_t, zpz_t, *others = e_list.iloc[i-1]

            #turn off the beamdump marker
            beamDumpOccured = False
            
        else:
            #unwrap df row for energy change
            e_t, ugap_t, crl_t, zpz_t, *others = e_list.iloc[i]
        
        if moveOptics: 
            yield from move_energy(e_t,ugap_t,zpz_t,crl_t,
                                   ignoreCRL= foilCalibScan, 
                                   ignoreZPZ = foilCalibScan)

        else: pass
        
        #open fast shutter to check if ic3 reading is satistactory
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
        yield from bps.sleep(3)
        
        #get ic3 value before peaking, e change
        ic3_ = sclr2_ch4.get()
        
        # if ic3 value is below the threshold, peak the beam
        if ic3_ < ic_3_init*0.9:
            
            if peakBeam: yield from peak_the_flux()
            fluxPeaked = True # for df record
        else:
            fluxPeaked = False
        
        #for df
        ic_3 = sclr2_ch4.get()
        ic_0 = sclr2_ch2.get()

        # move to particle location for alignemnt scan
        if doAlignScan:
        
            yield from bps.mov(zpssx, xcen)
            yield from bps.mov(zpssy, ycen)
        
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
            if alignX[-1]:
                yield from fly1d(dets,zpssx,alignX[0],alignX[1],alignX[2],alignX[3])
                xcen = return_line_center(-1,alignX[4],alignX[5])
                yield from bps.mov(zpssx, xcen)

            if alignY[-1]:
                yield from fly1d(dets,zpssy,alignY[0],alignY[1],alignY[2],alignY[3])
                ycen = return_line_center(-1,alignX[4],alignY[5])
                yield from bps.mov(zpssy, ycen)


        print(f'Current scan: {i+1}/{len(e_list)}')

        # do the fly2d scan

        if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

            if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.002) 
            #dead_time = 0.001 for 0.015 dwell

        else:

            if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(1)

        #close fast shutter
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        
        # after scan done go to 0,0 to rest
        if doAlignScan: 
            yield from bps.mov(zpssx, zpssx_i)
            yield from bps.mov(zpssy, zpssy_i)

        #ycen, xcen = return_center_of_mass_blurr(-1,'S') 
        # some cases use 2D mass center for alignemnt
        #print(ycen,xcen)

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
            for elem in pdfElem:
                insert_xrf_map_to_pdf(-1,elem)# plot data and add to pdf

        # save the DF in the loop so quitting a scan won't affect
        filename = f"HXN_nanoXANES_StartID{int(e_list['Scan ID'][0])}_{len(e_list)}_e_points.csv"
        e_list.to_csv(os.path.join(saveLogFolder, filename), float_format= '%.5f')

    #go back to max energy point if scans done reverese
    max_e_id = e_list['energy'].idxmax()
    e_max, ugap_max,  crl_max,zpz_max, *others = e_list.iloc[max_e_id]
    
    if not np.isclose(e_list['energy'].max(), e.position):
    
        yield from move_energy(e_max,ugap_max,zpz_max,crl_max,
                               ignoreCRL= foilCalibScan,
                               ignoreZPZ = foilCalibScan)
        
        yield from peak_the_flux()

    
    else: pass
        
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
    if pdfLog: save_page() #save the pdf



