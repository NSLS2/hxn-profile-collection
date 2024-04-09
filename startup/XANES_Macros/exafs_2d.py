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
from datetime import datetime
import pandas as pd
import xraylib, tqdm
import scipy.constants as consts


ZnEXAFS= {'high_e':9.7, 
          'low_e':9.6,
          'high_e_zpz1':50.92, 
          'zpz1_slope':-5.9,
          'energy':[]}

CrEXAFS= {'high_e':6.03, 
          'low_e':5.97,
          'high_e_zpz1':71.395, 
          'zpz1_slope':-5.9,
          'energy':[(5.97,5.98,0.005),(5.981,6.03,0.001), (6.032,6.046,0.005)],
          }

MnEXAFS= {'high_e':6.6, 
          'low_e':6.5,
          'high_e_zpz1':68.31, 
          'zpz1_slope':-5.9,
          'energy':[(6.520,6.530,0.005),(6.531,6.580,0.001),(6.585,6.601,0.005)],
          'mirror': 'Si'}

FeEXAFS= {'high_e':7.6, 
          'low_e':7.1,
          'high_e_zpz1':4.535, 
          'zpz1_slope':-5.04,
          'pre_edge':[(6.97,7.11,0.010)],
          'elem':"Fe", 
          'kmin':2, 
          'kmax':12, 
          'kstep':0.1}

LuL3EXAFS =  {'high_e':9.3, 
              'high_e_zpz1':-5.425, 
              'zpz1_slope':-5.04,
              'low_e':9.2,
            'pre_edge':[(9.150,9.200,0.005)],
            'elem':"Lu", 
            'kmin':2, 
            'kmax':12, 
            'kstep':0.1}

                                ######################################
                                ######### FUNCTIONS BELOW ############
                                ######################################
#copied from larch --modified

KTOE = 1.e20*consts.hbar**2 / (2*consts.m_e * consts.e) # 3.8099819442818976
ETOK = 1.0/KTOE

def etok(energy):
    """convert photo-electron energy to wavenumber"""
    if isinstance(energy, list):
        energy = np.array(energy)

    if energy < 0: return 0

    return np.around(np.sqrt(energy*ETOK),5)

def ktoe(k):
    """convert photo-electron wavenumber to energy"""
    if isinstance(k, list):
        k = np.array(k)
    return np.around(k*k*KTOE, 1)


def generateEPoints(ePointsGen,elem,kmin, kmax,kstep, reversed = True):

    """

    Generates a list of energy values from the given list

    input: Tuples in the format (start energy, end energy, energy resolution),
    example: [(9.645,9.665,0.005),(9.666,9.7,0.0006),(9.705,9.725,0.005)]
    
    if reversed is true the list will be transposed

    return : list of energy points

    """
    E0 = xraylib.EdgeEnergy(xraylib.SymbolToAtomicNumber(elem), xraylib.K_SHELL) + ktoe(kmin)*0.001 # kstart=2
    print(f"{E0 =}")

    krangeE = E0 + 0.001*ktoe(np.arange(0.4,kmax,kstep))
    print(f"{krangeE =}")

    e_points = []


    if isinstance(ePointsGen[0], tuple):
    
        for values in ePointsGen:
            #use np.arange to generate values and extend it to the e_points list
            e_points.extend(np.arange(values[0],values[1],values[2]))
            edgeStart = values[1]+0.001

        e_points.extend(np.arange(edgeStart,E0,0.001))
        e_points.extend(krangeE)


    elif isinstance(ePointsGen, list):
        e_points = ePointsGen

    else:
        print (" Unknown energy list format")

    if reversed:
        #retrun list in the reversted order
        return np.around(e_points[::-1],5)
    else:
        return np.around(e_points,5)


def generateEList(EXAFSParam = FeEXAFS, highEStart = False, startFrom = 0):

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
    elem = EXAFSParam['elem']
    kmin = EXAFSParam['kmin']
    kmax = EXAFSParam['kmax']
    kstep = EXAFSParam['kstep']

    e_list['energy'] = generateEPoints(EXAFSParam['pre_edge'],elem,kmin, kmax,kstep, reversed = highEStart)

    #read the paramer dictionary and calculate ugap list
    high_e, low_e = EXAFSParam['high_e'],EXAFSParam['low_e']

    #zone plate increament is very close to the theorticla value , same step as above for zp focus
    zpz1_ref, zpz1_slope = EXAFSParam['high_e_zpz1'],EXAFSParam['zpz1_slope']
    zpz1_list = zpz1_ref + (e_list['energy'] - high_e)*zpz1_slope
    e_list['ZP focus'] = zpz1_list

    #return the dataframe
    return e_list[startFrom:]


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


def peak_the_flux():

    cbpm_on(True)

    """ Scan the c-bpm set points to find IC3 maximum """

    print("IC3is below threshold; Peaking the beam.")
    yield from bps.sleep(2)
    yield from peak_bpm_y(-5,5,10)
    yield from bps.sleep(1)
    yield from peak_bpm_x(-10,10,10)
    yield from bps.sleep(1)
    yield from peak_bpm_y(-3,3,10)

def move_energy(e,zpz_ ):
    
    yield from bps.sleep(1)
    #tuning the scanning pv on to dispable c bpms
    #caput('XF:03IDC-ES{Status}ScanRunning-I', 1)
    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")
    yield from mov_zpz1(zpz_)
    yield from bps.sleep(4)

    
def zp_list_exafs2d(elemParam,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = True,
                    doAlignScan = True, alignX = (-2,2,100,0.1,'Fe',0.7, True),
                    alignY = (-2,2,100,0.1,'Fe',0.7, True), 
                    pdfElem = ('Fe','Cr'),doScan = True, moveOptics = True,pdfLog = True, 
                    foilCalibScan = False, peakBeam = False, startEPoint = 0,
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
           9. xcen, ycen; positions where alignemnt scan would be done. This number updates after each alignment scan
           10. Options for reginstration scans
           11. Options to save XRFs to pdf after each scan
           12. Options to do foil calibration scans
           13. Save important information in CSV format to selected forlder 
           14. The user can turn on and off alignemnt scans
    
    
    """   
    # marker to track beam dump             
    beamDumpOccured = False
                    
    e_list = generateEList(elemParam, highEStart =  highEStart, startFrom = startEPoint)

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
    yield from bps.sleep(3)
    
    #get the initial ic3 reading for peaking the beam
    ic_3_init =  sclr2_ch4.get()
     
    #close fast shutter after initial ic3 reading
    #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
    
    #remeber the start positions
    mot1_i = mot1.position
    mot2_i = mot2.position

    
    tot_time = (x_num*y_num*accq_t*len(e_list))/3600

    check = input(f"This plan takes about {tot_time*1.5 :.1f} hours \n"
                    f"Projected end time is {time.strftime('%Y-%m-%d, %A, %H:%M:%S', time.localtime(time.time()+tot_time*3600*1.5))}\n"
                    "continue (y/n)?")

    if check == "y":

        for i in tqdm.tqdm(range(len(e_list)),desc = 'Energy Scan'):

            #if beam dump occur turn the marker on
            if sclr2_ch2.get()<50000:
                beamDumpOccured = True
                

            #wait if beam dump occured beamdump
            yield from check_for_beam_dump(threshold=10000)
            
            if beamDumpOccured:
                #wait for about 3 minutes for all the feedbacks to kick in
                yield from bps.sleep(200)

                #redo the previous energy
                e_t, zpz_t, *others = e_list.iloc[i-1]

                #turn off the beamdump marker
                beamDumpOccured = False
                
            else:
                #unwrap df row for energy change
                e_t, zpz_t, *others = e_list.iloc[i]
            
            if moveOptics: 
                yield from move_energy(e_t,zpz_t)

            else: pass
            
            #open fast shutter to check if ic3 reading is satistactory
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
            yield from bps.sleep(5)
            
            #get ic3 value before peaking, e change
            ic3_ = sclr2_ch4.get()
            
            # if ic3 value is below the threshold, peak the beam
            if ic3_ < ic_3_init*0.9:
                
                #if peakBeam: yield from peak_the_flux()
                if peakBeam: yield from peak_xy_volt(2)
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
                if alignX[-1]:
                    yield from fly1d(dets_fs,zpssx,alignX[0],alignX[1],alignX[2],alignX[3])
                    xcen = return_line_center(-1,alignX[4],alignX[5])
                    yield from bps.movr(smarx, xcen*0.001)
                    print(f"zpssx centered to {xcen}")

                if alignY[-1]:
                    yield from fly1d(dets_fs,zpssy,alignY[0],alignY[1],alignY[2],alignY[3])
                    ycen = return_line_center(-1,alignX[4],alignY[5])
                    yield from bps.movr(smary, ycen*0.001)
                    print(f"zpssy centered to {ycen}")


            print(f'Current scan: {i+1}/{len(e_list)}')

            # do the fly2d scan
            

            if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

                if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t) 
                #dead_time = 0.001 for 0.015 dwell

            else:

                if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
            yield from bps.sleep(1)

            

            #close fast shutter
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
            
            # after scan done go to 0,0 to rest
            #if doAlignScan: 
                #yield from bps.mov(zpssx, zpssx_i)
                #yield from bps.mov(zpssy, zpssy_i)

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
            
            yield from peak_xy_volt(2)

        
        else: pass
            
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if pdfLog: save_page() #save the pdf

    else:
        return

def peak_b_bpm(bpm_name, start, end, n_steps):
    shutter_b_cls_status = caget('XF:03IDB-PPS{PSh}Sts:Cls-Sts')
    shutter_c_status = caget('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0')


    if shutter_b_cls_status == 0:

        caput('XF:03IDC-ES{Status}ScanRunning-I', 1)
        bpm_0 = caget(bpm_name)
        x = np.linspace(bpm_0+start,bpm_0+end,n_steps+1)
        y = np.arange(n_steps+1)
        #print(x)
        for i in range(n_steps+1):
            caput(bpm_name,x[i])
            if i == 0:
                yield from bps.sleep(5)
            else:
                yield from bps.sleep(2)

            if shutter_c_status == 0:
                y[i] = sclr2_ch2.get()

            else:
                y[i] = sclr2_ch4.get()


        peak = x[y == np.max(y)]
        caput(bpm_name,peak[0])
        yield from bps.sleep(2)

    else:
        print('Shutter B is Closed')


'''

    #plt.pause(5)
    #plt.close()


    #zp_list_xanes2d(FeXANES, dets1, zpssx,-4,4,80,zpssy, -4,4,80,0.05, highEStart=False, alignElem='S', alignX = (-5,5,100,0.05,0.5), alignY = (-5,5,100,0
   #...: .05,05), pdfElem=["Fe", "S"], saveLogFolder="\data\Staff\Ajith\2022Q2")"


   #zp_list_exafs2d(LuL3EXAFS,dets_fs,zpssx,-15,15,100,zpssy, -15,15,100,0.03,,\
   highEStart = False,alignX = (-10,10,100,0.05,'Au_M',0.7, True),
   alignY = (-10,10,100,0.05,'Au_M',0.5, True), pdfElem = ('Lu_L','Au_M'), 
   pdfLog = False, peakBeam = False,saveLogFolder = "/GPFS/XF03ID1/users/2022Q2/Tyson_2022Q2")



   zp_list_exafs2d(FeEXAFS,dets1,zpssx,-2.5,2.5,100,zpssy,-2.5,2.5,100,0.05,highEStart = False,
                    doAlignScan = True, alignX = (-4,4,100,0.05,'Cr',0.7, True),
                    alignY = (-4,4,100,0.05,'Cr',0.7, True), 
                    pdfElem = ['Fe','Cr'],pdfLog = True, 
                    , peakBeam = True, startEPoint = 0,
                    saveLogFolder = '/nsls2/data/hxn/legacy/users/2022Q3/Ajith_2022Q3/nano-EXAFS/')

    <zp_list_exafs2d(FeEXAFS,
                     dets1,
                     zpssx,
                     -2,
                     2,
                     80,
                     zpssy,
                     -2,
                     2,
                     80,
                     0.05,
                     highEStart = False, 
                     doAlignScan = True, 
                     alignX = (-4,4,100,0.05,'Cr',0.7, True), 
                     alignY = (-4,4,100,0.05,'Cr',0.7, True),
                     pdfElem = ['Fe','Cr'],
                     pdfLog = True, 
                     peakBeam = True, 
                     startEPoint = 0,
                     saveLogFolder = '/nsls2/data/hxn/legacy/users/2022Q3/Ajith_2022Q3/nano-EXAFS/')



    <zp_list_exafs2d(FeEXAFS,dets1,zpssx,-2,2,80,zpssy,-2,2,80,0.05,highEStart = False, doAlignScan = True, align
    ...: X = (-4,4,100,0.05,'Cr',0.7, True), alignY = (-4,4,100,0.05,'Cr',0.7, True),   pdfElem = ['Fe','Cr'],pdfLog = True, pea
    ...: kBeam = True, startEPoint = 0,saveLogFolder = '/nsls2/data/hxn/legacy/users/2022Q3/Ajith_2022Q3/nano-EXAFS/')


   '''

def run_exafs():

    yield from zp_list_exafs2d(FeEXAFS,
                            dets_fs,
                            zpssx,
                            -2.5,
                            1.8,
                            43,
                            zpssy,
                            -1.8,
                            2.5,
                            43,
                            0.1,
                            highEStart = True, 
                            doAlignScan = True, 
                            alignX = (-5,5,100,0.03,'Cr',0.5, True), 
                            alignY = (-5,5,100,0.03,'Cr',0.2, True),
                            pdfElem = ['Fe','Cr'],
                            pdfLog = True, 
                            peakBeam = True, 
                            startEPoint = 0,
                            saveLogFolder = '/nsls2/data/hxn/legacy/users/2023Q1/Ajith_2023Q1')