import numpy as np
import pandas as pd
import time
from datetime import datetime
import scipy.constants as consts


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

    return np.around(np.sqrt(energy*ETOK),2)

def ktoe(k):
    """convert photo-electron wavenumber to energy"""
    if isinstance(k, list):
        k = np.array(k)
    return np.around(k*k*KTOE, 1)


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
    yield from peak_bpm_y(-4,4,10)
    yield from bps.sleep(1)
    yield from peak_bpm_x(-10,10,6)
    yield from bps.sleep(1)
    yield from peak_bpm_y(-2,2,4)
    

def move_energy(e,zpz_ ):
    
    cbpm_on(False)
    yield from bps.sleep(1)

    #tuning the scanning pv on to dispable c bpms
    #caput('XF:03IDC-ES{Status}ScanRunning-I', 1)


    yield from Energy.move(e, moveMonoPitch=False, moveMirror = "ignore")
    yield from mov_zpz1(zpz_)
    yield from bps.sleep(4)

    cbpm_on(True)


class EXAFS():

    def __init__(self, param_file, edge = "K"):

        self.param_file = param_file
        self.edge = "K"

        self.edge_selection = {"K":xraylib.K_SHELL}

    @staticmethod
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

    
    def generateEList(self, highEStart = False, startFrom = 0):

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
        elem = self.param_file['elem']
        kmin = self.param_file['kmin']
        kmax = self.param_file['kmax']
        kstep = self.param_file['kstep']

        e_list['energy'] = generateEPoints(self.param_file['pre_edge'],elem,kmin, kmax,kstep, reversed = highEStart)

        #read the paramer dictionary and calculate ugap list
        high_e, low_e = self.param_file['high_e'],self.param_file['low_e']

        #zone plate increament is very close to the theorticla value , same step as above for zp focus
        zpz1_ref, zpz1_slope = self.param_file['high_e_zpz1'],self.param_file['zpz1_slope']
        zpz1_list = zpz1_ref + (e_list['energy'] - high_e)*zpz1_slope
        e_list['ZP focus'] = zpz1_list

        #return the dataframe
        return e_list[startFrom:]


    def run(self,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,highEStart = False,
            doAlignScan = True, alignX = (-2,2,100,0.1,'Fe',0.7, True),
            alignY = (-2,2,100,0.1,'Fe',0.7, True), 
            pdfElem = ('Fe','Cr'),doScan = True, moveOptics = True,pdfLog = True, 
            foilCalibScan = False, peakBeam = False, startEPoint = 0,
            saveLogFolder = '/home/xf03id/Downloads')


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
                        
        e_list = generateEList(self.param_file, highEStart =  False, startFrom = startEPoint)

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
        yield from bps.sleep(1)#time to quit if anything wrong
        
        #get intal ic1 value
        ic_0 = sclr2_ch2.get()
        
        #opening fast shutter for initial ic3 reading
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) 
        yield from bps.sleep(5)
        
        #get the initial ic3 reading for peaking the beam
        ic_3_init =  sclr2_ch4.get()
        
        #close fast shutter after initial ic3 reading
        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        
        #remeber the start positions
        mot1_i = mot1.position
        mot2_i = mot2.position
        tot_time_s = (x_num*y_num*accq_t*len(e_list))
        tot_time = tot_time_s/3600

        proj_end_time_s = time.time()+tot_time_s

        check = input(f"This plan takes about {tot_time*2 :.1f} hours \n", 
                      f"Projected end time is {time.strftime('%Y-%m-%d, %A, %H:%M:%S', time.localtime(proj_end_time_s*2))}\n"
                      "continue (y/n)?")

        if check == "y":

            for i in range (len(e_list)):

                #if beam dump occur turn the marker on
                if sclr2_ch2.get()<10000:
                    beamDumpOccured = True
                    cbpm_on(False)

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
                cbpm_on(False)

                if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

                    if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t) 
                    #dead_time = 0.001 for 0.015 dwell

                else:

                    if doScan: yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
                yield from bps.sleep(1)

                cbpm_on(True)

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
            
            yield from peak_the_flux()

        
        else: pass
        
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) 
        if pdfLog: save_page() #save the pdf

        else:
            return



class XANES():

        pass

fe_exafs = EXAFS(FeEXAFS)