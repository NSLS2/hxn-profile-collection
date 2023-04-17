print(f"Loading {__file__!r} ...")

import logging
import numpy as np
import pandas as pd
import pyqtgraph as pg
import matplotlib.pyplot as plt
from epics import caget, caput

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

wd = "/nsls2/data/hxn/shared/config/bluesky/profile_collection/startup/"

class HXNEnergy():
    
    def __init__(self, ugap_,bragg_e, dcm_pitch, ic1, calib_file_csv):

        self.ugap = ugap_
        self.bragg_e = bragg_e
        self.dcm_pitch = dcm_pitch
        self.ic1 = ic1
        self.calib_file = calib_file_csv
        self.df = pd.read_csv(self.calib_file)
        self.calibrate_optics(False)
        #self.calibrate_ugap()
        #self.calibrate_dcm_pitch()
        #self.calibrate_hfm_pitch()
        #self.calibrate_dcm_roll()

        self.defineFeedBackPVs()

    def defineFeedBackPVs(self):
        self.hcm_pf_sts = "XF:03IDC-CT{FbPid:01}PID:on"
        self.hfm_pf_sts = "XF:03IDC-CT{FbPid:02}PID:on"

        self.xbpm_x_rbv = "XF:03ID-BI{EM:BPM1}PosX:MeanValue_RBV"
        self.xbpm_y_rbv = "XF:03ID-BI{EM:BPM1}PosY:MeanValue_RBV"

        self.xbpm_x_val = "XF:03ID-BI{EM:BPM1}fast_pidX.VAL"
        self.xbpm_x_val = "XF:03ID-BI{EM:BPM1}fast_pidY.VAL"


    def get_values(self,targetE):

        u_gap, harmonics = self.gap(targetE, harmonicsChoice = -1)
        dcm_pitch = self.calculatePitch(targetE)
        dcm_roll = self.calculateRoll(targetE)
        hfm_pitch = self.calculateHFMPitch(targetE)
        mirror_coating = self.findAMirror(targetE)

        logger.info(f"{u_gap = :.1f}\n{harmonics = }\n{dcm_pitch = :.4f}\n"
                    f"{dcm_roll = :.4f}\n{hfm_pitch = :.4f}\n{mirror_coating = }")

        return u_gap, harmonics,dcm_pitch,dcm_roll,hfm_pitch,mirror_coating



    def calcGap(self,E,harmonics = 5, offset = 0):
        E1 = E/harmonics
        calc_gap =  np.polyval(self.ugap_coeffs, E1) + offset
        return (np.around(calc_gap,1))

    def gap(self,targetE, harmonicsChoice = -1):
        
        "Estimate a suitable gap for target energy"

        if harmonicsChoice == -1:
        
            harmonics = np.array([3,5,7,9]) #harmonics options

            if targetE>4.99 and targetE<25:     # if dialed a wrong number
                
                opt = np.array([self.calcGap(targetE, harmonics = hm) for hm in harmonics])
                idx = np.where(np.logical_and(opt>=6000, opt<=10000)) #6000 beacuse of the ugap scan limit
                gap =  opt[idx][-1] #last one has lowest gap in the list
                logger.info(f" estimated gap = {gap}")
                

                harmonic = harmonics[idx][-1]
                logger.info(f" Harmonic = {harmonic}")

                return gap, harmonic
            
            else:
                raise ValueError(" Requested Energy is out of range")
                return

        else:

            gap = self.calcGap(targetE, harmonics = int(harmonicsChoice))

            return gap, int(harmonicsChoice)
    
    
    def findAMirror(self,e):
    
        EnergyToMirror = {(5,8.5):'Si',(8.5,15):'Cr',(15,20):'Rh',(20,25):'Pt'}
    
        if not 4.99 <= e <= 25:
            raise ValueError (" Energy value out of range")
        
        else:
            for erange in EnergyToMirror.keys():
                if erange[0] <= e <= erange[1]:
                    return EnergyToMirror[erange]
                else:
                    pass
    

    def moveMirror(self,targetE, mirror = "auto"):
        
         MirrorPos = {'Si':(21,-4),'Cr':(5,15),'Rh':(30,-12),'Pt':(13.5,4)}         
         
         if mirror == "auto":
             foundMirror = self.findAMirror(targetE)
             positions = MirrorPos[foundMirror]
             logger.info(f"Moving to {foundMirror}")    
         else:
             positions = MirrorPos[mirror]
             logger.info(f"Moving to {mirror}")

         caput("XF:03IDA-OP{Mir:1-Ax:Y}Mtr",positions[0])
         #caput("XF:03IDA-OP{Mir:2-Ax:Y}Mtr",positions[1])
         yield from bps.mov(m1.y, positions[0], m2.y, positions[1] )

    def calibrate_optics(self, plot_after = False):

        en = self.df["energy"].to_numpy()

        adj_E = en/self.df["harmonic"].to_numpy()
        ugaps = self.df["ugap"].to_numpy()
        self.ugap_coeffs = np.polyfit(adj_E, ugaps, 3)

        dcm_p = self.df["dcmPitch"].to_numpy()
        self.dcm_p_coeffs = np.polyfit(en, dcm_p, 3)

        m2_p = self.df["hfmPitch"].to_numpy()
        self.m2p_coeffs = np.polyfit(en, m2_p, 3)

        dcm_r = self.df["dcmRoll"].to_numpy()
        self.dcm_r_coeffs = np.polyfit(en, dcm_r, 3)

        if plot_after:

            fig, axs = plt.subplots(2,2,figsize = (8,8))
            fig.subplots_adjust(hspace = 0.8, wspace = 0.8)
            fig.suptitle("Energy Calibration")
            fig.show()
            axs = axs.ravel()

            axs[0].scatter(adj_E,ugaps,label='ugap')
            axs[0].plot(adj_E, np.polyval(self.ugap_coeffs, adj_E),'r', label='ugap_fit')
            axs[0].set_title("First order E vs Ugap")
            axs[0].set_xlabel("First Order Energy (keV)")
            axs[0].set_ylabel("undulator Gap (um)")

            axs[1].scatter(en,dcm_p,label='dcm_p')
            axs[1].plot(en, np.polyval(self.dcm_p_coeffs, en),'r', label='dcm_p_fit')
            axs[1].set_title("E vs dcm_pitch")
            axs[1].set_xlabel("Energy (keV)")
            axs[1].set_ylabel("dcm_pitch")

            axs[2].scatter(en,m2_p,label='m2_p')
            axs[2].plot(en, np.polyval(self.m2p_coeffs, en),'r', label='m2_p_fit')
            axs[2].set_title("E vs m2_pitch")
            axs[2].set_xlabel("Energy (keV)")
            axs[2].set_ylabel("m2_pitch")
            
            axs[3].scatter(en,dcm_r,label='dcm_r')
            axs[3].plot(en, np.polyval(self.dcm_r_coeffs, en),'r', label='dcm_r_fit')
            axs[3].set_title("E vs dcm_roll")
            axs[3].set_xlabel("Energy (keV)")
            axs[3].set_ylabel("dcm_roll")
        


    '''
    def calibrate_ugap(self, plot_after = False):

        adj_E = self.df["energy"].to_numpy()/self.df["harmonic"].to_numpy()
        self.ugap_coeffs = np.polyfit(adj_E, self.df["ugap"].to_numpy(), 3)


        if plot_after:

            plt.figure()
            plt.scatter(adj_E,self.df["ugap"],label='data')
            plt.plot(adj_E, np.polyval(self.ugap_coeffs, adj_E),'r', label='fit')
            plt.xlabel("First Order Energy (keV)")
            plt.ylabel("undulator Gap (um)")
            plt.legend()
            plt.title(f"Undulator Calib_{pd.Timestamp.now().month}_{pd.Timestamp.now().year}")
            plt.show()

    def calibrate_dcm_pitch(self,plot_after = False):
        
        en, dcm_p = self.df["energy"].to_numpy(),self.df["dcmPitch"].to_numpy()
        self.dcm_p_coeffs = np.polyfit(en, dcm_p, 3)
        
        if plot_after:

            plt.figure()
            plt.scatter(self.df["energy"],self.df["dcmPitch"],label='dcm_pitch')
            plt.plot(self.df["energy"], np.polyval(self.dcm_p_coeffs, self.df["energy"]),'r', label='fit')
            plt.xlabel("Energy (keV)")
            plt.ylabel("undulator Gap (um)")
            plt.legend()
            plt.title(f"Undulator Calib_{pd.Timestamp.now().month}_{pd.Timestamp.now().year}")
            plt.show()

    def calibrate_hfm_pitch(self,plot_after = False):

        self.m2p_coeffs = np.polyfit(self.df["energy"].to_numpy(), self.df["hfmPitch"].to_numpy(), 3)


        if plot_after:

            plt.figure()
            plt.scatter(self.df["energy"],self.df["hfmPitch"],label='hfm_pitch')
            plt.plot(self.df["energy"], np.polyval(self.m2p_coeffs, self.df["energy"]),'r', label='fit')
            plt.xlabel("Energy (keV)")
            plt.ylabel("hfm_pitch (um)")
            plt.legend()
            plt.title(f"Undulator Calib_{pd.Timestamp.now().month}_{pd.Timestamp.now().year}")
            plt.show()

    def calibrate_dcm_roll(self, plot_after = False):
        
        self.dcm_r_coeffs = np.polyfit(self.df["energy"].to_numpy(), self.df["dcmRoll"].to_numpy(), 3)

        if plot_after:

            plt.figure()
            plt.scatter(self.df["energy"],self.df["dcmRoll"],label='dcm_roll')
            plt.plot(self.df["energy"], np.polyval(self.dcm_r_coeffs, self.df["energy"]),'r', label='fit')
            plt.xlabel("Energy (keV)")
            plt.ylabel("dcm_roll (um)")
            plt.legend()
            plt.title(f"Undulator Calib_{pd.Timestamp.now().month}_{pd.Timestamp.now().year}")
            plt.show()

    '''

    def calculatePitch(self,targetE, offset = 0):
        
        calc_pitch =  np.polyval(self.dcm_p_coeffs, targetE) + offset
        return np.around(calc_pitch,4) 

    def calculateRoll(self,targetE):
        calc_r =  np.polyval(self.dcm_r_coeffs, targetE)
        return (np.around(calc_r,4))

    def calculateHFMPitch(self,targetE):
        calc_m2p =  np.polyval(self.m2p_coeffs, targetE)
        return (np.around(calc_m2p,4))

    def move(self,targetE, harmChoice = -1, moveMonoPitch = True, moveMirror = "auto") :
        
        bbpm_auto = "XF:03ID{XBPM:17}AutoFbEn-Cmd"
        bbpm_x = "XF:03ID-BI{EM:BPM1}fast_pidX.FBON"
        bbpm_y = "XF:03ID-BI{EM:BPM1}fast_pidY.FBON"

        gap, hrm = Energy.gap(targetE,harmChoice)
        dcm_p_target = self.calculatePitch(targetE)
        hfm_p_target = self.calculateHFMPitch(targetE)
        dcm_r_target = self.calculateRoll(targetE)


        logger.info(f"{dcm_p_target = :.4f}, {hfm_p_target = :.2f}, {dcm_r_target = :.4f}")

        if gap<5200 and gap>10000:
            raise ValueError ("Incorrect gap calculation")
        else:
            #logger.info(f"Moving gap = {gap}")
            yield from bps.mov(ugap, gap)
            logger.info("Gap moved")
        
            logger.info(f"Mono Energy Target = {targetE}")
            yield from bps.mov(e,targetE)
            logger.info("Energy reached")

            if moveMonoPitch:

                    
                logger.info(f"Moving {dcm_p_target = :4f}")
                yield from bps.mov(dcm.p,dcm_p_target)
                logger.info(f"Moving {dcm_r_target = :4f}")
                yield from bps.mov(dcm.r,dcm_r_target)
                logger.info(f"Moving {hfm_p_target = :4f}")
                yield from bps.mov(m2.p, hfm_p_target)

                #change merlin energy
                caput("XF:03IDC-ES{Merlin:1}cam1:Acquire",1)
                caput("XF:03IDC-ES{Merlin:1}cam1:OperatingEnergy", targetE)
                

            if not moveMirror == "ignore":
                yield from self.moveMirror(targetE, moveMirror)
            
            logger.info("Energy change completed")




    
    def autoUGapCalibration(self,EStart = 8.5, EEnd = 9.5, EStep = 5):

        #Make sure to go to Pt coatings

        #caput("XF:03IDA-OP{Mir:2-Ax:Y}Mtr", 4)
        #caput("XF:03IDA-OP{Mir:1-Ax:Y}Mtr", 13.5)

        #change IC1 sensivity to 5 um
        #not that we are just caput the position of the value

                            #move out,FS if in 

        if caget('XF:03IDA-OP{FS:1-Ax:Y}Mtr.VAL')<-50:

            caput('XF:03IDA-OP{FS:1-Ax:Y}Mtr.VAL', -20.)
            yield from bps.sleep(10)
            #caput('XF:03IDC-OP{Stg:CAM6-Ax:X}Mtr.VAL', -50)



        caput("XF:03IDC-CT{SR570:1}sens_num.VAL",3)
        caput("XF:03IDC-CT{SR570:1}sens_unit.VAL",2)

        #open ssa2        
        caput('XF:03IDC-OP{Slt:SSA2-Ax:XAp}Mtr.VAL', 2)
        caput('XF:03IDC-OP{Slt:SSA2-Ax:YAp}Mtr.VAL', 2)



        #setm2pf to 10
        caput("XF:03IDA-OP{HFM:1-Ax:PF}Mtr.VAL", 10)

        ic1_init = sclr1_ch2.get()
        ugap_offset = 0
        df = pd.DataFrame(columns = ["Time Stamp","energy","harmonic","ugap","dcmPitch",'dcmRoll', "hfmPitch", "IC1"], dtype = "object")

        ePoints = np.linspace(EStart, EEnd, EStep+1)
        df["energy"] = ePoints

        print(df.head())
        yield from bps.sleep(2)

        for i in tqdm.tqdm(range(len(ePoints)),desc = 'Undulator Energy Calibration'):

            yield from check_for_beam_dump(1000)

        #for i in range(len(ePoints)):

            if sclr2_ch2.get() < ic1_init*0.25:
                raise RuntimeError ("Ion chamber value dropped; aborting calibration")

            if sclr2_ch2.get() < ic1_init*0.5:
                yield from bps.mov(ssa2.vgap,1)
                yield from fluxOptimizerScan(dcm.r,-0.03, 0.03, 12, ic = sclr2_ch2, moveToMax = True)
                yield from bps.mov(ssa2.vgap,2)

            target_e = df["energy"][i]

            gap_, hrm = Energy.gap(target_e)
            gap = gap_+ np.around(ugap_offset,1)

            if abs(gap-ugap.position)>2000 and gap<5200 and gap>10000:
                raise ValueError ("Incorrect gap calculation")
            else:
                logger.info(f"Moving gap = {gap}")
                yield from bps.mov(ugap, gap)
                logger.info("Gap moved")

            yield from bps.mov(e,target_e)
            yield from self.moveMirror(target_e)
            yield from bps.sleep(2)
            yield from Energy.fluxOptimizerScan(ugap,-40, 40, 40, ic = sclr2_ch2, moveToMax = True)
            #yield from fluxOptimizerScan(ugap,-5, 5, 10, ic = xbpm, moveToMax = True)
            #if i%2 == 0

            dcm_p_target = self.calculatePitch(target_e)
            yield from bps.mov(dcm.p,dcm_p_target)
            
            logger.info("performing m2_p course centering")
            yield from bps.mov(ssa2.hgap,2, ssa2.vgap,2)
            yield from Energy.fluxOptimizerScan(m2.p,-0.005, 0.005, 10, ic = sclr2_ch2, moveToMax = True)
            yield from Energy.fluxOptimizerScan(dcm.p,-0.01, 0.01, 10, ic = sclr2_ch2, moveToMax = True)

            '''
            yield from bps.mov(ssa2.hgap,0.1,ssa2.vgap,2 )
            yield from fluxOptimizerScan(m2.pf,-0.2, 0.2, 10, ic = sclr2_ch2, moveToMax = True)
            yield from bps.mov(ssa2.hgap,2.0, ssa2.vgap,0.1)
            yield from fluxOptimizerScan(dcm.p,-0.005, 0.005, 10, ic = sclr2_ch2, moveToMax = True)
            '''
            
            
            logger.info("optimize beam at ssa2")
            yield from find_beam_at_ssa2(500,2)
            m2_p = m2.p.position
            yield from bps.mov(m2.pf, 10)
            yield from bps.mov(m2.p,m2_p)
            yield from bps.sleep(5)


            df["Time Stamp"].at[i] = pd.Timestamp.now()
            df['harmonic'].at[i] = hrm
            df['ugap'].at[i] = ugap.position
            df['dcmPitch'].at[i] = dcm.p.position
            df['dcmRoll'].at[i] = dcm.r.position
            df['hfmPitch'].at[i] = m2.p.position
            df['IC1'].at[i] = sclr2_ch2.get()
            df.to_csv(wd+"ugap_calib_01302023.csv",float_format= '%.5f')
            plt.close('all')

            ugap_offset = ugap.position - gap_
            logger.info(f"Gap offset: {ugap_offset :.1f}")

        adj_E = df["energy"].to_numpy()/df["harmonic"].to_numpy()
        E_Ugap_fit = np.polyfit(adj_E, df["ugap"].to_numpy(),3)
        print(E_Ugap_fit)

        plt.figure()
        plt.scatter(adj_E,df["ugap"],label='data')
        plt.plot(adj_E, np.polyval(E_Ugap_fit, adj_E),'r', label='fit'+str(np.around(E_Ugap_fit,1)))
        plt.xlabel("First Order Energy (keV)")
        plt.ylabel("undulator Gap (um)")
        plt.legend()
        plt.title(f"Undulator Calib_{pd.Timestamp.now().month}_{pd.Timestamp.now().year}")
        plt.show()

    @staticmethod
    def fluxOptimizerScan(motor,rel_start, rel_end, steps, ic = sclr2_ch2, moveToMax = True):

  
        MtrPos = motor.position


        x = np.linspace(MtrPos+rel_start, MtrPos+rel_end, steps+1)
        y = np.arange(steps+1)

        
        for i in y:

            yield from bps.mov(motor, x[i])
            
            if motor == m2.p:
                yield from bps.sleep(4)

            else:

                yield from bps.sleep(2)
            
            if ic==xbpm:
                y[i] = caget("XF:03ID-BI{EM:BPM1}SumAll:MeanValue_RBV")

            else:
            
                y[i] = ic.get()


        peakPos = x[y == np.max(y)][-1]

        plt.figure()
        plt.title(motor.name)
        plt.plot(x,y)


        if moveToMax:

            yield from bps.mov(motor, peakPos)
        
        else:
            yield from bps.mov(motor, MtrPos)
            #print(peakPos)
            return peakPos

def foil_calib_scan(startE, endE,saveLogFolder):
    
    energies = np.arange(startE,endE,0.0005)
    
    print(len(energies))
    
    e_list = pd.DataFrame()
    e_list['TimeStamp'] = pd.Timestamp.now()
    e_list['energy'] = energies
    e_list['E Readback'] = energies
    e_list['IC3'] = sclr2_ch4.get()
    e_list['IC0'] = sclr2_ch2.get()

    print(e_list.head())
    
    time_ = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


    for i,en in tqdm.tqdm(enumerate(energies)):
        print (i/len(energies))

        yield from Energy.move(en, moveMonoPitch=False, moveMirror = "ignore")
        yield from bps.sleep(2)
        e_list['TimeStamp'].at[i] = pd.Timestamp.now()
        e_list['IC3'].at[i] = sclr2_ch4.get() 
        e_list['IC0'].at[i] = sclr2_ch2.get()
        e_list['E Readback'].at[i] = e.position #add real energy to the dataframe

        filename = f'HXN_nanoXANES_calib_{time_}.csv'
        #filename = f'HXN_nanoXANES_calib.csv'
        e_list.to_csv(os.path.join(saveLogFolder, filename), float_format= '%.5f')

        

    plt.figure()
    spec = -1*np.log(e_list['IC3'].to_numpy()/e_list['IC0'].to_numpy())
    plt.plot(e_list['E Readback'], spec)
    plt.plot(e_list['E Readback'], np.gradient(spec))
    plt.savefig(os.path.join(saveLogFolder, filename))
    plt.show()

def peak_hfm_pitch(fine = False, tweak_range = 0.005):

    if fine:
        yield from Energy.fluxOptimizerScan(m2.pf,-1*tweak_range ,tweak_range,10)
    else:
        yield from Energy.fluxOptimizerScan(m2.p,-1*tweak_range,tweak_range,10)


def peak_dcm_roll(tweak_range = 0.005):

    yield from Energy.fluxOptimizerScan(dcm.r,-1*tweak_range,tweak_range,10)

def center_ssa2(ic = sclr2_ch4):

    yield from Energy.fluxOptimizerScan(ssa2.vcen,-0.05,0.05,10, ic = ic)
    yield from Energy.fluxOptimizerScan(ssa2.hcen,-0.05,0.05,10, ic = ic)
    yield from Energy.fluxOptimizerScan(ssa2.vcen,-0.02,0.02,10, ic = ic)


def find_beam_at_ssa2(ic1_target_k = 500, max_iter = 3):
    
    #move out FS
    caput('XF:03IDA-OP{FS:1-Ax:Y}Mtr.VAL', -20.)
    yield from bps.sleep(10)
    caput("XF:03IDA-BI{FS:1-CAM:1}cam1:Acquire",0)

    #move in CAM06
    caput('XF:03IDC-OP{Stg:CAM6-Ax:X}Mtr.VAL', 0)
    caput("XF:03IDC-ES{CAM:06}cam1:Acquire",1)

    #b shutter open 
    caput("XF:03IDB-PPS{PSh}Cmd:Opn-Cmd", 1)

    #get ic1 sensitivity and unit
    ic_sens = caget("XF:03IDC-CT{SR570:1}sens_num.VAL")
    ic_unit = caget("XF:03IDC-CT{SR570:1}sens_unit.VAL")
    
    #change IC1 sensivity to 5 um
    #caput the position of the value
    caput("XF:03IDC-CT{SR570:1}sens_num.VAL",2)
    caput("XF:03IDC-CT{SR570:1}sens_unit.VAL",2)

    #close b shutter, so that first iter works
    caput("XF:03IDB-PPS{PSh}Cmd:Cls-Cmd", 1)
    yield from bps.sleep(2)

    iter=0

    while caget("XF:03IDC-ES{Sclr:2}_cts1.B")<ic1_target_k*1000 and iter<max_iter:
        
        #b shutter open 
        caput("XF:03IDB-PPS{PSh}Cmd:Opn-Cmd", 1)
        yield from bps.sleep(2)
        
        #fully open ssa2
        yield from bps.mov(ssa2.hgap, 2, ssa2.vgap, 2) 

        #hfm
        yield from peak_hfm_pitch(fine = False, tweak_range = 0.05)
        yield from peak_dcm_roll(0.05)

        if not caget("XF:03IDC-ES{Sclr:2}_cts1.B")<1000:

            yield from peak_hfm_pitch(fine = False, tweak_range = 0.02)
            yield from bps.mov(ssa2.hgap, 0.1)
            yield from peak_hfm_pitch(fine = False, tweak_range = 0.005)
            yield from peak_hfm_pitch(fine = True, tweak_range = 0.2)
        
        #dcm_roll
        yield from bps.mov(ssa2.hgap, 2 ,ssa2.vgap, 0.1)
        yield from peak_dcm_roll(0.05)
        
        #fully open ssa2
        yield from bps.mov(ssa2.hgap, 2, ssa2.vgap, 2)
        iter += 1

        plt.close("all")

    #change back to initial sensitivity 
    caput("XF:03IDC-CT{SR570:1}sens_num.VAL",ic_sens)


def find_beam_at_cam11():

    #close c shutter
    caput("XF:03IDC-ES{Zeb:2}:SOFT_IN:B0", 0)
    
    #cam06 out
    caput('XF:03IDC-OP{Stg:CAM6-Ax:X}Mtr.VAL', -50)
    caput("XF:03IDC-ES{CAM:06}cam1:Acquire",0)

    yield from bps.mov(ssa2.hgap, 2, ssa2.vgap, 2)
    yield from bps.mov(s5.hgap, 4, s5.vgap, 4)

    yield from go_det("cam11")

    zp_osa_pos = caget("XF:03IDC-ES{ANC350:5-Ax:1}Mtr.VAL")

    if zp_osa_pos<100:

        caput("XF:03IDC-ES{ANC350:5-Ax:1}Mtr.VAL", zp_osa_pos+2700)

    #open c shutter
    caput("XF:03IDC-ES{Zeb:2}:SOFT_IN:B0", 1)

    #move beam stop
    zp_bsx_pos = caget("XF:03IDC-ES{ANC350:8-Ax:1}Mtr.VAL")

    if zp_bsx_pos<20:

        caput("XF:03IDC-ES{ANC350:8-Ax:1}Mtr.VAL", zp_bsx_pos+100)
    






Energy = HXNEnergy(ugap,e,dcm.p, "ic3", wd+"ugap_calib.csv")







