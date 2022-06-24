import logging
import numpy as np
import pandas as pd
import pyqtgraph as pg
import matplotlib.pyplot as plt
from epics import caget, caput

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

class HXNEnergy():
    
    def __init__(self, ugap_,bragg_e, dcm_pitch, ic1, calib_file_csv):

        self.ugap = ugap_
        self.bragg_e = bragg_e
        self.dcm_pitch = dcm_pitch
        self.ic1 = ic1
        self.calib_file = calib_file_csv
        self.df = pd.read_csv(self.calib_file)
        self.calibrate_ugap()
        self.calibrate_dcm_pitch()

        self.defineFeedBackPVs()

    def defineFeedBackPVs(self):
        self.hcm_pf_sts = "XF:03IDC-CT{FbPid:01}PID:on"
        self.hfm_pf_sts = "XF:03IDC-CT{FbPid:02}PID:on"

        self.xbpm_x_rbv = "XF:03ID-BI{EM:BPM1}PosX:MeanValue_RBV"
        self.xbpm_y_rbv = "XF:03ID-BI{EM:BPM1}PosY:MeanValue_RBV"

        self.xbpm_x_val = "XF:03ID-BI{EM:BPM1}fast_pidX.VAL"
        self.xbpm_x_val = "XF:03ID-BI{EM:BPM1}fast_pidY.VAL"

    def calcGap(self,E,harmonics = 5, offset = 0):
        E1 = E/harmonics
        calc_gap =  np.polyval(self.ugap_coeffs, E1) + offset
        return (np.around(calc_gap,1))

    def gap(self,targetE, harmonicsChoice = "auto"):
        
        "Estimate a suitable gap for target energy"

        if harmonicsChoice == "auto":
        
            harmonics = np.array([3,5,7,9]) #harmonics options

            if targetE>5.8 and targetE<25:     # if dialed a wrong number
                
                opt = np.array([self.calcGap(targetE, harmonics = hm) for hm in harmonics])
                idx = np.where(np.logical_and(opt>=5600, opt<=10000)) #5400 beacuse of the ugap scan limit
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
    
        EnergyToMirror = {(5,9):'Si',(9,15):'Cr',(15,20):'Rh',(20,25):'Pt'}
    
        if not 5.5 <= e <= 25:
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

    def calculatePitch(self,targetE, offset = 0):
        
        calc_pitch =  np.polyval(self.dcm_p_coeffs, targetE) + offset
        return calc_pitch 

    def move(self,targetE, harmChoice = "auto", moveMonoPitch = True, moveMirror = "auto") :
        
        bbpm_auto = "XF:03ID{XBPM:17}AutoFbEn-Cmd"
        bbpm_x = "XF:03ID-BI{EM:BPM1}fast_pidX.FBON"
        bbpm_y = "XF:03ID-BI{EM:BPM1}fast_pidY.FBON"

        gap, hrm = Energy.gap(targetE,harmChoice)
        dcm_p_target = self.calculatePitch(targetE)
        logger.info(f"target pitch = {dcm_p_target}")

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

                if caget(bbpm_auto) or caget(bbpm_x) or caget(bbpm_y):
                    print(True)
                    caput(bbpm_auto,0)
                    yield from bps.sleep(2)
                    caput(bbpm_x, 0)
                    caput(bbpm_y,0)
                    yield from bps.sleep(2)
                    
                    caput(self.hcm_pf_sts,0)
                    yield from bps.sleep(2)
                    caput(self.hfm_pf_sts,0)
                    yield from bps.sleep(2)
                    logger.info(f"Moving Mono Pitch to {dcm_p_target}")
                    yield from bps.mov(dcm.p,dcm_p_target)
                    yield from bps.sleep(2)
                    caput(self.hcm_pf_sts,1)
                    yield from bps.sleep(2)
                    caput(self.hfm_pf_sts,1)
                        
                    yield from bps.sleep(2)
                    caput(bbpm_auto,1)
                    yield from bps.sleep(2)
                    caput(bbpm_x, 1)
                    yield from bps.sleep(2)
                    caput(bbpm_y,1)
                    yield from bps.sleep(5)

                else:

                    caput(self.hcm_pf_sts,0)
                    yield from bps.sleep(2)
                    caput(self.hfm_pf_sts,0)
                    yield from bps.sleep(2)
                    logger.info(f"Moving Mono Pitch to {dcm_p_target}")
                    yield from bps.mov(dcm.p,dcm_p_target)

            if not moveMirror == "ignore":
                yield from self.moveMirror(targetE, moveMirror)


            '''
            caput(self.hcm_pf_sts,1)
            yield from bps.sleep(2)
            caput(self.hfm_pf_sts,1)
            yield from bps.sleep(2)
            '''
            
            logger.info("Energy change completed")


    def calibrate_ugap(self):

        adj_E = self.df["energy"].to_numpy()/self.df["harmonic"].to_numpy()
        self.ugap_coeffs = np.polyfit(adj_E, self.df["ugap"].to_numpy(), 3)

    def calibrate_dcm_pitch(self):
        
        en, dcm_p = self.df["energy"].to_numpy(),self.df["dcmPitch"].to_numpy()
        self.dcm_p_coeffs = np.polyfit(en, dcm_p, 3)
        #logger.info(f" DCM_P Coefficients: {self.dcm_p_coeffs}")

    def calibrate_hfm_pitch(self):
        
        self.dcm_p_coeffs = np.polyfit(self.df["energy"].to_numpy(), self.df["hfmPitch"].to_numpy(), 3)

        pass

    def calibrate_dcm_roll(self):
        pass

    @staticmethod
    def autoUGapCalibration(EStart = 12, EEnd = 12.5, EStep = 10):

        #Make sure to go to Pt coatings

        caput("XF:03IDA-OP{Mir:2-Ax:Y}Mtr", 4)
        caput("XF:03IDA-OP{Mir:1-Ax:Y}Mtr", 13.5)

        #open ssa2        
        caput('XF:03IDC-OP{Slt:SSA2-Ax:XAp}Mtr.VAL', 2.5)
        caput('XF:03IDC-OP{Slt:SSA2-Ax:YAp}Mtr.VAL', 2.5)

        #move out,FS if in 

        caput('XF:03IDA-OP{FS:1-Ax:Y}Mtr.VAL', -20.)
        #caput('XF:03IDC-OP{Stg:CAM6-Ax:X}Mtr.VAL', -50)

        ic1_init = sclr1_ch2.get()
        ugap_offset = 0
        df = pd.DataFrame(columns = ["Time Stamp","energy","harmonic","ugap","dcmPitch",'dcmRoll', "hfmPitch", "IC1"], dtype = "object")

        ePoints = np.linspace(EStart, EEnd, EStep+1)
        df["energy"] = ePoints

        print(df.head())
        yield from bps.sleep(2)

        for i in range(len(ePoints)):

            if sclr2_ch2.get() < ic1_init*0.25:
                raise RuntimeError ("Ion chamber value dropped; aborting calibration")

            if sclr2_ch2.get() < ic1_init*0.5:
                yield from bps.mov(ssa2.vgap,1)
                yield from fluxOptimizerScan(dcm.r,-0.03, 0.03, 12, ic = sclr2_ch2, moveToMax = True)
                yield from bps.mov(ssa2.vgap,2.)

            gap_, hrm = Energy.gap(df["energy"][i])
            gap = gap_+ np.around(ugap_offset,1)

            if abs(gap-ugap.position)>2000 and gap<5200 and gap>10000:
                raise ValueError ("Incorrect gap calculation")
            else:
                logger.info(f"Moving gap = {gap}")
                yield from bps.mov(ugap, gap)
                logger.info("Gap moved")

            yield from bps.mov(e,(df["energy"][i]))
            yield from bps.sleep(2)
            yield from fluxOptimizerScan(ugap,-60, 60, 12, ic = sclr2_ch2, moveToMax = True)
            #yield from fluxOptimizerScan(ugap,-5, 5, 10, ic = xbpm, moveToMax = True)
            #if i%2 == 0
            logger.info("performing m2_p centering")
            yield from bps.mov(ssa2.hgap,1)
            yield from fluxOptimizerScan(m2.p,-0.005, 0.005, 10, ic = sclr2_ch2, moveToMax = True)
            yield from bps.mov(ssa2.hgap,2.5)
            yield from fluxOptimizerScan(dcm.p,-0.005, 0.005, 10, ic = sclr2_ch2, moveToMax = True)

            df["Time Stamp"].at[i] = pd.Timestamp.now()
            df['harmonic'].at[i] = hrm
            df['ugap'].at[i] = ugap.position
            df['dcmPitch'].at[i] = dcm.p.position
            df['dcmRoll'].at[i] = dcm.r.position
            df['hfmPitch'].at[i] = m2.p.position
            df['IC1'].at[i] = sclr2_ch2.get()
            df.to_csv("ugap_calib.csv",float_format= '%.5f')
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
        


Energy = HXNEnergy(ugap,e,dcm.p, "ic3", "/nsls2/data/hxn/shared/config/bluesky/profile_collection/startup/ugap_calib.csv")




