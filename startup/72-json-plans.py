print(f"Loading {__file__!r} ...")

import json
import numpy as np
import tqdm

det_dict = {"dets1":dets1,
            "dets2":dets2,
            "dets_fs":dets_fs}


#tomo_use_panda = False

def make_mll_tomo_plan(save_as = "/nsls2/data/hxn/legacy/user_macros/HXN_GUI/Scan/temp_files/mll_tomo_params.json" ):

    mll_tomo_scan = {

                    "flying_panda":True,
                    "angle_info":{"start":-90,
                                "end":90,
                                "angle_step":2},

                    "fly2d_scan":{'det':'dets1',
                                "x_start":-1,
                                "x_end":1,
                                "x_num":100,
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
                            "move_coarse":False,
                            "offset":0,
                            "zero_before_scan":false},

                    "yalign":{"do_align":True,
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5,
                            "move_coarse":False,
                            "offset":0,
                            "zero_before_scan":false},

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
                            "move_y":True,
                            "x_offest":0,
                            "y_offset":0,
                            "zero_before_scan":false
                            },

                    "stop_iter":False,
                    "add_angles":[ ],
                    "remove_angles":[-90,-91],
                    "stop_pdf":False,
                    "pdf_elems":["Au_L"],
                    "pause_scan":False,
                    "test":False,
                    "do_y_offset": True,
                    "ic_threshold":0.9,
                    "th_init": 58.99987099999999,
                    "y_init": -0.18154,
                    "scan_label":"HXN_Tomo_Scan",
                    "Notes": 'My sample # 10'

                }


    with open(save_as,"w") as fp:
            json.dump(mll_tomo_scan,fp, indent=6)

    fp.close()

def align_scan(mtr,start,end,num,exp,elem_, align_with="line_center",
               threshold = 0.5,move_coarse = False, neg_flag = False, offset = 0,
               tomo_use_panda= False, reset_piezos_to_zero = False):

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
    zp_flag = False
    
    if mtr.name.startswith('zp'):
        zp_flag = True
    
    uni_conv = 1
    if zp_flag:
        uni_conv = 1000


    if reset_piezos_to_zero:
            yield from piezos_to_zero(zp_flag = zp_flag)

    if not tomo_use_panda:
        yield from fly1d(dets_fs,
                        mtr,
                        start,
                        end,
                        num,
                        exp
                        )
    else:
        print("uses panda scan")
        yield from fly2dpd([fs,xspress3],
                        mtr,
                        start,
                        end,
                        num,
                        zpssx,
                        0,
                        0,
                        1,
                        exp,
                        position_supersample = 10
                        )
    if align_with == "line_center":
        xc = return_line_center(-1,elem_,threshold, neg_flag = neg_flag)

    elif align_with == "edge":
        xc,_ = erf_fit(-1,elem_,linear_flag=False)

    else:
        xc = mtr.position
    
    if tomo_use_panda:
        xc = xc + offset

    if move_coarse:
        yield from bps.movr(eval(fly_to_coarse[mtr.name]),xc/uni_conv)
        yield from piezos_to_zero(zp_flag = zp_flag)

    else:
        yield from bps.mov(mtr,xc)


def mll_tomo_2d_scan(angle,dets_,x_start,x_end,x_num,y_start,y_end,y_num,exp, 
                     tomo_use_panda = False):
    print("mll tomo 2d scan")

    if np.abs(angle) < 44.99:

        x_start_real = x_start / np.cos(angle * np.pi / 180.)
        x_end_real = x_end / np.cos(angle * np.pi / 180.)

        if not tomo_use_panda:
            yield from fly2d(dets_,
                            dssx,
                            x_start_real,
                            x_end_real,
                            x_num,
                            dssy,
                            y_start,
                            y_end,
                            y_num,
                            exp
                            )
        else:
            yield from fly2dpd([fs,eiger2,xspress3],
                            dssx,
                            x_start_real,
                            x_end_real,
                            x_num,
                            dssy,
                            y_start,
                            y_end,
                            y_num,
                            exp,
                            position_supersample = 10
                            )

    else:

        x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
        x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
        print(x_start_real,x_end_real)

        if not tomo_use_panda:
            yield from fly2d(dets_,
                            dssz,
                            x_start_real,
                            x_end_real,
                            x_num,
                            dssy,
                            y_start,
                            y_end,
                            y_num,
                            exp
                            )
        else:
            yield from fly2dpd([fs,eiger2,xspress3],
                            dssz,
                            x_start_real,
                            x_end_real,
                            x_num,
                            dssy,
                            y_start,
                            y_end,
                            y_num,
                            exp,
                            position_supersample = 10
                            )

def align_2d_com_scan(mtr1,x_s,x_e,x_n,mtr2,y_s,y_e,y_n,exp,elem_,
                      threshold,move_x=True,move_y=True, 
                      x_offset = 0, y_offset = 0,tomo_use_panda = False,):


    """
    @aaron removed 'angle' from first argument
    scan to align samples to field of view using using fly1d scan

    mtr1, mtr2--> scanning motor, dssx, dssy, dssz etc.
    xs,xe,xn,y_s,y_e,y_n,exp --> flyscan paramters
    elem_ --> element to use for alignemnt

    threshold --> threshold for center of mass
    move_x,move_y --> moves piezos to the com if true

    """
    #peform fly2d

    if not tomo_use_panda:
        yield from fly2d(dets_fs,
                        mtr1,
                        x_s,
                        x_e,
                        x_n,
                        mtr2,
                        y_s,
                        y_e,
                        y_n,
                        exp)
    else:
        yield from fly2dpd([fs,xspress3],
                        mtr1,
                        x_s,
                        x_e,
                        x_n,
                        mtr2,
                        y_s,
                        y_e,
                        y_n,
                        exp)
            

    #find com
    cx,cy = return_center_of_mass(-1,
                                elem_,
                                threshold
                                )
    
    cx = cx + x_offset
    cy = cy + y_offset

    #move if true
    if move_x:
        yield from bps.mov(mtr1,cx)
    if move_y:
        yield from bps.mov(mtr2,cy)

def mll_tomo_scan_to_loop(angle, tomo_params, ic_init, do_y_offset = True,tracking_file = None):

        #caput("XF:03IDC-ES{Merlin:2}HDF1:NDArrayPort","ROI1") #patch for merlin2 issuee

        #get parameters from json
        xalign = tomo_params["xalign"]
        yalign = tomo_params["yalign"]
        align_2d = tomo_params["align_2d_com"]
        image_scan = tomo_params["fly2d_scan"]
        dets = eval(image_scan["det"])
        elems_to_pdf = tomo_params["pdf_elems"]

        yield from bps.mov(dsth, angle)

        if do_y_offset:

            # precalculated y offset, mll only
            y_init = tomo_params["y_init"]
            th_init = tomo_params["th_init"]
            y_offset1 = sin_func(angle, 0.110, -0.586, 7.85,1.96)
            y_offset2 = sin_func(th_init, 0.110, -0.586, 7.85,1.96)
            yield from bps.mov(dssy,y_init+y_offset1-y_offset2)



        #yield from bps.mov(dssx,0,dssz,0)
        print(f'{tomo_params["flying_panda"] = }')
        #1d alignment sequence, based on angle x or z will be scanned
        if np.abs(angle) < 44.99:

            if xalign["do_align"]:
                yield from align_scan(dssx,
                                xalign["start"],
                                xalign["end"],
                                xalign["num"],
                                xalign["exposure"],
                                xalign["elem"],
                                xalign["center_with"],
                                xalign["threshold"],
                                xalign["move_coarse"],
                                xalign["neg_flag"],
                                xalign["offset"],
                                tomo_params["flying_panda"]
                                )

            #2d alignemnt using center of mass if condition is true
            elif align_2d["do_align"]:

                x_start_real = align_2d["x_start"] / np.cos(angle * np.pi / 180.)
                x_end_real = align_2d["x_end"] / np.cos(angle * np.pi / 180.)


                yield from align_2d_com_scan(dssx,
                                                x_start_real,
                                                x_end_real,
                                                align_2d["x_num"],
                                                dssy,
                                                align_2d["y_start"],
                                                align_2d["y_end"],
                                                align_2d["y_num"],
                                                align_2d["exposure"],
                                                align_2d["elem"],
                                                align_2d["threshold"],
                                                align_2d["move_x"],
                                                align_2d["move_y"],
                                                align_2d["x_offset"],
                                                align_2d["y_offset"],
                                                tomo_params["flying_panda"])

            else:
                pass

        else:

            if xalign["do_align"]:
                yield from align_scan(dssz,
                                xalign["start"],
                                xalign["end"],
                                xalign["num"],
                                xalign["exposure"],
                                xalign["elem"],
                                xalign["center_with"],
                                xalign["threshold"],
                                xalign["move_coarse"],
                                xalign["neg_flag"],
                                xalign["offset"],
                                tomo_params["flying_panda"]
                                )

            #2d alignemnt using center of mass if condition is true
            elif align_2d["do_align"]:

                x_start_real = align_2d["x_start"] / np.abs(np.sin(angle * np.pi / 180.))
                x_end_real = align_2d["x_end"] / np.abs(np.sin(angle * np.pi / 180.))

                yield from align_2d_com_scan(dssz,
                                                x_start_real,
                                                x_end_real,
                                                align_2d["x_num"],
                                                dssy,
                                                align_2d["y_start"],
                                                align_2d["y_end"],
                                                align_2d["y_num"],
                                                align_2d["exposure"],
                                                align_2d["elem"],
                                                align_2d["threshold"],
                                                align_2d["move_x"],
                                                align_2d["move_y"],
                                                align_2d["x_offset"],
                                                align_2d["y_offset"],
                                                tomo_params["flying_panda"]

                                                )
            else:
                pass

        #1d y alignemnt scan
        if yalign["do_align"]:
            yield from align_scan(dssy,
                                yalign["start"],
                                yalign["end"],
                                yalign["num"],
                                yalign["exposure"],
                                yalign["elem"],
                                yalign["center_with"],
                                yalign["threshold"],
                                yalign["move_coarse"],
                                xalign["neg_flag"],
                                yalign["offset"],
                                tomo_params["flying_panda"]
                )

        else:
            pass

        #2d scan sequence, based on angle x or z are scanned
        yield from mll_tomo_2d_scan(angle,
                                    dets,
                                    image_scan["x_start"],
                                    image_scan["x_end"],
                                    image_scan["x_num"],
                                    image_scan["y_start"],
                                    image_scan["y_end"],
                                    image_scan["y_num"],
                                    image_scan["exposure"],
                                    tomo_params["flying_panda"]
                                    )

        xspress3.unstage()

        #save images to pdf if
        if not tomo_params["stop_pdf"]:

            try:
                insert_xrf_map_to_pdf(-1,elems_to_pdf, "dsth", note = tomo_params["scan_label"])
                plt.close()

            except:
                pass
        if tracking_file is not None:
            flog = open(tracking_file,'a')
            flog.write("%d %.3f\n"%(db[-1].start['scan_id'],angle))
            flog.close()

def run_mll_tomo_json(path_to_json,tracking_file = None):


    """mll_tomo_scan by taking parameters from a json file,
    TODO add angles smartly

    """
    beamDumpOccured = False

    #open json file for angle info first
    with open(path_to_json,"r") as fp:
        tomo_params = json.load(fp)
    fp.close()
    print("json file loaded")

    #create angle list for iteration
    angle_info = tomo_params["angle_info"]
    print(angle_info)
    # angles = np.arange(angle_info["start"],
    #                     angle_info["end"]+angle_info["angle_step"],
    #                     angle_info["angle_step"]
    #                     )

    angles = np.linspace(angle_info["start"],
                        angle_info["end"],
                        int(1+abs(angle_info["end"] - angle_info["start"])/angle_info["angle_step"])
                        )
    if "split_tomo" in angle_info:
        split = int(angle_info["split_tomo"])
        angles0 = angles.copy()
        angles = angles0[0::split]
        for i in range(1,split):
            angles = np.concatenate((angles,angles0[i::split]))
    if "start_angle" in angle_info:
        angst = angle_info["start_angle"]
        for i,ang in enumerate(angles):
            if np.abs(ang-angst)<1e-3:
                angles = angles[i:]
                break

    angle_list = pd.DataFrame(columns= ['scan_id', 'angle'])

    angle_list["angle"] = angles

    print(angle_list)


    yield from bps.sleep(2)

    #get some initial parameters
    ic_0 = sclr2_ch2.get()
    th_init = dsth.position
    y_init = dssy.position
    do_y_offset = tomo_params["do_y_offset"]
    tomo_params["th_init"] = th_init
    tomo_params["y_init"] = y_init

    #opening fast shutter for initial ic3 reading
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
    yield from bps.sleep(2)

    logger.info("Reading IC3 value")

    #get the initial ic3 reading for peaking the beam
    ic_3_init =  sclr2_ch4.get()
    #open the json file to catch any updates

     #Ic values are useful for calibration

         #add real energy to the dataframe
    angle_list['energy_rbv'] = np.nan
    #add scan_id to the dataframe
    angle_list['scan_id'] = np.nan
    #record time
    angle_list['TimeStamp'] = pd.Timestamp.now()
    #record if peak beam happed before the scan
    angle_list['Peak Flux'] = False
    angle_list['IC3'] = ic_3_init
    angle_list['IC0'] = ic_0
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)


    #set the pause and stop inter keys to False before the loop
    #to reverse the abort scan and pause when using the gui
    tomo_params["stop_iter"] = False
    tomo_params["pause_scan"] = False

    with open(path_to_json,"w") as fp:
            json.dump(tomo_params,fp, indent=6)

    fp.close()
    
    fluxPeaked = False

    #loop with list of angles
    for n,angle in enumerate(tqdm.tqdm(angles,desc = 'MLL Tomo Scan')):
        yield from bps.sleep(1)

        #open the json file to catch any updates
        with open(path_to_json,"r") as fp:
            tomo_params = json.load(fp)
            fp.close()

        #stop data collection if necessary.user input taken
        if tomo_params["stop_iter"]:
            save_page()
            break

        while tomo_params["pause_scan"]:
            yield from bps.sleep(10) #check if this freezes the gui or not
            with open(path_to_json,"r") as fp:
                tomo_params = json.load(fp)
                fp.close()

            if not tomo_params["pause_scan"]:
                break

        if tomo_params["remove_angles"]==None:
            tomo_params["remove_angles"] = []

        if not tomo_params["test"]:
            if sclr2_ch2.get()<1000:
                beamDumpOccured = True
                yield from check_for_beam_dump()


            if beamDumpOccured:
                angle = angles[n-1]
                yield from bps.sleep(60)
                yield from recover_from_beamdump()
                beamDumpOccured = False

        #look for beam dump and ic3 threshold, ignores for code tests using json
        if not tomo_params["test"]:
            #opening fast shutter for initial ic3 reading
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
            yield from bps.sleep(2)
            if (sclr2_ch4.get() < (tomo_params["ic_threshold"]*ic_3_init)):
                 yield from peak_the_flux()
                 ic_0 = sclr2_ch2.get()
                 fluxPeaked = True

        if not angle in np.array(tomo_params["remove_angles"]):
            #tomo scan at a single angle
            yield from mll_tomo_scan_to_loop(angle, tomo_params,ic_0,do_y_offset = do_y_offset,tracking_file = tracking_file)

            last_sid = int(caget('XF:03IDC-ES{Status}ScanID-I'))

            #Add more info to the dataframe
            angle_list['energy_rbv'].at[n] = e.position #add real energy to the dataframe
            angle_list["angle"].at[n] = dsth.position
            angle_list['scan_id'].at[n] = int(last_sid) #add scan_id to the dataframe
            angle_list['TimeStamp'].at[n] = pd.Timestamp.now()
            ic_3 = sclr2_ch4.get()
            ic_0 = sclr2_ch2.get()
            angle_list['IC3'].at[n] = ic_3 #Ic values are useful for calibration
            angle_list['IC0'].at[n] = ic_0 #Ic values are useful for calibration
            angle_list['Peak Flux'].at[n] = fluxPeaked # recoed if peakflux was excecuted
            #angle_list['IC3_before_peak'].at[n] = ic3 #ic3 right after e change, no peaking
            fluxPeaked = False #reset
            
            #close c shutter
            caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)

            # save the DF in the loop so quitting a scan won't affect
            filename = f"{tomo_params.get('scan_label','')}_startID{int(angle_list['scan_id'][0])}.csv"
            angle_list.to_csv(os.path.join(tomo_params.get("save_log_to", "/data/users/current_user"), filename),
                              float_format= '%.5f')


        else:
            print(f"{angle} skipped")
            pass


    #TODO add angles to scan; need to be better
    #sort based on what current angle is
    if not tomo_params["add_angles"]==None:

        added_angles = tomo_params["add_angles"]

        for nn,angle in enumerate(tqdm.tqdm(added_angles,desc = 'MLL Tomo Scan; Additional Angles')):

            #open the json file to catch any updates
            with open(path_to_json,"r") as fp:
                tomo_params = json.load(fp)
                fp.close()

            #stop data collection if necessary.user input taken
            if tomo_params["stop_iter"]:
                save_page()
                break

            while tomo_params["pause_scan"]:
                yield from bps.sleep(10) #check if this freezes the gui or not
                with open(path_to_json,"r") as fp:
                    tomo_params = json.load(fp)
                    fp.close()

                if not tomo_params["pause_scan"]:
                    break

            if not angle in np.array(tomo_params["remove_angles"]):
                yield from mll_tomo_scan_to_loop(angle, tomo_params,ic_0,do_y_offset = do_y_offset)

                last_sid = int(caget('XF:03IDC-ES{Status}ScanID-I'))

                #Add more info to the dataframe
                angle_list['energy_rbv'].at[n+nn] = e.position #add real energy to the dataframe
                angle_list["angle"].at[n+nn] = dsth.position
                angle_list['scan_id'].at[n+nn] = int(last_sid) #add scan_id to the dataframe
                angle_list['TimeStamp'].at[n+nn] = pd.Timestamp.now()
                angle_list['IC3'].at[n+nn] = ic_3 #Ic values are useful for calibration
                angle_list['IC0'].at[n+nn] = ic_0 #Ic values are useful for calibration
                #angle_list['Peak Flux'].at[n] = fluxPeaked # recoed if peakflux was excecuted
                #angle_list['IC3_before_peak'].at[n] = ic3 #ic3 right after e change, no peaking
                #fluxPeaked = False #reset

                # save the DF in the loop so quitting a scan won't affect
                filename = f"hxn_zp_diff_{tomo_params.get('scan_label','')}_startID{int(angle_list['scan_id'][0])}.csv"
                angle_list.to_csv(os.path.join(tomo_params.get("save_log_to", "/data/users/current_user"), filename),
                                float_format= '%.5f')

            else:
                print(f"{angle} skipped")
                pass

    else:
        pass

    #save pdf
    save_page()

###################ZP Diffraction######################

def make_diff_plan(save_as = "/nsls2/data/hxn/legacy/user_macros/HXN_GUI/Scan/temp_files/diff_params.json" ):

    zp_diff_scan = {
                    "angle_info":{"th_motor":"zpsth",
                                  "start":70,
                                  "end":71,
                                  "num":20},

                    "fly2d_scan":{'det':'dets1',
                                "x_motor":"zpssx",
                                "x_start":-1,
                                "x_end":1,
                                "x_num":100,
                                "y_motor":"zpssy",
                                "y_start":-1,
                                "y_end":1,
                                "y_num":100,
                                "exposure":0.03},

                    "xalign":{"do_align":True,
                            "x_motor":"zpssx",
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5},

                    "yalign":{"do_align":True,
                            "y_motor":"zpssy",
                            "start":-2,
                            "end": 2,
                            "num": 100,
                            "exposure": 0.03,
                            "elem": "Fe",
                            "center_with":"line_center",
                            "threshold": 0.5},

                    "align_2d_com":{"do_align":False,
                            "x_motor":"zpssx",
                            "x_start":-2,
                            "x_end": 2,
                            "x_num": 100,
                            "y_motor":"zpssy",
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
                    "remove_angles":[],
                    "stop_pdf":False,
                    "pause_scan":False,
                    "test":False
                }


    with open(save_as,"w") as fp:
            json.dump(zp_diff_scan,fp, indent=6)

    fp.close()


def diff_scan_to_loop(angle, tomo_params, ic_init):

        #caput("XF:03IDC-ES{Merlin:2}HDF1:NDArrayPort","ROI1") #patch for merlin2 issuee

        #get parameters from json
        xalign = tomo_params["xalign"]
        yalign = tomo_params["yalign"]
        align_2d = tomo_params["align_2d_com"]
        image_scan = tomo_params["fly2d_scan"]
        dets = eval(image_scan["det"])


        yield from bps.mov(dsth, angle)


        #look for beam dump and ic3 threshold, ignores for code tests using json
        if not tomo_params["test"]:
            yield from check_for_beam_dump()

            while (sclr2_ch2.get() < (0.9*ic_init)):
                 yield from peak_bpm_y(-5,5,10)
                 yield from peak_bpm_x(-15,15,10)
                 ic_0 = sclr2_ch2.get()


        if xalign["do_align"]:
            yield from align_scan(eval(xalign["x_motor"]),
                            xalign["start"],
                            xalign["end"],
                            xalign["num"],
                            xalign["exposure"],
                            xalign["elem"],
                            xalign["center_with"],
                            xalign["threshold"],
                            )

        #2d alignemnt using center of mass if condition is true
        elif align_2d["do_align"]:


            yield from align_2d_com_scan(
                                        eval(align_2d["x_motor"]),
                                        align_2d["x_start"],
                                        align_2d["x_end"],
                                        align_2d["x_num"],
                                        eval(align_2d["y_motor"]),
                                        align_2d["y_start"],
                                        align_2d["y_end"],
                                        align_2d["y_num"],
                                        align_2d["exposure"],
                                        align_2d["elem"],
                                        align_2d["threshold"],
                                        align_2d["move_x"],
                                        align_2d["move_y"],
                                        align_2d["x_offset"],
                                        align_2d["y_offset"],
                                        tomo_params["flying_panda"])

        else:
            pass


        #1d y alignemnt scan
        if yalign["do_align"]:
            yield from align_scan(eval(yalign["y_motor"]),
                                yalign["start"],
                                yalign["end"],
                                yalign["num"],
                                yalign["exposure"],
                                yalign["elem"],
                                yalign["center_with"],
                                yalign["threshold"],
                )

        else:
            pass


        yield from fly2d(
                        dets,
                        eval(image_scan["x_motor"]),
                        image_scan["x_start"],
                        image_scan["x_end"],
                        image_scan["x_num"],
                        eval(image_scan["y_motor"]),
                        image_scan["y_start"],
                        image_scan["y_end"],
                        image_scan["y_num"],
                        image_scan["exposure"]
                        )
        #save images to pdf if
        if not tomo_params["stop_pdf"]:

            try:
                insert_xrf_map_to_pdf(-1,["Cu"], "dsth")
                plt.close()

            except:
                pass


def th_fly2d_json(path_to_json):

    #TODO make in generic for MLL and ZP

    beamDumpOccured = False

    #open json file for angle info first
    with open(path_to_json,"r") as fp:
        diff_params = json.load(fp)
    fp.close()

    print("json file loaded")
    th_motor = eval(diff_params["angle_info"]["th_motor"])
    th_start = diff_params["angle_info"]["start"]
    th_end = diff_params["angle_info"]["end"]
    num = diff_params["angle_info"]["num"]

    """relative theta"""

    #yield from shutter('open')
    init_th = th_motor.position
    th_step = (th_end - th_start) / num
    th_pos = np.linspace(init_th + th_start,init_th + th_end, num+1)
    ic_0 = sclr2_ch4.get()


    for i in tqdm.tqdm(range(num + 1), desc = 'Theta Scan'):


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
            yield from bps.sleep(60)
            yield from recover_from_beamdump()
            beamDumpOccured = False

        yield from bps.mov(th_motor,th_pos[i]) #TODO make in generic for MLL and ZP

        xalign = diff_params["xalign"]
        yalign = diff_params["yalign"]
        align_2d = diff_params["align_2d_com"]
        image_scan = tomo_params["fly2d_scan"]
        dets = eval(image_scan["det"])

        if xalign["do_align"]:
            yield from align_scan(eval(xalign["x_motor"]),
                            xalign["start"],
                            xalign["end"],
                            xalign["num"],
                            xalign["exposure"],
                            xalign["elem"],
                            xalign["center_with"],
                            xalign["threshold"],
                            )

        #2d alignemnt using center of mass if condition is true
        elif align_2d["do_align"]:


            yield from align_2d_com_scan(
                                        eval(align_2d["x_motor"]),
                                        align_2d["x_start"],
                                        align_2d["x_end"],
                                        align_2d["x_num"],
                                        eval(align_2d["y_motor"]),
                                        align_2d["y_start"],
                                        align_2d["y_end"],
                                        align_2d["y_num"],
                                        align_2d["exposure"],
                                        align_2d["elem"],
                                        align_2d["threshold"],
                                        align_2d["move_x"],
                                        align_2d["move_y"],
                                        align_2d["x_offset"],
                                        align_2d["y_offset"],
                                        tomo_params["flying_panda"])

        else:
            pass


        #1d y alignemnt scan
        if yalign["do_align"]:
            yield from align_scan(eval(yalign["y_motor"]),
                                yalign["start"],
                                yalign["end"],
                                yalign["num"],
                                yalign["exposure"],
                                yalign["elem"],
                                yalign["center_with"],
                                yalign["threshold"],)

        else:
            pass


        yield from fly2d(
                    dets,
                    eval(image_scan["x_motor"]),
                    image_scan["x_start"],
                    image_scan["x_end"],
                    image_scan["x_num"],
                    eval(image_scan["y_motor"]),
                    image_scan["y_start"],
                    image_scan["y_end"],
                    image_scan["y_num"],
                    image_scan["exposure"]
                    )
        #save images to pdf if
        if not tomo_params["stop_pdf"]:

            try:
                insert_xrf_map_to_pdf(-1,["Cu"], th_motor)
                plt.close()

            except:
                pass





