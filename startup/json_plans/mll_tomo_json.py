import json
import numpy as np
import tqdm

det_dict = {"dets1":dets1,
            "dets2":dets2,
            "dets_fs":dets_fs}


def align_scan(mtr, start,end,num,exp,elem_, align_with="line_center", threshold = 0.5):

    """
    scan to align samples to field of view using using fly1d scan 

    mtr--> scanning motor, dssx, dssy, dssz etc.
    start,end,num,exp --> flyscan paramters
    elem_ --> element to use for alignemnt
    align_with --> choose bettween "edge" or "line_center"
    threshold --> threshold for line centering
    
    """

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
        xc = mtr.position

    yield from bps.mov(mtr,xc)

def mll_tomo_2d_scan(angle,dets_,x_start,x_end,x_num,y_start,y_end,y_num,exp):
    print("mll tomo 2d scan")

    if np.abs(angle) < 44.99:
                
        x_start_real = x_start / np.cos(angle * np.pi / 180.)
        x_end_real = x_end / np.cos(angle * np.pi / 180.)

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

        x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
        x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
        print(x_start_real,x_end_real)

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

def align_2d_com_scan(angle,mtr1,x_s,x_e,x_n,mtr2,y_s,y_e,y_n,exp,elem_,threshold,move_x=True,move_y=True):


    """
    scan to align samples to field of view using using fly1d scan 

    mtr1, mtr2--> scanning motor, dssx, dssy, dssz etc.
    xs,xe,xn,y_s,y_e,y_n,exp --> flyscan paramters
    elem_ --> element to use for alignemnt

    threshold --> threshold for center of mass
    move_x,move_y --> moves piezos to the com if true

    """
    #peform fly2d
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

    #find com
    cx,cy = return_center_of_mass(-1,
                                elem_,
                                threshold
                                )
    #move if true
    if move_x:
        yield from bps.mov(mtr1,cx)
    if move_y:
        yield from bps.mov(mtr2,cy)



def mll_tomo_json(path_to_json):


    """mll_tomo_scan by taking parameters from a json file,
    TODO add angles smartly
    
    """

    #open json file for angle info first
    with open(path_to_json,"r") as fp:
        tomo_params = json.load(fp)
    fp.close()
    print("json file loaded")

    #create angle list for iteration
    angle_info = tomo_params["angle_info"]
    print(angle_info)
    angles = np.linspace(angle_info["start"], 
                         angle_info["end"],
                         angle_info["num"]+1
                         )
    print(angles)
    
    #get some initial parameters 
    ic_0 = sclr2_ch2.get()
    th_init = dsth.position
    y_init = dssy.position


    #loop with list of angles
    for angle in tqdm.tqdm(angles,desc = 'MLL Tomo Scan'):
        
        #open the json file to catch any updates 
        with open(path_to_json,"r") as fp:
            tomo_params = json.load(fp)
            fp.close()

        #get parameters from json
        xalign = tomo_params["xalign"]
        yalign = tomo_params["yalign"]
        align_2d = tomo_params["align_2d_com"]
        image_scan = tomo_params["fly2d_scan"]
        dets = det_dict[image_scan["det"]]

        #stop data collection if necessary.user input taken 
        if tomo_params["stop_iter"]:
            save_page()
            break
            
        #ignore an angle using json 
        if not angle in np.array(tomo_params["remove_angles"]):
            yield from bps.mov(dsth, angle)

            # precalculated y offset, mll only 
            #y_offset1 = sin_func(angle, 0.110, -0.586, 7.85,1.96)
            #y_offset2 = sin_func(th_init, 0.110, -0.586, 7.85,1.96)
            #yield from bps.mov(dssy,y_init+y_offset1-y_offset2)


            #look for beam dump and ic3 threshold, ignores for code tests using json
            if not tomo_params["test"]:
                yield from check_for_beam_dump()

                while (sclr2_ch2.get() < (0.95*ic_0)):
                    yield from peak_bpm_y(-5,5,10)
                    yield from peak_bpm_x(-15,15,10)
                    ic_0 = sclr2_ch4.get()
            

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
                                                 align_2d["move_y"],)

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
                                                 align_2d["move_y"]
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
                                        image_scan["exposure"]
                                        )

            merlin1.unstage()
            xspress3.unstage()

            #save images to pdf if
            if not tomo_params["stop_pdf"]:

                try:
                    insert_xrf_map_to_pdf(-1,["Au_L"], "dsth")
                
                except:
                    pass

            while tomo_params["pause_scan"]:
                yield from bps.sleep(10)
                with open(path_to_json,"r") as fp:
                    tomo_params = json.load(fp)
                    fp.close() 

                if not tomo_params["pause_scan"]:   
                    break

        else:
            print(f"{angle} skipped")
            pass

    
    #save pdf
    save_page()
