
import json


mll_tomo_scan = {   
                "angle_info":{"start":-90, 
                               "end":90, 
                               "num":90}, 

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
                           "threshold": 0.5},
                
                "yalign":{"do_align":True,
                           "start":-2,
                           "end": 2,
                           "num": 100,
                           "exposure": 0.03,
                           "elem": "Fe",
                           "center_with":"line_center",
                           "threshold": 0.5},

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
                           "move_y":True},
                
                "stop_iter":False,
                "add_angles":[91,92],
                "remove_angles":[-90,-91],
                "stop_pdf":False,
                "pause_scan":False,
                "test":False 
            }


with open("/data/user_macros/HXN_GUI/Scan/json_plans/mll_tomo/mll_tomo_params.json","w") as fp:
        json.dump(mll_tomo_scan,fp, indent=6)

fp.close()