from ophyd.controls import PVPositioner

# Undulator

ugap = PVPositioner('SR:C3-ID:G1{IVU20:1-Mtr:2}Inp:Pos',
                    readback='SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.RBV',
                    act='SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go', act_val=1,
                    stop='SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.STOP', stop_val=1,
                    done='SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.DMOV', done_val=1,
                    put_complete=False,
                    settle_time=1.5,
                    name='ugap')

# Front End Slits (Primary Slits)

fe_tb = PVPositioner('FE:C03A-OP{Slt:1-Ax:T}Mtr.VAL',
                     readback='FE:C03A-OP{Slt:1-Ax:T}Mtr.RBV',
                     stop='FE:C03A-OP{Slt:1-Ax:T}Mtr.STOP',
                     stop_val=1, put_complete=True,
                     name='fe_tb')

fe_bb = PVPositioner('FE:C03A-OP{Slt:2-Ax:B}Mtr.VAL',
                     readback='FE:C03A-OP{Slt:2-Ax:B}Mtr.RBV',
                     stop='FE:C03A-OP{Slt:2-Ax:B}Mtr.STOP',
                     stop_val=1, put_complete=True,
                     name='fe_bb')

fe_ib = PVPositioner('FE:C03A-OP{Slt:2-Ax:I}Mtr.VAL',
                     readback='FE:C03A-OP{Slt:2-Ax:I}Mtr.RBV',
                     stop='FE:C03A-OP{Slt:2-Ax:I}Mtr.STOP',
                     stop_val=1, put_complete=True,
                     name='fe_ib')

fe_ob = PVPositioner('FE:C03A-OP{Slt:1-Ax:O}Mtr.VAL',
                     readback='FE:C03A-OP{Slt:1-Ax:O}Mtr.RBV',
                     stop='FE:C03A-OP{Slt:1-Ax:O}Mtr.STOP',
                     stop_val=1, put_complete=True,
                     name='fe_ob')


