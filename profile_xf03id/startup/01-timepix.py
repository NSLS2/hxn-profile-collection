from ophyd.controls.areadetector.detectors import (AreaDetector, ADSignal)


class TimepixDetector(AreaDetector):
    _html_docs = []

    tpx_corrections_dir = ADSignal('TPXCorrectionsDir', string=True)
    tpx_dac = ADSignal('TPXDAC_RBV', rw=False)
    tpx_dac_available = ADSignal('TPX_DACAvailable')
    tpx_dac_file = ADSignal('TPX_DACFile', string=True)
    tpx_dev_ip = ADSignal('TPX_DevIp', has_rbv=True)
    _tpx_extended_frame = ADSignal('TPX_ExtendedFrame', has_rbv=True)
    _tpx_extended_frame_no = ADSignal('TPX_ExtendedFrameNo')
    _tpx_extended_frame_yes = ADSignal('TPX_ExtendedFrameYes')

    @property
    def tpx_extended_frame(self):
        return self._tpx_extended_frame

    @tpx_extended_frame.setter
    def tpx_extended_frame(self, value):
        if value:
            self._tpx_extended_frame_yes.put(1)
        else:
            self._tpx_extended_frame_no.put(1)

    tpx_frame_buff_index = ADSignal('TPXFrameBuffIndex')
    tpx_hw_file = ADSignal('TPX_HWFile', string=True)
    tpx_initialize = ADSignal('TPX_Initialize', has_rbv=True)
    tpx_load_dac_file = ADSignal('TPXLoadDACFile')
    tpx_num_frame_buffers = ADSignal('TPXNumFrameBuffers', has_rbv=True)
    tpx_pix_config_file = ADSignal('TPX_PixConfigFile', string=True)
    tpx_reset_detector = ADSignal('TPX_resetDetector')

    tpx_raw_image_number = ADSignal('TPXImageNumber')
    tpx_raw_prefix = ADSignal('TPX_DataFilePrefix', string=True)
    tpx_raw_path = ADSignal('TPX_DataSaveDirectory', string=True)

    _tpx_save_to_file = ADSignal('TPX_SaveToFile', has_rbv=True)
    _tpx_save_to_file_no = ADSignal('TPX_SaveToFileNo')
    _tpx_save_to_file_yes = ADSignal('TPX_SaveToFileYes')

    @property
    def tpx_save_raw(self):
        return self._tpx_save_to_file

    @tpx_save_raw.setter
    def tpx_save_raw(self, value):
        if value:
            self._tpx_save_to_file_yes.put(1)
        else:
            self._tpx_save_to_file_no.put(1)

    tpx_start_sophy = ADSignal('TPX_StartSoPhy', has_rbv=True)
    tpx_status = ADSignal('TPXStatus_RBV', rw=False)
    tpx_sync_mode = ADSignal('TPXSyncMode', has_rbv=True)
    tpx_sync_time = ADSignal('TPXSyncTime', has_rbv=True)
    tpx_system_id = ADSignal('TPXSystemID')
    tpx_trigger = ADSignal('TPXTrigger')
