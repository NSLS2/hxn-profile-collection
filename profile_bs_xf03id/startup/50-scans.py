# vim: sw=4 ts=4 sts expandtab smarttab

# TODO: figure this out
olog_client = None

from bluesky.standard_config import *
import bluesky
from bluesky.run_engine import DocumentNames
from hxnfly.bs import (BSFlyScan1D, Flyer)


class HxnScanNumberPrinter:
    def __init__(self):
        self._last_start = None

    def __call__(self, name, doc):
        if name == DocumentNames.start:
            self._last_start = doc
        if self._last_start is None:
            return
        if name in (DocumentNames.start, DocumentNames.stop):
            print('Scan ID: {scan_id} [{uid}]'.format(**self._last_start))


gs.RE.subscribe('start', HxnScanNumberPrinter())


default_detectors = [det_sclr1, det_beamstatus,
                     sclr1_ch2, sclr1_ch3, sclr1_ch4,
                     det1_Pt, det2_Pt, det3_Pt, det1_Al, det2_Al, det3_Al,
                     det1_Si, det2_Si, det3_Si, det1_Gd, det2_Gd, det3_Gd,
                     det1_YSZ, det2_YSZ, det3_YSZ, det1_Ca, det2_Ca, det3_Ca,
                     det1_Ti, det2_Ti, det3_Ti, det1_V, det2_V, det3_V,
                     det1_Cr, det2_Cr, det3_Cr, det1_Mn, det2_Mn, det3_Mn,
                     det1_Fe, det2_Fe, det3_Fe, det1_Co, det2_Co, det3_Co,
                     det1_Ni, det2_Ni, det3_Ni, det1_Cu, det2_Cu, det3_Cu,
                     det1_Ce, det2_Ce, det3_Ce, det1_Au, det2_Au, det3_Au,
                     ssx_rbv, ssy_rbv, ssz_rbv, t_base, t_sample, t_vlens,
                     t_hlens]


user_detectors = [xspress3.filestore, timepix1.filestore]

beamline_config_pvs = [ssx_rbv, ssy_rbv, ssz_rbv, t_base, t_sample, t_vlens]
project_info = 'project_information'
scan_metadata = {'sample': {'owner': 'bnl',
                            'name': 'sample_name'}}
scan_owner = 'xf03id'

def synchronize(detectors, integration_time):
    # EpicsScaler calls it preset_time;
    # AreaDetector calls it exposure_time.
    integration_time = float(integration_time)

    attrs = ['preset_time', 'exposure_time']
    for det in detectors:
        for attr in attrs:
            try:
                setattr(det, attr, integration_time)
            except AttributeError:
                pass
            # TODO Raise better.

    # For a step scan using the zebra or scaler 1, the scaler is the master
    # Scaler LNE -> Zebra IN1_TTL -> Pulse1 -> fanout
    if zebra in detectors or det_sclr1 in detectors:
        sclr1._preset_time.put(integration_time)
        zebra.preset_time = integration_time

    for tpx in [timepix1, timepix2]:
        if tpx.filestore in detectors:
            tpx.acquire_time.put(integration_time)
            tpx.acquire_period.put(integration_time)

    for merlin in [merlin1, ]:
        if merlin.filestore in detectors:
            merlin.acquire_time.put(integration_time)
            merlin.acquire_period.put(integration_time + 0.001)


def movr_hth(angle):
    radian = angle*pi/180.0
    correction = -1.*tan(radian)*34376.6
    movr(hth, angle)
    movr(hx,correction)


def movr_ssx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)

def movr_sx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_sz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_ssz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)
