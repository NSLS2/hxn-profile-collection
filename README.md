# HXN Profile Collection

Ophyd/Bluesky startup configuration for NSLS-II's Hard X-ray Nanoprobe (HXN) beamline.

## Overview

This repository contains the complete beamline control system setup, including:

- **Motor & Detector Definitions**: Ophyd devices for beamline hardware
- **Plan Functions**: High-level experiment execution routines
- **Data Acquisition**: Flyscan and step-scan implementations
- **Alignment Tools**: Multi-Layer Lens (MLL) and Zone Plate (ZP) optimization
- **Data Processing**: Real-time analysis and visualization

## Quick Start

### Environment Setup

```bash
# Install dependencies using pixi
pixi install

# Activate environment
pixi shell
```

### Running Experiments

The startup sequence loads all beamline components in order:

1. **00-startup.py** - Core initialization
2. **01-parameters.py** - Beamline parameters
3. **10-optics.py** - Optical elements (monochromator, mirrors)
4. **13-mll.py** - Multi-Layer Lens system
5. **15-zp.py** - Zone Plate optics
6. **20-detectors.py** - Detector configuration
7. **50-scans.py** - Scan utilities
8. **72-json-plans.py** - JSON plan execution
9. **90-alignment.py** - Alignment procedures

### Example: Running a Tomography Plan

```python
# Define parameters in JSON file (see startup/72-json-plans.py for format)
# Then execute:
RE(run_tomo_json('/path/to/tomo_params.json'))
```

## Directory Structure

Key startup files (loaded in order):

```
startup/
├── 00-startup.py              # Core initialization, imports, logging
├── 01-parameters.py           # Beamline configuration, constants, metadata
├── 02-tiled_writer.py         # Data writing, Tiled catalog, AD integration
├── 03-user_setup.py           # User environment, proposals, home directory
├── 10-optics.py               # Monochromator, mirrors, energy control, e device
├── 11-machine.py              # Machine parameters, beamline state
├── 12-endstation.py           # Sample stages, aligners, diffraction geometry
├── 13-mll.py                  # Multi-Layer Lens motors, vz, hz, vchi
├── 14-dutycycling.py          # Duty cycle management, beam conditioning
├── 15-zp.py                   # Zone Plate motors, zp alignment, focus
├── 16-ptycho.py               # Ptychography setup, scanning probe
├── 20-detectors.py            # Detector definitions, baseline hardware
├── 20-eiger.py                # Eiger detector, AD plugin, HDF5 output
├── 21-xspress3.py             # XSpressMC detector, fluorescence, energy resolution
├── 22-scalers.py              # Scaler configuration, pulse counting, I0/It
├── 23-interferometers.py      # Interferometer setup, laser feedback
├── 43-bechxn.py               # BECHXN detector, back-scatter analysis
├── 50-scans.py                # Scan utilities, step scans, normalization
├── 51-flyscan.py              # Flyscan implementation, continuous motion
├── 52-suspenders.py           # Scan suspenders, safety, beam monitoring
├── 53-viewer2d.py             # 2D visualization, live display, matplotlib
├── 61-nano-es.py              # Energy scanning, edge mapping, spectroscopy
├── 62-nano-liveplots.py       # Live plotting, real-time visualization
├── 63-nano-monitor-patch.py   # Monitor patches, callback functions
├── 64-nano-merlin.py          # Merlin detector, photon counting
├── 65-nano-eiger*.py          # Eiger detector variants, frame grabbing, data pipeline
├── 66-nano-zebra.py           # Zebra timing control, synchronization
├── 67-nano-flyscans.py        # Flyscan utilities, trajectory generation
├── 67-nano-panda*.py          # PandA I/O controller, real-time triggering
├── 68-nano-functions*.py      # Nanoprobe utility functions, analysis
├── 69-nano-fly2d.py           # 2D flyscan, raster scanning
├── 70-users_plans.py          # User-defined plans, custom workflows
├── 71-energy.py               # Energy control utilities, calibration
├── 72-json-plans.py           # JSON-based plan execution, tomo/diffraction
├── 73-qserver-plans.py        # Qserver integration, remote execution
├── 73-xanes_plans.py          # XANES scan plans, absorption spectroscopy
├── 74-qserver_automap_plans.py # Automated mapping plans, grid scanning
├── 74-zp-tomo-json.py         # ZP tomography from JSON, zone plate reconstruction
├── 76-spectro-tomo.py         # Spectroscopy tomography, energy-resolved imaging
├── 80-uscan.py                # User scan utilities, scripting, batch operations
├── 82-createPDF.py            # PDF generation, report creation
├── 90-alignment.py            # Alignment procedures, focus optimization, rotation
├── 91-timescan.py             # Time-dependent scans, temporal resolution
├── 92-eiger2.py               # Eiger2 detector, next-gen detector control
├── 93-scanplan-panda.py       # PandA scan plans, hardware-synchronized acquisition
├── 94-user-panda-utils.py     # PandA utilities, configuration, triggering
└── 95-xspress3-panda.py       # Xspress3 with PandA, detector triggering
```

## License

See [LICENSE](LICENSE) file.

