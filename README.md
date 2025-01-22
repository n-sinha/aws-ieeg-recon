# iEEG Electrode Reconstruction Pipeline

This pipeline provides tools for reconstructing intracranial EEG (iEEG) electrode locations from pre-implant MRI and post-implant CT scans.

## Prerequisites

- Python 3.8 or higher
- FSL (FMRIB Software Library)
- FreeSurfer
- ITK-SNAP with Greedy registration tool

### System Dependencies

1. FSL: [Installation Guide](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
2. FreeSurfer: [Installation Guide](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
3. ITK-SNAP with Greedy: [Installation Guide](http://www.itksnap.org/pmwiki/pmwiki.php?n=Main.Downloads)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/ieeg-reconstruction.git
cd ieeg-reconstruction
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Configuration

1. Copy the example configuration file:
```bash
cp code/config_ieeg_recon.example.json code/config_ieeg_recon.json
```

2. Edit `code/config_ieeg_recon.json` with your system paths:
```json
{
    "fsl": {
        "dir": "/path/to/fsl",
        "output_type": "NIFTI_GZ"
    },
    "itksnap": {
        "dir": "/path/to/itksnap"
    },
    "freesurfer": {
        "home": "/path/to/freesurfer",
        "subjects_dir": "/path/to/subjects_dir"
    }
}
```

## Data Structure

The pipeline expects data in BIDS format:

```
data/
├── sub-01/
│   ├── ses-preimplant01/
│   │   └── anat/
│   │       └── sub-01_ses-preimplant01_T1w.nii.gz
│   └── ses-implant01/
│       ├── ct/
│       │   └── sub-01_ses-implant01_ct.nii.gz
│       └── ieeg/
│           └── sub-01_ses-implant01_electrodes.txt
└── ...
```

## Quick Start Guide

1. **Prepare Your Data**
   - Organize your data in BIDS format as shown above
   - Ensure your electrode coordinates file is space-separated with format:
     ```
     electrode_name x y z
     ```

2. **Generate File Paths**
   ```bash
   python code/utilities.py
   ```
   This creates `code/file_paths.csv` containing paths to all required files.

3. **Run the Pipeline**
   ```bash
   python code/ieeg_recon.py code/file_paths.csv 0 --qa-viewer freeview
   ```

## Detailed Usage

### Command Line Interface

```bash
python code/ieeg_recon.py <file_paths.csv> <subject_index> [options]
```

#### Required Arguments
- `file_paths.csv`: Path to CSV containing file locations
- `subject_index`: Index of subject to process (0-based)

#### Optional Arguments
```
--config PATH          Path to config file
--modules {1,2}       Modules to run (default: both)
--skip-existing       Skip if output exists
--reg-type {gc,g,gc_noCTthereshold}
                     Registration type (default: gc)
--qa-viewer {freeview,freeview_snapshot,itksnap}
                     QA visualization tool (default: freeview)
```

### Registration Types

1. `gc` (Default): Greedy registration with image centering
   - Best for most cases
   - Uses mutual information metric
   - Automatically centers images before registration

2. `g`: Standard greedy registration
   - Use when `gc` fails
   - No pre-centering of images

3. `gc_noCTthereshold`: No CT thresholding
   - Use when CT values are already clean
   - Skips negative value removal

### Quality Assurance Options

1. `freeview` (Interactive):
   ```bash
   python code/ieeg_recon.py file_paths.csv 0 --qa-viewer freeview
   ```
   - Opens interactive FreeSurfer viewer
   - Best for detailed inspection

2. `freeview_snapshot`:
   ```bash
   python code/ieeg_recon.py file_paths.csv 0 --qa-viewer freeview_snapshot
   ```
   - Generates static images
   - Creates sagittal, coronal, axial, and 3D views

3. `itksnap`:
   ```bash
   python code/ieeg_recon.py file_paths.csv 0 --qa-viewer itksnap
   ```
   - Opens ITK-SNAP viewer
   - Useful for registration quality check

## Pipeline Modules

### Module 1: Coordinate Export
```bash
python code/ieeg_recon.py file_paths.csv 0 --modules 1
```
- Extracts electrode coordinates
- Outputs:
  - `electrode_names.txt`: List of electrode labels
  - `electrodes_inCTvox.txt`: Voxel coordinates
  - `electrodes_inCTmm.txt`: World coordinates (mm)

### Module 2: Registration & Visualization
```bash
python code/ieeg_recon.py file_paths.csv 0 --modules 2
```
- Registers CT to MRI
- Transforms coordinates
- Creates visualization files
- Outputs:
  - `ct_to_mri.nii.gz`: Registered CT
  - `electrodes_inMRI.nii.gz`: Electrode spheres
  - `electrodes_inMRImm.txt`: MRI-space coordinates
  - QA visualizations

## Output Structure

```
output/
└── sub-XX/
    └── ieeg_recon/
        ├── module1/
        │   ├── electrode_names.txt
        │   ├── electrodes_inCTvox.txt
        │   └── electrodes_inCTmm.txt
        └── module2/
            ├── ct_to_mri.nii.gz
            ├── electrodes_inMRI.nii.gz
            ├── electrodes_inMRI_freesurferLUT.txt
            ├── electrodes_inMRImm.txt
            ├── electrodes_inMRIvox.txt
            └── QA/
                ├── QA_registration_sagittal.png
                ├── QA_registration_coronal.png
                ├── QA_registration_axial.png
                └── QA_registration_3D.png
```

## Common Issues and Solutions

1. **Registration Failure**
   - Try different registration types
   - Check input image quality
   - Verify CT threshold values

2. **Missing Dependencies**
   - Ensure FSL, FreeSurfer, and ITK-SNAP are in PATH
   - Verify config.json paths
   - Check environment variables

3. **File Format Issues**
   - Ensure BIDS compliance
   - Check electrode coordinate file format
   - Verify NIfTI file integrity

## Support

For issues and questions:
1. Check the [Common Issues](#common-issues-and-solutions) section
2. Open an issue on GitHub
3. Include relevant error messages and logs

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.