#%% 
import os
import json
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path
from scipy.spatial import cKDTree
import argparse

#%% 
class IEEGRecon:
    def __init__(self, pre_implant_mri, post_implant_ct, ct_electrodes, output_dir, config_path=None):
        """
        Initialize IEEGRecon with required paths
        
        Args:
            pre_implant_mri (str): Path to pre-implant MRI
            post_implant_ct (str): Path to post-implant CT
            ct_electrodes (str): Path to electrode coordinates CSV
            output_dir (str): Output directory path (required)
            config_path (str, optional): Path to config JSON file
        """
        # Set main properties
        self.preImplantMRI = pre_implant_mri
        self.postImplantCT = post_implant_ct
        self.postImplantCT_electrodes = ct_electrodes
        self.output = Path(output_dir)  # Convert to Path object
        
        # Setup environment variables and paths
        self._setup_environment(config_path)

    def _setup_environment(self, config_path=None):
        """Setup environment variables and paths from config file"""
        if config_path is None:
            config_path = Path(__file__).parent / 'config_ieeg_recon.json'
        
        try:
            with open(config_path) as f:
                config = json.load(f)
            
            # Set paths
            self.fslLoc = config['fsl']['dir']
            self.itksnap = config['itksnap']['dir']
            self.freeSurfer = config['freesurfer']['home']
            self.freeSurferDir = config['freesurfer']['subjects_dir']
            
            # Set environment variables
            env_vars = {
                'FSLDIR': config['fsl']['dir'],
                'FSLOUTPUTTYPE': config['fsl']['output_type'],
                'FREESURFER_HOME': config['freesurfer']['home'],
                'SUBJECTS_DIR': config['freesurfer']['subjects_dir']
            }
            
            for key, value in env_vars.items():
                os.environ[key] = value
            
            # Run FreeSurfer setup if available
            freesurfer_setup = Path(self.freeSurfer) / 'SetUpFreeSurfer.sh'
            if freesurfer_setup.exists():
                subprocess.run(['sh', str(freesurfer_setup)], check=True)
                
        except (FileNotFoundError, KeyError) as e:
            print(f"Warning: Configuration error - {str(e)}")
            print("Using default None values for paths")
            self.fslLoc = self.itksnap = self.freeSurfer = self.freeSurferDir = None

    def module1(self):
        """
        Module1: exports electrode coordinates of post implant CT in voxel and
        native space. Outputs of this module goes in output:ieeg_recon/module1 folder
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module1'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export electrode coordinates in CT space in mm and vox
        elecCTvox = pd.read_csv(self.postImplantCT_electrodes,  sep=r'\s+', header=None)
            
        # Write electrode names
        with open(output_dir / 'electrode_names.txt', 'w') as f:
            f.write('\n'.join(elecCTvox.iloc[:, 0].tolist()))
        
        # Write electrode coordinates
        np.savetxt(
            output_dir / 'electrodes_inCTvox.txt',
            elecCTvox.iloc[:, 1:4].values.astype(int),  # Convert to integers
            delimiter=' ',
            fmt='%d'  # Use integer format
        )

        # Load CT data
        ct_img = nib.load(self.postImplantCT)
        ct_data = ct_img.get_fdata()
        ct_affine = ct_img.affine

        # Convert electrode coordinates from voxel to world space
        elec_vox = elecCTvox.iloc[:, 1:4].values
        elec_homog = np.hstack((elec_vox, np.ones((elec_vox.shape[0], 1))))  # Add homogeneous coordinate
        elecCTmm = np.dot(ct_affine, elec_homog.T).T[:, :3]  # Transform and remove homogeneous coordinate

        # Save world coordinates
        np.savetxt(
            output_dir / 'electrodes_inCTmm.txt',
            elecCTmm,
            delimiter=' ',
            fmt='%.2f'  # Use float format
        )

    def module2(self, reg_type, skip_existing=False):
        """
        Module2: Outputs go in output:ieeg_recon/module2 folder
        
        Args:
            reg_type (str): Registration type - 'gc', 'g', or 'gc_noCTthereshold'
            skip_existing (bool): If True, skip processing if output files exist
        
        Returns:
            dict: Paths to output files
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module2'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Define output file locations
        file_locations = {
            'ct_to_mri': output_dir / 'ct_to_mri.nii.gz',
            'electrodes_inMRI': output_dir / 'electrodes_inMRI.nii.gz',
            'electrodes_inMRI_freesurferLUT': output_dir / 'electrodes_inMRI_freesurferLUT.txt',
            'electrodes_inMRImm': output_dir / 'electrodes_inMRImm.txt',
            'electrodes_inMRIvox': output_dir / 'electrodes_inMRIvox.txt'
        }

        # Check if files exist and skip if requested
        if skip_existing and all(path.exists() for path in file_locations.values()):
            return file_locations

        # Remove negative values from CT image
        cmd = ['fslmaths', 
               self.postImplantCT, 
               '-thr', '0', 
               str(output_dir / 'ct_thresholded.nii.gz')]
        subprocess.run(cmd, check=True)

        # Handle different registration types
        if reg_type == 'gc':
            self._run_greedy_centered_registration(output_dir)
        elif reg_type == 'g':
            self._run_flirt_greedy_registration(output_dir)
        elif reg_type == 'gc_noCTthereshold':
            self._run_greedy_centered_no_threshold(output_dir)
        else:
            raise ValueError("Registration type must be 'gc', 'g', or 'gc_noCTthereshold'")

        # Apply registration to electrode coordinates
        self._transform_electrode_coordinates(output_dir)
        
        # Create electrode spheres
        self._create_electrode_spheres(output_dir)

        return file_locations

    def _run_greedy_centered_registration(self, output_dir):
        """Run greedy registration with image centering"""
        print('Running greedy registration with image centering')
        
        # Run greedy registration
        subprocess.run([
            f"{self.itksnap}/greedy",
            "-d", "3",
            "-i", self.preImplantMRI,
            str(output_dir / 'ct_thresholded.nii.gz'),
            "-o", str(output_dir / 'ct_to_mri.mat'),
            "-a", "-dof", "6",
            "-m", "NMI",
            "-ia-image-centers",
            "-n", "100x100x0x0",
            "-jitter", "0.9",
            "-search", "1000", "10", "20"
        ], stdout=open(output_dir / 'greedy.log', 'w'), check=True)

        # Convert transform to FSL format
        subprocess.run([
            f"{self.itksnap}/c3d_affine_tool",
            "-ref", self.preImplantMRI,
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            str(output_dir / 'ct_to_mri.mat'),
            "-ras2fsl",
            "-o", str(output_dir / 'ct_to_mri_xform.txt')
        ], check=True)

        # Remove temporary mat file
        (output_dir / 'ct_to_mri.mat').unlink()

        # Apply transform
        cmd = ['flirt',
               '-in', str(output_dir / 'ct_thresholded.nii.gz'),
               '-ref', self.preImplantMRI,
               '-init', str(output_dir / 'ct_to_mri_xform.txt'),
               '-out', str(output_dir / 'ct_to_mri.nii.gz'),
               '-applyxfm']
        subprocess.run(cmd, check=True)

    def _transform_electrode_coordinates(self, output_dir):
        """Apply registration transform to electrode coordinates"""
        # Transform mm coordinates
        subprocess.run(["img2imgcoord",
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            "-dest", str(output_dir / 'ct_to_mri.nii.gz'),
            "-xfm", str(output_dir / 'ct_to_mri_xform.txt'),
            "-mm", str(Path(self.output) / 'ieeg_recon/module1/electrodes_inCTmm.txt')
        ], stdout=open(output_dir / 'electrodes_inMRImm.txt', 'w'), check=True)

        # Transform voxel coordinates
        subprocess.run(["img2imgcoord",
            "-src", str(output_dir / 'ct_thresholded.nii.gz'),
            "-dest", str(output_dir / 'ct_to_mri.nii.gz'),
            "-xfm", str(output_dir / 'ct_to_mri_xform.txt'),
            "-vox", str(Path(self.output) / 'ieeg_recon/module1/electrodes_inCTvox.txt')
        ], stdout=open(output_dir / 'electrodes_inMRIvox.txt', 'w'), check=True)

    def _create_electrode_spheres(self, output_dir):
        """Create spheres for electrodes in registered space"""
        # Load registered CT data
        ct_img = nib.load(output_dir / 'ct_to_mri.nii.gz')
        ct_data = ct_img.get_fdata()
        ct_affine = ct_img.affine

        # Create blank image
        blank_data = np.zeros_like(ct_data)
        vox_coords = np.array(np.where(blank_data == 0)).T
        
        # Convert to world coordinates
        vox_homog = np.hstack((vox_coords, np.ones((vox_coords.shape[0], 1))))
        world_coords = np.dot(ct_affine, vox_homog.T).T[:, :3]

        # Load electrode coordinates
        electrodes_mm = np.loadtxt(output_dir / 'electrodes_inMRImm.txt', skiprows=1)
        electrode_names = np.loadtxt(
            Path(self.output) / 'ieeg_recon/module1/electrode_names.txt',
            dtype=str
        )

        # Create FreeSurfer LUT
        n_electrodes = len(electrode_names)
        lut_data = {
            'index': np.arange(1, n_electrodes + 1),
            'names': electrode_names,
            'R': np.full(n_electrodes, 90),
            'G': np.full(n_electrodes, 150),
            'B': np.full(n_electrodes, 60),
            'alpha': np.zeros(n_electrodes)
        }
        pd.DataFrame(lut_data).to_csv(
            output_dir / 'electrodes_inMRI_freesurferLUT.txt',
            sep=' ',
            header=False,
            index=False
        )

        # Find points within 2mm of each electrode
        tree = cKDTree(world_coords)
        dist, idx = tree.query(electrodes_mm, k=1)
        
        # Create electrode map
        electrode_data = blank_data.copy()
        mask = dist <= 2
        for i, (valid, coord) in enumerate(zip(mask, vox_coords[idx]), 1):
            if valid:
                electrode_data[tuple(coord)] = i

        # Save electrode map
        nib.save(
            nib.Nifti1Image(electrode_data, ct_affine),
            output_dir / 'electrodes_inMRI.nii.gz'
        )

    def module2_QualityAssurance(self, file_locations, imageviewer):
        """
        Generate quality assurance visualizations for module2 results
        
        Args:
            file_locations (dict): Dictionary containing paths to module2 output files
            imageviewer (str): Type of viewer to use - 'freeview_snapshot', 'freeview', or 'itksnap'
        """
        # Create output directory
        output_dir = Path(self.output) / 'ieeg_recon' / 'module2'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Only launch viewer if not set to 'none'
        if imageviewer.lower() != 'none':
            if imageviewer == 'freeview':
                # Open interactive freeview session
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "sagittal"
                ], check=True)
            elif imageviewer == 'freeview_snapshot':
                # Generate sagittal view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "sagittal",
                    "-ss", str(output_dir / "QA_registation_sagittal.png")
                ], check=True)

                # Generate coronal view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "coronal",
                    "-ss", str(output_dir / "QA_registation_coronal.png")
                ], check=True)

                # Generate axial view
                subprocess.run([
                    "freeview",
                    "-v", self.preImplantMRI,
                    f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}",
                    "-viewport", "axial",
                    "-ss", str(output_dir / "QA_registation_axial.png")
                ], check=True)

                # Generate 3D view
                subprocess.run([
                    "freeview",
                    "-v", f"{file_locations['ct_to_mri']}:colormap=heat",
                    f"{file_locations['electrodes_inMRI']}:colormap=lut:lut={file_locations['electrodes_inMRI_freesurferLUT']}:isosurface=on",
                    "-viewport", "3d", "-view", "anterior",
                    "-ss", str(output_dir / "QA_registation_3D.png")
                ], check=True)
            elif imageviewer == 'itksnap':
                # Open interactive ITK-SNAP session
                subprocess.run([
                    f"{self.itksnap}/itksnap",
                    "-g", self.preImplantMRI,
                    "-o", file_locations['ct_to_mri']
                ], check=True)
            else:
                raise ValueError(f"Unknown imageviewer option: {imageviewer}")

#%% 
if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='iEEG Reconstruction Pipeline')
    parser.add_argument('filepath', type=str, help='Path to the file_paths.csv')
    parser.add_argument('index', type=int, help='Index of the subject to process')
    parser.add_argument('--config', type=str, default=None, 
                       help='Path to config file (default: config_ieeg_recon.json in script directory)')
    parser.add_argument('--modules', type=str, nargs='+', default=['1', '2'],
                       choices=['1', '2'], help='Modules to run (default: 1 2)')
    parser.add_argument('--skip-existing', action='store_true',
                       help='Skip processing if output files exist')
    parser.add_argument('--reg-type', type=str, default='gc', choices=['gc', 'g', 'gc_noCTthereshold'],
                       help='Registration type (default: gc)')
    parser.add_argument('--qa-viewer', type=str, default='freeview',
                       choices=['freeview', 'freeview_snapshot', 'itksnap', 'none'],
                       help='Quality assurance viewer (default: freeview)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Read filepath CSV
    filepath = pd.read_csv(args.filepath)
    
    # Get subject data using index
    if args.index >= len(filepath):
        raise ValueError(f"Index {args.index} is out of range. File contains {len(filepath)} subjects.")
    
    # Set paths for the selected subject
    pre_implant_mri = filepath['t1w'][args.index]
    post_implant_ct = filepath['ct'][args.index]
    ct_electrodes = filepath['electrodes'][args.index]
    output_dir = Path(filepath['sub'][args.index])
    
    # Set config path
    if args.config is None:
        config_path = Path(__file__).parent / 'config_ieeg_recon.json'
    else:
        config_path = Path(args.config)
    
    # Initialize reconstruction object
    print(f"Processing subject at index {args.index}")
    recon = IEEGRecon(
        pre_implant_mri=pre_implant_mri,
        post_implant_ct=post_implant_ct,
        ct_electrodes=ct_electrodes,
        output_dir=output_dir,
        config_path=config_path
    )
    
    # Run selected modules
    if '1' in args.modules:
        print("Running Module 1...")
        recon.module1()
    
    file_locations = None
    if '2' in args.modules:
        print("Running Module 2...")
        file_locations = recon.module2(args.reg_type, skip_existing=args.skip_existing)
        
        print("Output files:")
        for name, path in file_locations.items():
            print(f"{name}: {path}")
        
        recon.module2_QualityAssurance(file_locations, args.qa_viewer)
    
    print("Processing complete!")

# %%
