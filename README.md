# Workflow scripts for PP, PN, and NN Bidentate Ligand Calculations

Scripts constructed for bisphosphine ligand libray construction (with Jamie Cadge), modified for PN, NN ligands

Series of Python scripts to run the conformational search and generate for bisphosphine ligands. Complexes currently supported:
- Palladium(II) dichloride complexes, (D^D)PdCl2
- Zinc(II) dichloride complexes, (D^D)ZnCl2

## 1. Setup

### 1.1. Python environment

**Using conda:**

Create a new conda environment using the environment.yml file:

`conda env create -n bidentate_calculations --file environment.yml`

**Using pip:**

Clone the repository by running:

`git clone  .....`

Then, install the required packages:

`pip install -r requirements.txt`

### 1.2 Other system requirements

- [CREST 2.12](https://crest-lab.github.io/crest-docs/page/installation/install_basic.html) (Need to specify path or module to CREST executable in the settings/crest_settings.yml).
- [xTB 6.4.0](https://xtb-docs.readthedocs.io/en/latest/setup.html#setup-and-installation) (Need to specify path or module to xTB executable in crest_settings.yml).
- [OpenBabel 2.4.1 or higher](https://anaconda.org/conda-forge/openbabel)
  (To use as part of the conda environment, install separately to the environment using `conda install -c conda-forge openbabel`).

### 1.3 PP_conformer settings

All job submission and directory path settings are stored in the settings/crest_settings.yml file and should be adjusted 
based on system settings:

- **crest_executable**: Location or module containing the CREST executable
- **xtb_executable**: Location or module containing the xTB executable
- **crest_scratch_dir**: Filepath to  the system scratch directory (where files will be moved to for calculation)
- **crest_output_subdir**: Sub-directory for the CREST output
- **crest_results_dir**: Directory where the CREST outputs will be moved to after calculation [do not change]
- **dft_inputs_dir**: Directory for the generated Gaussian16 optimization .com files [do not change]
- **dft_sp_inputs_dir**: Directory for the generated Gaussian16 single point energy .com files [do not change]
- **slurm_partition**: Where the HPC partition is specified for submission *via* SLURM task scheduler
- **slurm_account**: SLURM task scheduler account name
- **slurm_walltime**: Max wall time for the conformer generation
- **slurm_nproc**: Number of processors to use for the conformer generation

*NB: SLURM settings specified in the PP_conformer_settings.yml file are for the conformational searches only. DFT jobs 
should be run using separate submission scripts for your HPC cluster.*

## 2. Find conformers

Run using the find_conformers.py script. To find conformers for all .xyz files in the main directory, run:

`python find_conformers.py` or `python find_conformers.py --all` or `python find_conformers.py -a`

To find conformers for an individual .xyz file, run:

`python find_conformers.py <filename>.xyz`

For help, type: `python find_conformers.py --help` or `python find_conformers.py -h`

This script runs the CREST conformational search at the GFN2-xTB//GFN-FF level. Constraints at the metal center are 
automatically applied. Currently supported metal centers are PdCl2 and ZnCl2.

## 3. Analyze conformers

This script analyzes the conformers generated from the find_conformers.py script:

- Checks for normal termination of the CREST program. Also checks for any structural changes using the before and after 
molecular formulae as well as the D–D (donor atom-donor atom) bond distance.
- Splits the clustered conformer ensemble into individual files, where the suffix "_1" is given to the lowest energy 
conformer 
and so on. 
- Selects conformers based on D-M–D (donor atom-metal–donor atom) bite angle (y-equidistant) and the lowest energy GFN2-xTB conformer. 
Gaussian16 optimization .com input files are generated for the selected conformers (PBE-D3(BJ)/def2SVP level of theory).

Conformers may be analyzed individually by specifying the Ligand ID as an argument:

`python analyze_conformers.py <Ligand ID>`

Alternatively, a list of Ligand IDs may be specified in a text file (crest_submissions.txt) and analyzed in batch:

`python analyze_conformers.py --all` This is the default option.

For help, type `python analyze_conformers.py --help` or `python analyze_conformers.py -h`

Generated Gaussian16 .com files should be run using a separate submission script.

> **Note:**
> D-D bond distances and D-M-D are determined using MORFEUS (for details, see [here](https://digital-chemistry-laboratory.github.io/morfeus/)). Conformers from the CREST output 
> .xyz file are split using OpenBabel 2.4.1 (for details, see "Other system requirements" section).

## 4. Generate Gaussian16 single point energy inputs

Run using the g16_single_points.py script. To generate Gaussian16 single point energy files (PBE0-D3(BJ)/def2TZVP level)
from all finished optimization jobs in the dft_inputs folder, run:

`python g16_single_points.py --generate`

For help, type: `python g16_single_points.py --help` or `python g16_single_points.py -h`

This script is run on all Gaussian16 outputs in the dft_inputs directory. New single point energy input files will be 
placed in the dft_sp_inputs directory:

- Gaussian16 optimization outputs are processed and checked for normal termination and imaginary frequencies.
Files with these errors are placed in separate directories.
- For files with normal termination and no imaginary frequencies, Gaussian16 single point energy calculation input files
are created for metal complexes. Free ligand calculations are generated from the PdCl2 complexes only.

The Gaussian16 outputs from single point energy calculations can also be processed in the g16_single_points.py script.
To do this, run:

`python g16_single_points.py --process`

The outputs are processed in a manner similar to the Gaussian16 DFT optimization jobs.

## Acknowledgments

We would like to acknowledge the following for providing elements of code which were adapted for this project:
- Jordan Dotson
- Lucy van Dijk
- Brittany Haas
- Melissa Hardy
- Lucas Karas
