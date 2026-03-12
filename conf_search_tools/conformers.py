"""
GENERATION 2 BISPHOSPHINE LIBRARY WORKFLOW
-Conformer Functions Script-
"""

import os
import yaml
import shutil
import subprocess
import time
from morfeus.io import read_xyz
import numpy as np
from conf_search_tools import utils
import pandas as pd
from pathlib import Path

# define settings variables from crest_settings.yml file
with open("./settings/crest_settings.yml") as f:
    settings = yaml.safe_load(f)

CREST_METHOD = "gfn2//gfnff"
CREST_ENERGY_WINDOW = 5
CONSTRAINT_FORCE_CONSTANT = 1.0

# Donor atom reference distances based on value from J. Chem. Soc., Perkin Trans. 2, 1987, S1-S9
REFERENCE_DISTANCES = {
    ('P', 'P'): 2.3,
    ('P', 'N'): 1.7,
    ('N', 'P'): 1.7,
    ('N', 'N'): 1.3  # N-N bond length dependent on torsion angle, this value is approximated from P,P
    # and P,N bond lengths (some jobs reasonably produce bond lengths shorter than this.
}


def write_constraints_file(xyz_file: str, force_constant: float = CONSTRAINT_FORCE_CONSTANT):
    """
    Function to write a constraints file for CREST conformer generation
    Constraints are set for the bite angle and four angles between the metal center and the ligand atoms
    If the structure contains ferrocene, additional constraints are set for the Fe-C bond lengths
    :param xyz_file: Name of the .xyz file to write constraints for
    :param force_constant: Force constant for the constraints, default is 1.0
    """
    start_line = "$constrain"
    end_line = "$end"

    # determine whether the structure contains ferrocene or not (use MORFEUS read_xyz)
    # set ferrocene = True if Fe is found in the list of elements
    elements, coordinates = read_xyz(xyz_file)
    ferrocene = False
    if 'Fe' in elements:
        ferrocene = True

    # determine atom coordinates for metal center
    (metal, donor1, donor2, donor_elements, metal_ligand_atom1, metal_ligand_atom2,
     metal_ligand_elements) = utils.find_metal_and_donor_atoms(xyz_file)
    bite_angle_line = f"  angle: {donor1}, {metal}, {donor2}, auto"
    angle1_line = f"  angle: {donor1}, {metal}, {metal_ligand_atom1}, auto"
    angle2_line = f"  angle: {donor1}, {metal}, {metal_ligand_atom2}, auto"
    angle3_line = f"  angle: {donor2}, {metal}, {metal_ligand_atom1}, auto"
    angle4_line = f"  angle: {donor2}, {metal}, {metal_ligand_atom2}, auto"

    # create a new .inp file and write constraints to it
    # xyz_filename = xyz_file.name
    try:
        with open(xyz_file.split(".")[0] + ".inp", "w") as file:
            file.write(f"{start_line}\n")
            file.write(f"  force constant={force_constant}\n")
            file.write(f"{bite_angle_line}\n")
            file.write(f"{angle1_line}\n")
            file.write(f"{angle2_line}\n")
            file.write(f"{angle3_line}\n")
            file.write(f"{angle4_line}\n")

            if ferrocene is True:
                fe, carbon_list = utils.find_ferrocene_atoms(xyz_file)
                for carbon in carbon_list:
                    file.write(f"  distance: {fe}, {carbon}, auto\n")
            else:
                pass

            file.write(f"{end_line}")

    except FileNotFoundError:
        raise Exception(f"Error writing constraints file for {xyz_file}")


def setup_slurm_submission():
    """
    Function to set up the CREST slurm submission script
    :return:
    """
    # name of the SLURM submission script
    slurm_filename = "crest_bisphosphine.slurm"

    submit_template = f"""#!/bin/bash
#SBATCH --partition={settings["slurm_partition"]}
#SBATCH --account={settings["slurm_account"]}
#SBATCH --time={settings["slurm_walltime"]}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={settings["slurm_nproc"]}
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N

env | grep SLURM

{settings['xtb_executable']}
{settings['crest_executable']}

file=$1
TMPDIR={settings["crest_scratch_dir"]}
mkdir -p $TMPDIR
WORKDIR=$PWD

cp $WORKDIR/${{file}}.xyz $TMPDIR/
cp $WORKDIR/${{file}}.inp $TMPDIR/
cd $TMPDIR

ulimit -s unlimited
export OMP_STACKSIZE=30GB
export OMP_NUM_THREADS={settings["slurm_nproc"]},1
export MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS={settings["slurm_nproc"]}

## TODO: change level of theory here:

crest ${{file}}.xyz --cinp ${{file}}.inp --{CREST_METHOD} --ewin {CREST_ENERGY_WINDOW} --cluster --noreftopo --keepdir > ${{file}}_crestlog.log
mv crest_conformers.xyz ${{file}}_crest_conformers.xyz
mv crest_clustered.xyz ${{file}}_crest_clustered.xyz
mv crest_best.xyz ${{file}}_lec.xyz

mkdir $WORKDIR/${{file}}_output
cp -r * $WORKDIR/${{file}}_output
rm -rf $TMPDIR
        """

    # check if the SLURM submission script already exists
    if not os.path.isfile(slurm_filename):
        # create the file
        with open(slurm_filename, "w") as file:
            file.write(submit_template)
        print(f"{slurm_filename} submission script created")
    else:
        print(f"{slurm_filename} already exists")


def execute_crest():
    """
    Function to execute the CREST conformational search using the SLURM submission script
    :return:
    """
    setup_slurm_submission()
    command = "bash conf_search_tools/run_crest.sh"
    subprocess.run(command, shell=True)
    # remove the SLURM submission script after execution
    os.remove("./crest_submission.slurm")


def copy_crest_files(xyz_file: str, results_dir: Path):
    """
    Function to copy CREST conformer files to the results directory
    Conformer files include:
    (i) All conformers generated by CREST
    (ii) Clustered conformers
    (iii) Lowest energy conformer
    :param xyz_file: Name of the .xyz file to copy CREST conformers for
    :param results_dir: Path to the results directory
    """
    xyz_filename = xyz_file[:-4]
    crest_path = Path(f"{xyz_filename}_output")
    all_conformers_path = crest_path / f"{xyz_filename}_crest_conformers.xyz"
    clustered_path = crest_path / f"{xyz_filename}_crest_clustered.xyz"
    lec_path = crest_path / f"{xyz_filename}_lec.xyz"
    os.makedirs(results_dir, exist_ok=True)

    try:
        shutil.copy(all_conformers_path, results_dir)
    except FileNotFoundError:
        print(f"ERROR! {xyz_filename} CREST conformers file not found")

    try:
        shutil.copy(clustered_path, results_dir)
    except FileNotFoundError:
        print(f"ERROR! {xyz_filename} clustered CREST conformers file not found")

    try:
        shutil.copy(lec_path, results_dir)
    except FileNotFoundError:
        print(f"ERROR {xyz_filename} CREST lowest energy conformer file not found")
    

def split_conformers(xyz_file: str, results_dir: Path):
    """
    Function to split up all CREST conformer .xyz files (clustered only)
    in the crest results directory into individual .xyz files.
    Requires OpenBabel 2.4.1 to be loaded
    :param xyz_file: Name of the .xyz file to split conformers for
    :param results_dir: Path to the results directory
    """
    xyz_filename = xyz_file[:-4]
    clustered_filepath = results_dir / f"{xyz_filename}_crest_clustered.xyz"
    conformer_dir = results_dir / "split_conformers"
    os.makedirs(conformer_dir, exist_ok=True)

    try:
        shutil.copy(clustered_filepath, conformer_dir)
    except FileNotFoundError:
        print(
            f"ERROR! Copying {xyz_filename} clustered conformers to split_conformers directory failed - "
            f"file not found")

    new_clustered_filepath = conformer_dir / clustered_filepath.name

    obabel_command = f"obabel {new_clustered_filepath} -O {conformer_dir / xyz_filename}_xtb_.xyz -m"

    print(f"Splitting clustered conformers using OpenBabel for {xyz_filename}")
    subprocess.run(obabel_command, shell=True)

    os.remove(new_clustered_filepath)
    files = os.listdir(conformer_dir)
    xyz_files = [file for file in files if xyz_filename in file]
    for file in xyz_files:
        conformer_number = file.split('_')[-1].split('.')[0]
        conformer_number_zfilled = conformer_number.zfill(4)
        new_filename = f"{xyz_filename}_xtb_{conformer_number_zfilled}.xyz"
        os.rename(conformer_dir / file, conformer_dir / new_filename)


def compress_crest_output(xyz_file: str):
    """
    Function to compress the CREST output folder for archiving
    Then removes the CREST output folder
    :param xyz_file: Name of the .xyz file to compress CREST output for
    """
    xyz_filename = xyz_file[:-4]
    output_path = Path(f"{xyz_filename}_output")
    archived_path = f"{xyz_filename}_output"

    shutil.make_archive(archived_path, 'zip', output_path)
    print(f"Output for {xyz_filename} compressed to {archived_path}")
    time.sleep(20)  # add a wait of 20 seconds to ensure the zip file is created before removing the directory

    remove_command = f"rm -rf {output_path}"
    subprocess.run(remove_command, shell=True)



def remove_files(xyz_file: str):
    """
    Function to remove unecessary files generated from CREST
    :param xyz_file: xyz file name
    """
    xyz_filename = xyz_file[:-4]

    # remove .slurm files
    slurm_file = f"{xyz_filename}_crest.slurm"
    if os.path.exists(slurm_file):
        os.remove(slurm_file)
    else:
        print(f"Cannot remove {slurm_file}, does not exist!")

    # remove .inp files
    inp_file = f"{xyz_filename}.inp"
    if os.path.exists(inp_file):
        os.remove(inp_file)
    else:
        print(f"Cannot remove {inp_file}, does not exist!")


def select_lec(df: pd.DataFrame, column: str) -> str:
    """
    Returns the lowest energy conformer from a pre-generated Pandas dataframe
    :param df: DataFrame containing conformer data
    :param column: Column heading containing the energies
    :return: Lowest energy conformer
    """
    lec_df = df.sort_values(column)
    indices = np.linspace(0, len(lec_df) - 1, 1, dtype=int)
    lec = lec_df.iloc[indices]

    return lec


def select_equidistant_values(df, column: str, y: int = 10) -> object:
    """
    Selects conformers based on a y-equidistant distribution of values from a pre-generated Pandas dataframe
    :param df: DataFrame containing conformer data
    :param column: Column heading containing the feature to select conformers by
    :param y: Number of conformers to select, default is 10
    :return: Selected conformers
    """
    sorted_df = df.sort_values(column)
    # print(sorted_df)  # TODO: do we want an option here which prints a csv with the conformer details
    indices = np.linspace(0, len(sorted_df) - 1, y, dtype=int)
    selected_values = sorted_df.iloc[indices]

    return selected_values


def crest_end(xyz_file: str, results_dir: Path):
    """
    This function runs after a completed CREST conformer search

    The following checks are performed:
    (i) Normal termination of CREST
    (ii) Molecular formula check between input and output files
    # (iii) Donor-atom distance check

    Then the following actions are performed:
    (i) Copy conformer files to the results directory
    (ii) Split clustered conformers into individual XYZ coordinate files
    (iii) Compress CREST output folder for archiving
    (iv) Remove unecessary files generated for CREST

    :param xyz_file: Name of the .xyz file to check CREST output for
    :param results_dir: Path to the results directory
    """
    xyz_filename = xyz_file[:-4]
    subdirectory = Path(f"{xyz_filename}_output")


    # log_file = f"{xyz_filename}_crestlog.log"
    # log_filepath = subdirectory / log_file

    # with open(log_filepath, "r") as file:
    #     line_list = list(file)
    #     line_list.reverse()

        # for line in line_list:
        #     normal_term = " CREST terminated normally."
        #     if line.find(normal_term) != 1:
        #         print(f"CREST conformational search for {xyz_file} terminated normally.")
    copy_crest_files(xyz_file, results_dir)
    split_conformers(xyz_file, results_dir)
    compress_crest_output(xyz_file)
    remove_files(xyz_file)
        #         break
        #     else:
        #         print(f"Error in termination of CREST log file for {xyz_file}, check!")
        #         break

    # results_directory = subdirectory
    lec_file = f"{xyz_filename}_lec.xyz"
    lec_filepath = results_dir / lec_file

    input_molecular_formula = utils.get_molecular_formula(xyz_file)
    lec_output_molecular_formula = utils.get_molecular_formula(lec_filepath)

    if input_molecular_formula == lec_output_molecular_formula:
        print(f"Input and output molecular formulae match for LEC!")
    else:
        print(f"Molecular formulae for input {input_molecular_formula} and output {lec_output_molecular_formula} "
              f"structures do not match, check!")

    split_conf_directory = results_dir / "split_conformers"
    files = [file for file in split_conf_directory.iterdir() if file.is_file()]
    split_xyz_files = [file for file in files if file.suffix == '.xyz']

    # for file in split_xyz_files:
    #     if xyz_filename in file.name:
    #         donor_distance, donor_elements = utils.get_donor_atom_distance(file)
    #         donor_atom1, donor_atom2 = donor_elements
    #         reference_distance = REFERENCE_DISTANCES.get(donor_elements, None)
    #         if donor_distance > reference_distance:
    #             print(f"{donor_atom1}-{donor_atom2} distance in {file} ok.")
    #         else:
    #             print(f"WARNING! Short {donor_atom1}–{donor_atom2} distance found for {file} "
    #                   f"({donor_distance:.2f} angstroms), check!")


# def select_files_for_dft(xyz_file_list: list) -> list:
#     """
#     Function to select conformers for DFT calculation
#     :param xyz_file_list: List of .xyz files to select conformers from
#     :return: List of selected conformers for DFT calculation
#     """

#     xyz_files_for_dft = []
#     # takes all ligand conformers for one ID and adds them to a new pandas dataframe
#     conformers = [file for file in xyz_file_list]
#     df = pd.DataFrame({'xyz_file': conformers})
#     # extracts xTB energies from XYZ coordinate files and adds them to the dataframe
#     energies = [utils.get_xtb_energy(conformer) for conformer in conformers]
#     df['xtb_energy'] = energies
#     # # extracts ligand bite angle using MORFEUS and adds them to the dataframe
#     bite_angles = [utils.get_bite_angle(conformer) for conformer in conformers]
#     df['bite_angle'] = bite_angles
#     # selects lowest energy conformer
#     lec_selected = select_lec(df, 'xtb_energy')
#     lec_selected_list = lec_selected["xyz_file"].tolist()
#     # selects equidistant conformers based on bite angle values
#     eq_selected = select_equidistant_values(df, 'xtb_energy')
#     eq_selected_list = eq_selected['xyz_file'].tolist()
#     # combines all selected conformers into one list
#     all_selected = eq_selected_list + lec_selected_list
#     print(all_selected)
#     selected = list(set(all_selected))
#     xyz_files_for_dft.extend(selected)

#     return xyz_files_for_dft


def select_files_for_dft(xyz_file_list: list) -> list:
    """
    Function to select conformers for DFT calculation
    :param xyz_file_list: List of .xyz files to select conformers from
    :return: List of selected conformers for DFT calculation
    """

    xyz_files_for_dft = []
    # takes all ligand conformers for one ID and adds them to a new pandas dataframe
    conformers = [file for file in xyz_file_list]
    df = pd.DataFrame({'xyz_file': conformers})
    # extracts xTB energies from XYZ coordinate files and adds them to the dataframe
    energies = [utils.get_xtb_energy(conformer) for conformer in conformers]
    df['xtb_energy'] = energies
    # extracts ligand bite angle using MORFEUS and adds them to the dataframe
    bite_angles = [utils.get_bite_angle(conformer) for conformer in conformers]
    df['bite_angle'] = bite_angles
    # selects lowest energy conformer
    lec_selected = select_lec(df, 'xtb_energy')
    lec_selected_list = lec_selected["xyz_file"].tolist()
    # selects equidistant conformers based on bite angle values
    eq_selected = select_equidistant_values(df, 'bite_angle')
    eq_selected_list = eq_selected['xyz_file'].tolist()
    # combines all selected conformers into one list
    all_selected = eq_selected_list + lec_selected_list
    selected = list(set(all_selected))
    xyz_files_for_dft.extend(selected)

    return xyz_files_for_dft