"""
GENERATION 2 BISPHOSPHINE LIBRARY WORKFLOW
-DFT Functions Script-
Parts of this code are adapted from scripts by Jordan Dotson, Lucy van Dijk, Lucas Karas, Brittany Haas and Melissa
Hardy
"""
import os
import re
import yaml
import shutil
from pathlib import Path
from pp_conformers import utils

# Get Gaussian16 DFT settings from the gaussian_settings.yml file
with open("./settings/gaussian_settings.yml", "r") as file:
    settings = yaml.safe_load(file)

# Gaussian16 DFT calculation settings:
OPT_PROCS = settings['optimization_procs']
OPT_MEM = settings['optimization_memory']
SP_PROCS = settings['single_point_procs']
SP_MEM = settings['single_point_memory']

# Gaussian16 DFT method definitions:
OPT_FUNCTIONAL = "PBEPBE"
SP_FUNCTIONAL = "PBE1PBE"
OPT_BASIS_SET = "def2SVP"
SP_BASIS_SET = "def2TZVP"
DISPERSION = "EmpiricalDispersion=GD3BJ"

# heavy atoms included - where there are def2-ECPs available
DEF2_HEAVY_ATOMS = ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
                    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"]

# Gaussian16 DFT property keywords:
NBO = "Pop=NBO7"
HIRSHFELD = "Pop=Hirshfeld"
NMR = "NMR"
VOLUME = "Volume"

# single point job setup:
sp_1_properties = f"{NMR}"
sp_2_metal_properties = f"{NBO} {VOLUME}"
sp_2_nom_properties = f"{NBO} {HIRSHFELD} {VOLUME}"


def write_opt_com_template(xyz_file: Path):
    """
    This function generates the .com templates for the metal complex geometry optimization calculation
    :param xyz_file: Path to the .xyz file
    """
    xyz_filename = xyz_file.name
    cat_name = xyz_filename[:-4]
    cat_name_split = cat_name.split("_")

    opt_com_filename = f"{cat_name_split[0]}_opt_{cat_name_split[2]}.com"
    geometry = utils.geom_from_xyz(xyz_file)

    # make metal optimization .com template
    with open(xyz_file.parent / opt_com_filename, "w") as new_com:
        new_com.write("#Put Keywords Here, check Charge and Multiplicity.\n")
        new_com.write("\n")
        new_com.write(f"{cat_name_split[0]} OPT \n")
        new_com.write("\n")
        new_com.write("0 1 \n")
        for line in geometry:
            new_com.write(line)
        new_com.write("\n")


def write_spc_com_template(logfile: Path):
    """
    Function to generate the .com templates for the metal and free ligand single point energy calculation
    :param logfile: Path to the .log file for the DFT optimization calculation
    """
    log_filename = logfile.name
    cat_name = log_filename[:-4]
    cat_name_split = cat_name.split("_")

    spe_com_filename = f"{cat_name_split[0]}_spc_{cat_name_split[2]}.com"
    spe_nom_com_filename = f"{cat_name_split[0]}_spc-nom_{cat_name_split[2]}.com"

    geometry = utils.geom_from_gaussian(logfile)

    # make metal single point .com template
    with open(logfile.parent / spe_com_filename, 'w') as new_com:
        new_com.write("#Put Keywords Here, check Charge and Multiplicity.\n")
        new_com.write("\n")
        new_com.write(f"{cat_name_split[0]} SPE \n")
        new_com.write("\n")
        new_com.write("0 1 \n")

        for a_num in range(len(geometry)):
            new_com.write(f"{str(geometry[a_num][0])}   {str(geometry[a_num][1])}   {str(geometry[a_num][2])}   "
                          f"{str(geometry[a_num][3])}\n")

    # make no metal (nom) single point, for Pd-containing structures only
    atom_list = []

    for a_num in range(len(geometry)):
        atom_list.append(geometry[a_num][0])

    if 'Pd' in atom_list:
        with open(logfile.parent / spe_nom_com_filename, 'w') as new_com:
            new_com.write("#Put Keywords Here, check Charge and Multiplicity.\n")
            new_com.write("\n")
            new_com.write(f"{cat_name_split[0]} SPE NoM \n")
            new_com.write("\n")
            new_com.write("0 1 \n")
            for a_num in range(len(geometry)):
                if geometry[a_num][0] == "Pd" or geometry[a_num][0] == "Cl":
                    pass
                else:
                    new_com.write(
                       f"{str(geometry[a_num][0])}  {str(geometry[a_num][1])}   {str(geometry[a_num][2])}   "
                       f"{str(geometry[a_num][3])}\n")


def write_g16_opt_com(comfile: Path):
    """
    This function writes the required level of theory and to the G16 .com file template (geom optimization)
    :param comfile: Path to the .com file
    """
    with open(comfile, "r") as file:
        comfile_name = comfile.name
        name_without_extension = comfile_name[:-4]
        lines = file.readlines()
        start_index = 0
        for i, line in enumerate(lines):
            pattern = r"^[a-zA-Z ].*\s\s\s.*\d$"
            match = re.match(pattern, line)
            if match:
                start_index = i-1
                break
        end_index = 0
        for i, line in reversed(list(enumerate(lines))):
            pattern_end = r"^[a-zA-Z ].*\s\s\s.*\d$"
            match = re.match(pattern_end, line)
            if match:
                end_index = i+1
                break
        # removes all lines before the start_index and after the end_index, defined above
        lines = lines[start_index:end_index]

        # creates a list of unique atoms for the gen section
        all_atoms = []
        for index, line in enumerate(lines):
            if index != 0:
                atom = line.split(" ")[0]
                all_atoms.append(atom)

        unique_atoms = list(set(all_atoms))
        remove_item = "\n"

        if remove_item in unique_atoms:
            unique_atoms.remove(remove_item)

        unique_atoms_list = " ".join(unique_atoms)

        # checking for need for ECP
        need_ecp = []
        for line in lines:
            for i in DEF2_HEAVY_ATOMS:
                if i in line:
                    need_ecp.append(i)

    # write the optimization job file
    with open(comfile, "w") as file:
        file.write(f"%NProcShared={OPT_PROCS}\n%Mem={OPT_MEM}\n")
        file.write(f"%chk={name_without_extension}.chk\n")
        file.write(f"#N OPT {DISPERSION} {OPT_FUNCTIONAL}/gen pseudo=read FREQ=NORAMAN Integral(Ultrafine)")
        file.write(f"\n\n{name_without_extension}\n\n")

        file.writelines(lines)
        file.write("\n")

        file.write(f"-{unique_atoms_list} 0\n{OPT_BASIS_SET}\n****\n\n")
        if need_ecp:
            for atom in need_ecp:
                file.write(f"-{atom} 0\n{OPT_BASIS_SET} \n")
            file.write("\n\n\n")
        print(f"Sucessfully converted {name_without_extension} to OPT G16 com file!")


def write_g16_spc_com(comfile: Path):
    """
    This function writes the required level of theory and properties to the G16 .com file template (single point
    energy calculation)
    :param comfile: Path to the .com file
    """
    with open(comfile, "r") as file:
        lines = file.readlines()

    # Remove blank lines from the end of the file
    while lines and lines[-1].strip() == "":
        lines.pop()

    filename = comfile.name
    name_without_extension = filename[:-4]

    charge_spin = str(lines[4])
    lines = lines[5:]

    # creates a list of unique atoms for the gen section
    all_atoms = []
    for index, line in enumerate(lines):
        atom = line.split(" ")[0]
        all_atoms.append(atom)

    unique_atoms = list(set(all_atoms))
    remove_item = "\n"

    if remove_item in unique_atoms:
        unique_atoms.remove(remove_item)

    unique_atoms_list = " ".join(unique_atoms)

    # checking for need for ECP
    need_ecp = []
    for line in lines:
        for i in DEF2_HEAVY_ATOMS:
            if i in line:
                need_ecp.append(i)

    # write the file (2 single point jobs)
    with open(comfile, "w") as file:
        # start of first job, first single point
        file.write(f"%NProcShared={SP_PROCS}\n%Mem={SP_MEM}\n")
        file.write(f"%chk={name_without_extension}.chk\n")
        # made change here to reflect changing the order of single point jobs around
        if 'Pd' in unique_atoms or 'Zn' in unique_atoms:
            file.write(f"#N {SP_FUNCTIONAL}/gen pseudo=read {DISPERSION} {sp_1_properties}")
            file.write(f"\n\n{name_without_extension}\n\n")
            file.write(f"{charge_spin}")
            file.writelines(lines)
            file.write("\n")
            file.write(f"-{unique_atoms_list} 0\n{SP_BASIS_SET}\n****\n\n")
            if need_ecp:
                for atom in need_ecp:
                    file.write(f"-{atom} 0\n{SP_BASIS_SET} \n")
                file.write("\n")
        else:
            file.write(f"#N {SP_FUNCTIONAL}/{SP_BASIS_SET} {DISPERSION} {sp_1_properties}")
            file.write(f"\n\n{name_without_extension}\n\n")
            file.write(f"{charge_spin}")
            file.writelines(lines)
            file.write("\n")
            if need_ecp:
                for atom in need_ecp:
                    file.write(f"-{atom} 0\n{SP_BASIS_SET} \n")
                file.write("\n")

        # start of Link1 section, 2nd single point job
        file.write(f"--Link1--\n%NProcShared={SP_PROCS}\n%Mem={SP_MEM}\n")
        file.write(f"%chk={name_without_extension}.chk\n")
        # second part of the change made for the order of single point jobs
        if 'Pd' in unique_atoms or 'Zn' in unique_atoms:
            file.write(f"#N Guess=Read Geom=Check {SP_FUNCTIONAL}/gen pseudo=read {DISPERSION} {sp_2_metal_properties}")
            file.write(f"\n\n{name_without_extension}\n\n")
            file.write(f"{charge_spin}")
            file.write("\n")
            file.write(f"-{unique_atoms_list} 0\n{SP_BASIS_SET}\n****\n\n")
            if need_ecp:
                for atom in need_ecp:
                    file.write(f"-{atom} 0\n{SP_BASIS_SET}\n\n")
            file.write('\n\n\n')
        else:
            file.write(f"#N Guess=Read Geom=Check {SP_FUNCTIONAL}/{SP_BASIS_SET} {DISPERSION} {sp_2_nom_properties}")
            file.write(f"\n\n{name_without_extension}\n\n")
            file.write(f"{charge_spin}")
            file.write("\n")
            if need_ecp:
                for atom in need_ecp:
                    file.write(f"-{atom} 0\n{SP_BASIS_SET}\n\n")
            file.write('\n\n\n')

        print(f"Sucessfully converted {name_without_extension} to SPE G16 com file!")

# TODO: if we want to add generation of a .wfn file, need to add the <filename>.wfn line at the
#  bottom of the file and need the output=wfn command in the # line


def move_all_dft_inputs(results_dir: Path, dft_inputs_dir: Path):
    """
    This function moves all of the generated .com files into a separate directory.
    :param results_dir: Path to the results directory
    :param dft_inputs_dir: Path to the DFT inputs directory
    """
    results_dir = results_dir
    com_start_dir = results_dir / "split_conformers"
    com_target_dir = dft_inputs_dir
    os.makedirs(com_target_dir, exist_ok=True)

    files = os.listdir(com_start_dir)
    com_files = [file for file in files if file.endswith('.com')]

    for file in com_files:
        source_path = os.path.join(com_start_dir, file)
        destination_path = os.path.join(com_target_dir, file)
        shutil.move(source_path, destination_path)

    print('All .com files moved successfully')


def normal_termination(logfile: Path):
    """
    Function to check for normal termination in the Gaussian .log files
    Modified from Brittany Haas and Melissa Hardy's Gaussian post-processing script.
    :param logfile: Path to log file to check
    """
    log_filename = logfile.name
    log_stem = log_filename.split(".")[0]
    normalt = 0
    with open(logfile, "r") as inp:
        lines = inp.readlines()
        for line in lines:
            if "Normal termination" in line:
                normalt += 1

    # if the log files do not have 2 normal termination lines, they are moved to the term_error directory along
    # with their respective .com files

    term_error_dir = logfile.parent.parent / "term_error"

    if normalt != 2:
        print(f"{log_stem} did not terminate normally. \n")
        shutil.move(logfile, term_error_dir / f"{log_stem}.log")
        shutil.move(logfile.parent.parent / f"{log_stem}.com", term_error_dir / f"{log_stem}.com")
    else:
        print(f"{log_stem} terminated normally.")


def lowest_frequencies(logfile: Path):
    """
    Function to check for imaginary frequencies in the Gaussian .log files
    Modified from Brittany Haas and Melissa Hardy's Gaussian post-processing script.
    :param logfile: Path to log file to check
    """
    log_filename = logfile.name
    log_stem = log_filename.split(".")[0]

    imag_freq_dir = logfile.parent.parent / "imag"

    with open(logfile, "r") as inp:
        lines = inp.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if " and normal coordinates:" in line:
                lowfreq = lines[i + 3]
                if float(lowfreq[18:27]) <= 0:
                    print(f"{logfile}: {lowfreq[18:27]}")
                    shutil.move(logfile, imag_freq_dir / f"{log_stem}.log")
                    shutil.move(logfile.parent.parent / f"{log_stem}.com", imag_freq_dir / f"{log_stem}.com")


def g16_log_check(log_directory: Path):
    """
    This function parses the above two functions over all log files.
    :param log_directory: Path to the directory containing the log files
    """
    files = [file for file in log_directory.iterdir() if file.is_file()]
    log_files = [file for file in files if file.suffix == ".log"]

    for log_file in log_files:
        normal_termination(log_file)

    files = [file for file in log_directory.iterdir() if file.is_file()]
    log_files = [file for file in files if file.suffix == ".log"]
    for log_file in log_files:
        lowest_frequencies(log_file)


def process_opt_output(dft_output_dir: Path, dft_sp_dir: Path):
    """
    This function processes the output of the DFT optimization calculations.
    :param dft_output_dir: Path to the directory containing the DFT outputs
    :param dft_sp_dir: Path to the directory containing the DFT single point inputs
    """
    # directory for jobs which require re-submission
    os.makedirs(dft_output_dir / "resubmit", exist_ok=True)
    # directory for .log files with normal termination
    os.makedirs(dft_output_dir / "logs", exist_ok=True)
    # directory for .chk files where .logs have normal termination
    os.makedirs(dft_output_dir / "chks", exist_ok=True)
    # directory for .log and .com files which have a termination error
    os.makedirs(dft_output_dir / "term_error", exist_ok=True)
    # directory for .log and .com files which have at least 1 imaginary frequency
    os.makedirs(dft_output_dir / "imag", exist_ok=True)

    files = [file for file in dft_output_dir.iterdir() if file.is_file()]
    com_files = [file for file in files if file.suffix == ".com"]

    com_filename_only = []
    for i in com_files:
        file = i.name
        filename = file[:-4]
        com_filename_only.append(filename)

    for i in com_filename_only:
        com_file_path = dft_output_dir / f"{i}.com"
        log_file_path = dft_output_dir / f"{i}.log"

        # check if there is a log file for the com file and move to appropriate directory
        if log_file_path.exists():
            source_path = log_file_path
            destination_path = dft_output_dir / "logs" / f"{i}.log"
            shutil.move(source_path, destination_path)
        else:
            source_path = com_file_path
            destination_path = dft_output_dir / "resubmit" / f"{i}.com"
            shutil.move(source_path, destination_path)

    g16_log_check(dft_output_dir / "logs")

    log_file_dir = dft_output_dir / "logs"
    files = [file for file in log_file_dir.iterdir() if file.is_file()]
    log_files = [file for file in files if file.suffix == ".log"]

    for i in log_files:
        write_spc_com_template(i)
        log_filename = i.name
        chk_file = f"{log_filename[:-4]}.chk"
        chk_destination_dir = dft_output_dir / "chks"
        shutil.move(dft_output_dir / chk_file, chk_destination_dir / chk_file)

    files = [file for file in log_file_dir.iterdir() if file.is_file()]
    new_com_files = [file for file in files if file.suffix == ".com"]
    for i in new_com_files:
        write_g16_spc_com(i)

    # make directory in home for the single point jobs
    os.makedirs(dft_sp_dir, exist_ok=True)

    for i in new_com_files:
        shutil.move(i, dft_sp_dir / i.name)

    all_files = [file for file in dft_output_dir.iterdir() if file.is_file()]
    chks_to_remove = [file for file in all_files if file.suffix == ".chk"]

    for chk_to_remove in chks_to_remove:
        os.remove(chk_to_remove)


def process_spc_output(dft_output_dir: Path):
    """
    This function processes the output of the DFT single point calculations.
    :param dft_output_dir: Path to the directory containing the DFT outputs
    """
    # setup directories for .log file processing
    # directory for jobs which require re-submission
    os.makedirs(dft_output_dir / "resubmit", exist_ok=True)
    # directory for .log files with normal termination
    os.makedirs(dft_output_dir / "logs", exist_ok=True)
    # directory for .chk files where .logs have normal termination
    os.makedirs(dft_output_dir / "chks", exist_ok=True)
    # directory for .log and .com files which have a termination error
    os.makedirs(dft_output_dir / "term_error", exist_ok=True)
    # directory for .log and .com files which have at least 1 imaginary frequency
    os.makedirs(dft_output_dir / "imag", exist_ok=True)

    files = [file for file in dft_output_dir.iterdir() if file.is_file()]
    com_files = [file for file in files if file.suffix == ".com"]

    # create a list of the filenames without extension
    com_filename_only = []
    for i in com_files:
        file = i.name
        filename = file[:-4]
        com_filename_only.append(filename)

    for i in com_filename_only:
        com_file_path = dft_output_dir / f"{i}.com"
        log_file_path = dft_output_dir / f"{i}.log"

        # check if there is a log file for the com file and move to appropriate directory
        if log_file_path.exists():
            source_path = log_file_path
            destination_path = dft_output_dir / "logs" / f"{i}.log"
            shutil.move(source_path, destination_path)
        else:
            source_path = com_file_path
            destination_path = dft_output_dir / "resubmit" / f"{i}.com"
            shutil.move(source_path, destination_path)

    g16_log_check(dft_output_dir / "logs")

    log_file_dir = dft_output_dir / "logs"
    files = [file for file in log_file_dir.iterdir() if file.is_file()]
    log_files = [file for file in files if file.suffix == ".log"]

    for i in log_files:
        log_filename = i.name
        chk_file = f"{log_filename[:-4]}.chk"
        chk_destination_dir = dft_output_dir / "chks"
        shutil.move(dft_output_dir / chk_file, chk_destination_dir / chk_file)

    all_files = [file for file in dft_output_dir.iterdir() if file.is_file()]
    chks_to_remove = [file for file in all_files if file.suffix == ".chk"]

    for chk_to_remove in chks_to_remove:
        os.remove(chk_to_remove)

def format_spc_input(dft_sp_dir: Path):
    """
    This function formats Gaussian16 single-point input files from coms
    :param dft_sp_dir: Path to the directory containing the DFT single point inputs
    """

    files = [file for file in dft_sp_dir.iterdir() if file.is_file()]
    com_files = [file for file in files if file.suffix == ".com"]

    for i in com_files:
        write_g16_spc_com(i)
