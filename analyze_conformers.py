import argparse
import warnings
from pathlib import Path
from pp_conformers import conformers as conf
from pp_conformers import rmsd
from pp_conformers import dft
from pp_conformers import utils

warnings.filterwarnings("ignore")

DESCRIPTION = """


        ╔════════════════════════════════════════════════════════════════════╗
        ║   ____   ____    _      ___  ____   ____      _     ____  __   __  ║
        ║  |  _ \ |  _ \  | |    |_ _|| __ ) |  _ \    / \   |  _ \ \ \ / /  ║
        ║  | |_) || |_) | | |     | | |  _ \ | |_) |  / _ \  | |_) | \ V /   ║
        ║  |  __/ |  __/  | |___  | | | |_) ||  _ <  / ___ \ |  _ <   | |    ║
        ║  |_|    |_|     |_____||___||____/ |_| \_\/_/   \_\|_| \_\  |_|    ║
        ║                                                                    ║
        ╚════════════════════════════════════════════════════════════════════╝

                          -- BISPHOSPHINE LIGAND LIBRARY --
                               Conformer Analysis Tool


"""


def analyze_conformer(ligand_id: str, results_dir: Path, dft_inputs_dir: Path):
    """
    Analyze the conformers of a ligand from the CREST conformer search

    Performs the following steps:
    1. Analyze CREST output
    2. Run RMSD analysis
    3. Select conformers for DFT geometry refinement
    4. Convert conformers to Gaussian16 .com files
    5. Move .com files to DFT input directory

    :param ligand_id: Ligand ID
    :param results_dir: Path to the directory containing CREST results
    :param dft_inputs_dir: Path to the directory containing DFT input files
    """
    print(f"---{ligand_id} START---")
    print(f"\n*** Analyzing CREST output for ligand {ligand_id} ***")
    conf.crest_end(f"{ligand_id}.xyz", results_dir)

    print(f"\n*** Running RMSD analysis for ligand {ligand_id} ***")
    try:
        rmsd.run_conformer_rmsd(results_dir / f"{ligand_id}_crest_conformers.xyz")
        print("RMSD analysis successful!")
    except FileNotFoundError:
        print(f"ERROR! Could not find {ligand_id}_crest_conformers.xyz")

    print(f"\n*** Selecting {ligand_id} conformers for DFT geometry refinement ***")
    split_conformers_dir = results_dir / "split_conformers"
    files = [file for file in split_conformers_dir.iterdir() if file.is_file()]
    conformer_file_list = []
    for file in files:
        if ligand_id in file.name:
            conformer_file_list.append(file)
    print(f"Found {len(conformer_file_list)} total clustered conformers")

    xyz_files_for_dft = conf.select_files_for_dft(conformer_file_list)
    print(f"Selected {len(xyz_files_for_dft)} conformers for DFT geometry refinement")

    print(f"\n*** Converting {len(xyz_files_for_dft)} {ligand_id} conformers to Gaussian16 .com files ***")
    for file in xyz_files_for_dft:
        dft.write_opt_com_template(file)

    all_files = [file for file in split_conformers_dir.iterdir() if file.is_file()]
    for file in all_files:
        if {ligand_id} and ".com" in file.name:
            dft.write_g16_opt_com(file)

    dft.move_all_dft_inputs(results_dir, dft_inputs_dir)
    print(f"---{ligand_id} END---\n\n")


def get_run_conformers(crest_jobs_file: str) -> list:
    """
    Get the list of ligand IDs to analyze from the CREST jobs file
    :param crest_jobs_file: CREST jobs file
    :return: List of ligand IDs
    """
    with open(crest_jobs_file, "r") as file:
        try:
            crest_jobs = [line.strip() for line in file]
        except FileNotFoundError:
            print(f"ERROR! Cannot find {crest_jobs_file} file!")

    return crest_jobs


def analyze_all_conformers(crest_jobs_file: str, results_dir: Path, dft_inputs_dir: Path):
    """
    Analyze all conformers in the CREST jobs file
    :param crest_jobs_file: CREST jobs file
    :param results_dir: Path to the directory containing CREST results
    :param dft_inputs_dir: Path to the directory containing DFT input files
    """
    ligands = get_run_conformers(crest_jobs_file)
    for ligand_id in ligands:
        analyze_conformer(ligand_id=ligand_id,
                          results_dir=results_dir,
                          dft_inputs_dir=dft_inputs_dir)


def get_args() -> argparse.Namespace:
    """
    Get the arguments for the command line interface
    :return: Arguments
    """
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, 2, 40),
        allow_abbrev=False
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--ligand-id',
                       dest='ligand_id',
                       help='Ligand ID of conformers to analyze',
                       metavar='<ligand_id>'
                       )

    group.add_argument('-a', '--all',
                       dest='all',
                       action='store_true',
                       help='Analyze all conformers in conformers.txt file'
                       )

    parser.add_argument('--crest-results-dir',
                        dest='crest_results_dir',
                        help='Directory containing CREST results',
                        metavar='<crest_results_dir>',
                        default='crest_results'
                        )

    parser.add_argument('--crest-jobs-file',
                        dest='crest_jobs_file',
                        help='File containing CREST job names',
                        metavar='<crest_jobs_file>',
                        default='conformers.txt'
                        )

    parser.add_argument('--dft-opt-dir',
                        dest='dft_inputs_dir',
                        help='Directory containing DFT geometry optimization input files',
                        metavar='<dft_opt_dir>',
                        default='dft_inputs'
                        )

    args = parser.parse_args()

    return args


def main():
    """
    Main function for analyzing conformers
    """
    args = get_args()

    # crest_xtb_loaded = utils.is_crest_loaded()
    # openbabel_loaded = utils.is_openbabel_loaded()
    # if not crest_xtb_loaded:    
    #     exit(1)
    # if not openbabel_loaded:
        # exit(1)

    ligand_id = args.ligand_id
    run_all = args.all
    crest_results_dir = Path(f"./{args.crest_results_dir}")
    crest_jobs_file = args.crest_jobs_file
    dft_inputs_dir = Path(f"./{args.dft_inputs_dir}")

    if run_all:
        analyze_all_conformers(crest_jobs_file=crest_jobs_file,
                               results_dir=crest_results_dir,
                               dft_inputs_dir=dft_inputs_dir)

    elif ligand_id:
        analyze_conformer(ligand_id=ligand_id,
                          results_dir=crest_results_dir,
                          dft_inputs_dir=dft_inputs_dir)


if __name__ == '__main__':
    main()
