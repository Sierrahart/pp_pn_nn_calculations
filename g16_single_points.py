import argparse
from pathlib import Path
from conf_search_tools import dft

# TODO make a function here which analyzes the sp outputs - will need to change the name of the script

DESCRIPTION = """


    ╔══════════════════════════════════════════════════════════════════════════════╗
    ║   ____    ___    ____      _      ___  ____   ____      _     ____  __   __  ║
    ║  | __ )  |_ _|  |  _ \    | |    |_ _|| __ ) |  _ \    / \   |  _ \ \ \ / /  ║
    ║  |  _ \   | |   | | \ |   | |     | | |  _ \ | |_) |  / _ \  | |_) | \ V /   ║
    ║  | |_) |  | |   | |_/ |   | |___  | | | |_) ||  _ <  / ___ \ |  _ <   | |    ║
    ║  |____/  |___|  |____/    |_____||___||____/ |_| \_\/_/   \_\|_| \_\  |_|    ║
    ║                                                                              ║
    ╚══════════════════════════════════════════════════════════════════════════════╝

                          -- BIDENTATE LIGAND LIBRARY --
                Gaussian Single Point Job Generation and Analysis Tool

"""




def generate_spe_coms(dft_opt_dir: Path, dft_sp_dir: Path):
    """
    Generate Gaussian16 single-point input files from optimization .log files.
    :param dft_opt_dir: Path to directory containing DFT optimization input files.
    :param dft_sp_dir: Path to directory containing DFT single-point input files.
    """
    dft.process_opt_output(dft_opt_dir, dft_sp_dir)


def process_spe_coms(dft_sp_dir: Path):
    """
    Process Gaussian16 single-point output files.
    :param dft_sp_dir: Path to directory containing DFT single-point input files.
    """
    dft.process_spc_output(dft_sp_dir)

def format_spe_coms(dft_sp_dir: Path):
    """
    Formate Gaussian16 single-point input files.
    :param dft_sp_dir: Path to directory containing DFT single-point input files.
    """
    dft.format_spc_input(dft_sp_dir)
    
def get_args() -> argparse.Namespace:
    """
    Get command-line arguments.
    :return: Arguments
    """
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, 2, 40),
        allow_abbrev=False
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--generate',
                       action='store_true',
                       help='Generate all Gaussian16 single-point jobs from optimization .log files.'
                       )

    group.add_argument('--process',
                       action='store_true',
                       help='Process Gaussian16 single-point jobs.'
                       )
    
    group.add_argument('--format',
                       action='store_true',
                       help='Format Gaussian16 single-point jobs.'
                       )

    parser.add_argument('--dft-opt-dir',
                        dest='dft_inputs_dir',
                        help='Directory containing DFT geometry optimization input files.',
                        metavar='<dft_opt_dir>',
                        default='dft_inputs'
                        )

    parser.add_argument('--dft-sp-dir',
                        dest='dft_sp_dir',
                        help='Directory containing DFT single-point input files.',
                        metavar='<dft_sp_dir>',
                        default='dft_sp_inputs'
                        )

    args = parser.parse_args()

    return args


def main():
    """
    Main function.
    """
    args = get_args()

    generate = args.generate
    process = args.process
    format = args.format
    dft_opt_dir = Path(f"./{args.dft_inputs_dir}")
    dft_sp_dir = Path(f"./{args.dft_sp_dir}")

    if generate:
        generate_spe_coms(dft_opt_dir, dft_sp_dir)
    elif process:
        process_spe_coms(dft_sp_dir)
    elif format:
        format_spe_coms(dft_sp_dir)


if __name__ == '__main__':
    main()
