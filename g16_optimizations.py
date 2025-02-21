import argparse
from pathlib import Path
from pp_conformers import dft

# TODO make a function here which analyzes the sp outputs - will need to change the name of the script

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
                 Gaussian16 Single-Point Correction Calculation Tool


"""


def generate_opt_coms(dft_opt_dir: Path):
    """
    Generate Gaussian16 single-point input files from optimization .log files.
    :param dft_opt_dir: Path to directory containing DFT optimization input files.
    :param dft_sp_dir: Path to directory containing DFT single-point input files.
    """
    for file in dft_opt_dir.iterdir():
        if file.suffix == '.com':
            dft.write_g16_opt_com(file)


# def process_opt_coms(dft_sp_dir: Path):
#     """
#     Process Gaussian16 single-point output files.
#     :param dft_sp_dir: Path to directory containing DFT single-point input files.
#     """
#     dft.process_spc_output(dft_sp_dir)

    
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
                       help='Generate all Gaussian16 optimization jobs from unformatted .com files.'
                       )

    parser.add_argument('--directory',
                        dest='dft_inputs_dir',
                        help='Directory containing unformatted DFT geometry optimization input files.',
                        metavar='<dft_opt_dir>',
                        default='dft_inputs'
                        )

    args = parser.parse_args()

    return args


def main():
    """
    Main function.
    """
    args = get_args()
    print(args)
    generate = args.generate
    dft_opt_dir = Path(f"./{args.dft_inputs_dir}")

    if generate:
        generate_opt_coms(dft_opt_dir)


if __name__ == '__main__':
    main()
