import os
import argparse
from conf_search_tools import conformers as conf
from conf_search_tools import utils

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
                              Conformer Analysis Tool

"""



# def verify_xyz_filename(xyz_file: str) -> bool:
#     """
#     Verify that the filename is valid and follows the format: <ligand_id>-<metal_fragment>.xyz
#     :param xyz_file: str
#     :return: bool
#     """
#     verify_return = utils.verify_xyz_filename(xyz_file)

#     return verify_return


# def xyz_file_rename(xyz_file: str) -> str:
#     """
#     Rename the .xyz file to the correct format: <ligand_id>-<metal_fragment>.xyz
#     :param xyz_file: str
#     :return: str
#     """
#     rename_return = utils.xyz_file_rename(xyz_file)

#     return rename_return


def run_one_file(xyz_file: str):
    """
    Run a conformational search on a single .xyz file
    First write the constraints file, then execute the CREST job
    :param xyz_file: str
    """
    print(f"Running conformational search on {xyz_file}")
    conf.write_constraints_file(xyz_file)
    conf.execute_crest()


def run_all_files(cur_dir: str):
    """
    Run a conformational search on all .xyz files in the specified directory
    :param cur_dir: str
    """
    print(f"Running conformational searches on all files in {cur_dir}")
    all_files = os.listdir(cur_dir)

    for file in all_files:
        if ".xyz" in file:
            # not_valid_name = verify_xyz_filename(file)
            # if not_valid_name:
                # print(f"\033[91m[ERROR!] The filename {file} is not valid.\033[0m")
            # else:
            conf.write_constraints_file(file)
    conf.execute_crest()


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

    group.add_argument('--xyz-file',
                       dest='xyz_file',
                       help='Initial XYZ coordinate file to run',
                       metavar='<xyz_filename>'
                       )

    group.add_argument('-a', '--all',
                       dest='all',
                       action='store_true',
                       help='Run all XYZ coordinate files in directory'
                       )

    parser.add_argument('--xyz-verify',
                        dest='xyz_verify',
                        help='Verify XYZ coordinate file name',
                        metavar='<xyz_filename>'
                        )

    parser.add_argument('--xyz-dir',
                        dest='xyz_dir',
                        help='Directory containing XYZ coordinate files',
                        metavar='<xyz_directory>',
                        default=os.getcwd()
                        )

    args = parser.parse_args()

    return args


def main():
    """
    Main function to run the conformational search
    """
    args = get_args()

    xyz_file = args.xyz_file
    # xyz_file_to_verify = args.xyz_verify
    run_all = args.all
    xyz_dir = args.xyz_dir

    if run_all:
        run_all_files(xyz_dir)

    elif xyz_file:
        if os.path.isfile(xyz_file):
            # not_valid_name = verify_xyz_filename(xyz_file)
            # if not_valid_name:
                # print(f"\033[91m[ERROR!] The filename {xyz_file} is not valid.\033[0m")
            # else:
            run_one_file(xyz_file)
        else:
            print(f"\033[91m[ERROR!] {xyz_file} does not exist!\033[0m")

    # elif xyz_file_to_verify:
        # if os.path.isfile(xyz_file_to_verify):
            # not_valid_name = verify_xyz_filename(xyz_file_to_verify)
# 
            # if not_valid_name:
                # print(f"\033[91m[ERROR!] The filename {xyz_file_to_verify} is not valid.\033[0m")
                # input("\nPress enter to rename the file...")
                # new_name = xyz_file_rename(xyz_file_to_verify)
                # os.rename(xyz_file_to_verify, new_name)
                # print(f"\033[92mFile successfully renamed to {new_name}!\033[0m")
                # run_search = input(f"Would you like to run the conformational search on {new_name}? [y/n]: ")
                # if run_search.lower() == 'y':
                    # run_one_file(new_name)
                # elif run_search.lower() == 'n':
                    # print("Exiting...")
        #             exit()
        #     else:
        #         print("\033[92m[INFO] The filename is valid.\033[0m")
        #         run_search = input(f"Would you like to run the conformational search on {xyz_file_to_verify}? [y/n]: ")
        #         if run_search.lower() == 'y':
        #             run_one_file(xyz_file_to_verify)
        #         elif run_search.lower() == 'n':
        #             print("Exiting...")
        #             exit()
        # else:
        #     print(f"\033[91m[ERROR!] {xyz_file_to_verify} does not exist!\033[0m")


if __name__ == '__main__':
    main()
