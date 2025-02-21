"""
GENERATION 2 BISPHOSPHINE LIBRARY WORKFLOW
-Conformer RMSD Analysis Functions-
"""

import MDAnalysis as mda
import pandas as pd
import yaml
from pathlib import Path
from MDAnalysis.analysis import rms


def run_conformer_rmsd(xyz_file: Path):
    """
    Function for running RMSD analysis on the CREST-generated conformational ensembles
    Incorporates descriptors into a YAML file
    Only for Pd complexes
    :param xyz_file: Path to the xyz file
    """

    # extract ligand name from CREST file
    xyz_filename = xyz_file.name
    ligand_name_extract = xyz_filename.split("_")
    ligand_name = ligand_name_extract[0]

    # open xyz coordinate file and convert to MD analysis Universe
    u = mda.Universe(xyz_file)

    # gather file information
    ligand_id = ligand_name_extract
    filename = xyz_file
    no_atoms = len(u.atoms)
    no_frames = len(u.trajectory)
    file_info_list = [str(ligand_name), filename, no_atoms, no_frames]

    # run rmsd trajectory analysis
    u.trajectory[0]
    rmsd_analysis = rms.RMSD(u,
                             select="all and not (element Pd or element Cl or element H)",
                             groupselections=['all', 'element P', 'element C', 'element O'],
                             )

    rmsd_analysis.run()

    # generate dataframe
    df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                      columns=['all', 'all2', 'P', 'C', 'O'],
                      index=rmsd_analysis.results.rmsd[:, 1])
    df.index.name = 'Frame'

    # calculate mean, min, max, sd for each of the group selections
    all2_mean = df['all2'].mean()
    p_mean = df['P'].mean()
    c_mean = df['C'].mean()
    o_mean = df['O'].mean()

    all2_min = df['all2'].min()
    p_min = df['P'].min()
    c_min = df['C'].min()
    o_min = df['O'].min()

    all2_max = df['all2'].max()
    p_max = df['P'].max()
    c_max = df['C'].max()
    o_max = df['O'].max()

    all2_sd = df['all2'].std()
    p_sd = df['P'].std()
    c_sd = df['C'].std()
    o_sd = df['O'].std()

    # compile descriptors
    all2 = [float(all2_mean), float(all2_min), float(all2_max), float(all2_sd)]
    p = [float(p_mean), float(p_min), float(p_max), float(p_sd)]
    c = [float(c_mean), float(c_min), float(c_max), float(c_sd)]
    o = [float(o_mean), float(o_min), float(o_max), float(o_sd)]

    # make descriptor dictionary
    desc_dict = {
        'file_info': file_info_list,
        'all': all2,
        'P': p,
        'C': c,
        'O': o,
        'descriptor_headers': ["mean", "min", "max", "std_dev"],
    }

    # generate YAML file containing descriptors
    with open(xyz_file.parent / f"{ligand_name}_rmsd_analysis.yml", "w") as file:
        desc_file = yaml.dump(desc_dict, file)
