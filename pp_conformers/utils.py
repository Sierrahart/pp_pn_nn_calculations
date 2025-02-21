"""
GENERATION 2 BISPHOSPHINE LIBRARY WORKFLOW
-General Utility Functions Script-
"""

from morfeus.io import read_xyz
from morfeus import BiteAngle
from morfeus.utils import get_connectivity_matrix
import numpy as np
from pathlib import Path

# define lists of accepted donor atoms and metal atoms Current reference bond lengths only exist for PP, PN and NN
# ligands in conformers.py script.
DONOR_ATOMS = ['P', 'N', 'O']
METAL_ATOMS = ['Pd', 'Zn', 'Co', 'Ni']
METAL_LIGAND = ['Cl']


def verify_xyz_filename(filename: str) -> bool:
    """
    This function verifies the filename of an XYZ file to ensure it is in the correct format.
    The correct format is:
    ppXXXXXX-METALLIGAND.xyz
    where XXXXXX is a 6-digit ligand number, METAL is a valid metal and LIGAND is a valid metal ligand
    combination, e.g. PdCl2.

    Filenames are verified by checking the following:
    - The filename ends with .xyz
    - The filename starts with pp
    - The filename contains a hyphen to separate the metal fragment from the ID
    - The filename contains a valid metal fragment
    - The filename contains a valid metal ligand fragment
    - The filename contains a 6 digit ligand number

    :param filename: Name of the file to be verified
    :return: Boolean
    """

    if not filename.endswith(".xyz"):
        raise Exception("\033[91mFile is not an XYZ file. Please provide an XYZ file.\033[0m")

    error = False

    filename_no_ext = filename.split(".")[0]

    if not filename_no_ext.startswith("pp"):
        print(f"\033[91m[ERROR!] File {filename} does not start with 'pp'.\033[0m")
        error = True

    if "-" not in filename_no_ext:
        print(f"\033[91m[ERROR!] File {filename} does not contain a hyphen to separate the metal "
              f"fragment from the ID.\033[0m")
        error = True

    valid_metal = False
    for metal in METAL_ATOMS:
        if metal in filename_no_ext:
            valid_metal = True
            break
    if not valid_metal:
        print(f"\033[91m[ERROR!] File {filename} does not contain a valid metal fragment ({METAL_ATOMS}).\033[0m")
        error = True

    valid_metal_ligand = False
    for ligand in METAL_LIGAND:
        if ligand + '2' in filename_no_ext:
            valid_metal_ligand = True
            break
    if not valid_metal_ligand:
        print(f"\033[91m[ERROR!] File {filename} does not contain a valid metal ligand fragment "
              f"({[ligand + '2' for ligand in METAL_LIGAND]}).\033[0m")
        error = True

    ligand_number = filename_no_ext[2:].split("-")[0]
    if not ligand_number.isdigit():
        print(f"\033[91m[ERROR!] {ligand_number} is not a valid ligand number.\033[0m")
        error = True
    if len(ligand_number) != 6:
        print(f"\033[91m[ERROR!] {ligand_number} is not a valid 6 digit ligand number.\033[0m")
        error = True

    return error


def xyz_file_rename(filename: str) -> str:
    """
    This function renames an XYZ file to a standard format.
    :param filename: Name of the file to be renamed
    :return: New filename
    """
    print(f"Renaming {filename} to a standard format...")
    while True:
        new_ligand_id = input(f"Enter a new ligand ID for {filename}: ")
        new_ligand_id_zfilled = new_ligand_id.zfill(6)
        if new_ligand_id_zfilled.isdigit() and len(new_ligand_id_zfilled) == 6:
            break
        else:
            print(f"\033[91m[ERROR!] Ligand ID {new_ligand_id} is not a valid ID number. Please enter a "
                  f"valid ID number.\033[0m")

    use_defaults = input("Would you like to use the default metal center and ligand? (PdCl2) [y/n]: ")
    if use_defaults.lower() == "y":
        new_metal_center = "Pd"
        new_metal_ligand = "Cl"
    else:
        while True:
            new_metal_center = input(f"Enter a new valid metal center for {filename}: ")
            if new_metal_center in METAL_ATOMS:
                break
            else:
                print(f"\033[91m[ERROR!] Ligand fragment {new_metal_center} is not a valid metal center. "
                      f"Please enter a valid metal center.\033[0m")

        while True:
            new_metal_ligand = input(f"Enter a new valid metal ligand for {filename}: ")
            if new_metal_ligand in METAL_LIGAND:
                break
            else:
                print(f"\033[91m[ERROR!] Ligand fragment {new_metal_ligand} is not a valid metal ligand. "
                      f"Please enter a valid metal ligand.\033[0m")

    new_filename = f"pp{new_ligand_id_zfilled}-{new_metal_center}{new_metal_ligand}2.xyz"
    print(f"\033[92mOld filename {filename} will be renamed to {new_filename}\033[0m")
    input("\nPress enter to confirm...")

    return new_filename


def find_metal_and_donor_atoms(xyz_file: Path) -> tuple:
    """
    This function uses MORFEUS to determines the atom labels and corresponding coordinates for the metal atom
    and the two donor atoms.
    Metals currently supported: Pd, Zn, Co, Ni
    Donor atoms currently supported: P
    Looks also for Cl atoms bound to the complex also
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Tuple containing the metal atom, donor atoms, donor elements, metal ligand atoms and metal ligand elements
    """
    elements, coordinates = read_xyz(xyz_file)

    m_index = None
    for metal in METAL_ATOMS:
        if metal in elements:
            m_index = np.where(elements == metal)[0][0]
            break  # loop breaks once a metal atom is identified.

    metal_atom = m_index + 1  # add 1 to the index to remove the 0-indexing

    try:
        if xyz_file.startswith('bpmo'):
            xyz_file = xyz_file
    except:
        xyz_file = str(xyz_file).split('/')[-1]

    if xyz_file.startswith('bpmo'):
        connectivity_matrix = get_connectivity_matrix(coordinates, elements, radii_type='pyykko', scale_factor=1.4)
        connected_atoms = np.where(connectivity_matrix[m_index, :])[0]

        phosphine_atoms = []
        oxide_atom = []
        phosphine_oxide_distances = []
        donor_atoms = []
        donor_atom_type = []
        metal_ligand_atoms = []
        metal_ligand_type = []

        for atom in connected_atoms:
            if elements[atom] == 'P':
                phosphine_atoms.append(atom)
            if elements[atom] == 'O':
                oxide_atom.append(atom)
            if elements[atom] in METAL_LIGAND:
                metal_ligand_atoms.append(atom)
                metal_ligand_type.append(elements[atom])
        
        if len(phosphine_atoms) > 1:
            for phosphine_atom in phosphine_atoms:
                x_distance = coordinates[phosphine_atom][0] - coordinates[oxide_atom[0]][0]
                y_distance = coordinates[phosphine_atom][1] - coordinates[oxide_atom[0]][1]
                z_distance = coordinates[phosphine_atom][2] - coordinates[oxide_atom[0]][2]
                phosphine_oxide_distance = np.sum([x_distance**2, y_distance**2, z_distance**2])
                phosphine_oxide_distances.append(np.sqrt(phosphine_oxide_distance))

            if phosphine_oxide_distances[0] > phosphine_oxide_distances[1]:
                donor_atoms.append(phosphine_atoms[0])
                donor_atoms.append(oxide_atom[0])
                donor_atom_type.append(elements[phosphine_atoms[0]])
                donor_atom_type.append(elements[oxide_atom[0]])
            else:
                donor_atoms.append(phosphine_atoms[1])
                donor_atoms.append(oxide_atom[0])
                donor_atom_type.append(elements[phosphine_atoms[1]])
                donor_atom_type.append(elements[oxide_atom[0]])
        
        elif len(phosphine_atoms) == 1:
            donor_atoms.append(phosphine_atoms[0])
            donor_atoms.append(oxide_atom[0])       
            donor_atom_type.append(elements[phosphine_atoms[0]])
            donor_atom_type.append(elements[oxide_atom[0]])     

    else:
        connectivity_matrix = get_connectivity_matrix(coordinates, elements, radii_type='pyykko', scale_factor=1.4)
        connected_atoms = np.where(connectivity_matrix[m_index, :])[0]

        donor_atoms = []
        donor_atom_type = []
        metal_ligand_atoms = []
        metal_ligand_type = []

        for atom in connected_atoms:
            if elements[atom] in DONOR_ATOMS:
                donor_atoms.append(atom)
                donor_atom_type.append(elements[atom])
            if elements[atom] in METAL_LIGAND:
                metal_ligand_atoms.append(atom)
                metal_ligand_type.append(elements[atom])

    if len(donor_atoms) != 2:
        print(f'{xyz_file} has {len(donor_atoms)} donor atoms.')

    # add 1 to the indices to remove the 0-indexing
    donor_atom1 = donor_atoms[0] + 1
    donor_atom2 = donor_atoms[1] + 1
    donor_elements = (donor_atom_type[0], donor_atom_type[1])
    metal_ligand_atom1 = metal_ligand_atoms[0] + 1
    metal_ligand_atom2 = metal_ligand_atoms[1] + 1
    metal_ligand_elements = (metal_ligand_type[0], metal_ligand_type[1])

    return (metal_atom, donor_atom1, donor_atom2, donor_elements, metal_ligand_atom1, metal_ligand_atom2,
            metal_ligand_elements)


def find_ferrocene_atoms(xyz_file: Path) -> tuple:
    """
    This function uses MORFEUS to determine the atom labels and corresponding coordinates for the Fe atom and the
    carbon atoms in a ferrocene ligand.
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Tuple containing the Fe atom and carbon atoms
    """
    elements, coordinates = read_xyz(xyz_file)

    fe_index = np.where(elements == 'Fe')[0][0]
    fe_atom = fe_index + 1

    connectivity_matrix = get_connectivity_matrix(coordinates, elements, radii_type='pyykko', scale_factor=1.2)
    connected_atoms = np.where(connectivity_matrix[fe_index, :])[0]

    cp_atoms = []
    for atom in connected_atoms:
        if elements[atom] == 'C':
            cp_atoms.append(atom)
    if len(cp_atoms) < 10:
        raise Exception('Ferrocene fragment found is invalid (has less than 10 carbons), check.')

    cp_carbons = []
    for atom in cp_atoms:
        carbon = atom + 1
        cp_carbons.append(carbon)

    return fe_atom, cp_carbons


def geom_from_xyz(xyz_file: Path) -> list:
    """
    This function extracts coordinates from an XYZ file and returns as a list.
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: List of coordinates
    """
    coordinates = []

    with open(xyz_file, "r") as f:
        next(f)
        next(f)
        for line in f:
            coordinates.append(line)

    return coordinates


def get_molecular_formula(xyz_file: Path) -> str:
    """
    This function determines the molecular formula of a complex based on the coordinates in an xyz file.
    Used in comparing before and after CREST conformational searches to check for atom changes.
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Molecular formula of the complex
    """
    coordinates = geom_from_xyz(xyz_file)

    atom_list = []
    for coordinate in coordinates:
        atom_list.append(coordinate.split("  ")[0])

    atoms = []
    for atom in atom_list:
        nospace = "".join(filter(lambda x: x != " ", atom))
        atoms.append(nospace)

    atom_counts = {}
    for atom in atoms:
        if atom in atom_counts:
            atom_counts[atom] += 1
        else:
            atom_counts[atom] = 1

    sorted_elements = sorted(atom_counts.items(), key=lambda x: x[0])

    molecular_formula = ""
    for element, count in sorted_elements:
        if count > 1:
            molecular_formula += element + str(count)
        else:
            molecular_formula += element

    return molecular_formula


def get_xtb_energy(xyz_file: str) -> float:
    """
    This function gets the xTB energy from an XYZ coordinate file (generated from CREST).
    Extracts from the XYZ title line.
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Energy of the complex
    """
    with open(xyz_file, "r") as f:
        energy = float(f.readlines()[1])

    return energy


def get_bite_angle(xyz_file: Path) -> float:
    """
    This function calculates the bite angle of a complex using MORFEUS.
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Bite angle of the complex
    """
    elements, coordinates = read_xyz(xyz_file)
    (metal_atom, donor_atom1, donor_atom2, donor_elements, metal_ligand_atom1, metal_ligand_atom2,
     metal_ligand_element) = find_metal_and_donor_atoms(xyz_file)

    ba = BiteAngle(coordinates, metal_atom, donor_atom1, donor_atom2)
    bite_angle = ba.angle

    return bite_angle


def get_donor_atom_distance(xyz_file: Path) -> tuple:
    """
    This function calculates the distance between the two donor atoms in a complex.
    This serves a check to make sure the no major structural changes occured before and after CREST conformational
    search.
    Ligand donor atoms and metal atoms must be in the lists at the top of the script. Accepted bond lengths between
    donor atoms have only been incorporated for P,P P,N and N,N donor atoms (conformers.py script)
    :param xyz_file: Name of the XYZ file to be analyzed
    :return: Tuple containing the distance between the donor atoms and the donor elements
    """

    elements, coordinates = read_xyz(xyz_file)
    (metal_atom, donor_atom1, donor_atom2, donor_elements, metal_ligand_atom1, metal_ligand_atom2,
     metal_ligand_elements) = find_metal_and_donor_atoms(xyz_file)

    # calculate distance between donor atoms
    a = coordinates[donor_atom1]
    b = coordinates[donor_atom2]
    ba = np.array(a) - np.array(b)
    distance = round(np.linalg.norm(ba), 5)

    return distance, donor_elements


def geom_from_gaussian(log_file: Path) -> list:
    """
    This function extracts coordinates from a Gaussian log file and returns as a list.
    :param log_file: Path to the Gaussian log file to be analyzed
    :return: List of coordinates
    """
    with open(log_file, 'r') as file:
        loglines = file.readlines()

    starts = [i for i, line in enumerate(loglines) if "1\\1\\" in line]
    ends = [i for i, line in enumerate(loglines) if "@" in line]

    streams = [
        "".join(loglines[i][1:-1] for i in range(start, end + 1)).split("\\")
        for start, end in zip(starts, ends)
    ]

    coordinates = []

    for item in streams[-1][16:]:
        if item == "":
            break
        coordinates.append([item.split(",")[0], float(item.split(",")[-3]), float(item.split(",")[-2]),
                            float(item.split(",")[-1])])

    return coordinates
