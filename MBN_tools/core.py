"""
MBN_tools package

Core tools for working with MBN Explorer IO files and simulations.
"""

import numpy as np
import warnings
from typing import Union, Optional, List, Dict, Tuple


def create_structured_array(atoms: List[str], frames: Union[np.ndarray, List[np.ndarray]]) -> np.ndarray:
    """
    Create a single 3D structured numpy array from a list of atom types and a nested list/array of coordinates for multiple frames.
    
    Parameters:
        atoms (list of str): A list of atom type labels.
        frames (numpy array or list): A list or array of 3D coordinates for multiple frames 
                                      (shape: MxNx3, where M is the number of frames, and N is the number of atoms).
        
    Returns:
        numpy.ndarray: A 3D structured array with 'atoms' and 'coordinates' fields, one for each frame.
    """
    
    frames = np.array(frames)
    
    dtype = [('atoms', 'U10'), ('coordinates', 'f8', 3)]
    
    structured_array = np.array(
        [[(atom, coord) for atom, coord in zip(atoms, frame)] for frame in frames],
        dtype=dtype
    )
    
    return structured_array


def run_task(task_file: str, MBN_path: str, show_output: bool = False) -> Tuple[str, str]:
    """
    Run an MBN Explorer simulation using a specified Task file and return the standard output and error.

    Parameters:
        task_file (str): Path to the Task file.
        MBN_path (str): Path to the MBN Explorer executable.
        show_output (bool, optional): If True, prints the simulation output to the screen. Default is False.

    Returns:
        tuple: (stdout, stderr) where stdout and stderr are strings containing the standard output and error of the simulation.
    """
        
    import subprocess
    result = subprocess.run(MBN_path + ' -t ' + task_file, capture_output=True, text=True, creationflags=0x08000000)

    if show_output:
        print('Output:')
        if len(result.stdout) > 0:
            print(result.stdout)
        else:
            print(result.stderr)

    return result.stdout, result.stderr


def read_task(task_file: str, flatten: bool = False) -> Dict[str, Dict[str, str]]:
    """
    Read a Task file and split it into a dictionary of options and their parameters.

    Parameters:
        task_file (str): Path to the Task file.
        flatten (bool, optional): If True, flattens the nested dictionary into a single dictionary. Default is False.

    Returns:
        dict: A dictionary of Task file options. If flatten=True, a flattened dictionary ignoring Task file sections.
    """
    
    file_options = {}
    with open(task_file) as f:
        data = f.read().split('\n')
    
    section_options = {}
    for line in data:
        try:
            if line[0] ==';' and len(line)>1:
                if len(section_options)>0:
                    file_options[section] = section_options
                    section_options = {}
                section = line.replace(';','').strip()
            if '=' in line:
                line=line.split(' = ')
                section_options[line[0].strip()] = str(line[1].strip())
        except IndexError:
            pass
    file_options[section] = section_options
    
    if flatten:
        file_options_flat = {}
        for key in file_options.keys():
            file_options_flat.update(file_options[key])
        return file_options_flat
    else:
        return file_options
    

def write_task(task_file: str, file_options: Dict[str, Dict[str, str]]) -> None:
    """
    Write a dictionary of Task file options to a Task file.

    Parameters:
        task_file (str): Path to the Task file.
        file_options (dict): A dictionary of Task file options, where the keys are section headers and the values are dictionaries of parameters.
    
    Returns:
        None: Writes the Task file.
    """
    
    with open(task_file, 'w') as f:
        f.write(';\n')
        for key in file_options:
            f.write('\n; '+key+'\n')
            for sub_key in file_options[key]:
                if sub_key == 'Random':
                    f.write('\n'+'{:<31}'.format(sub_key)+'= '+file_options[key][sub_key]+'\n')
                else:
                    f.write('{:<31}'.format(sub_key)+'= '+file_options[key][sub_key]+'\n')


def read_xyz(xyz_file: str) -> np.ndarray:
    """
    Read an XYZ file and return its contents as a structured array.

    Parameters:
        xyz_file (str): Path to the XYZ file.

    Returns:
        numpy.ndarray: A structured array with fields 'atoms' and 'coordinates'.
    """
    
    with open(xyz_file, 'r') as f:
        data = [i.split() for i in f.read().split('\n')[2:-1]]
        atoms = [i[0] for i in data]
        coordinates = [[[float(i[1]), float(i[2]), float(i[3])] for i in data]]

    xyz_data = create_structured_array(atoms, coordinates)

    return xyz_data


def write_xyz(coords: np.ndarray, xyz_file: str, frame: int = 0) -> None:
    """
    Write coordinates to an XYZ file from a structured array.

    Parameters:
        coords (numpy.ndarray): A structured array of coordinates.
        xyz_file (str): Path to the XYZ file.
        frame (int, optional): The simulation frame to use if multiple are specified. Defaults to the first frame.
    
    Returns:
        None: Writes the XYZ file.
    """
    
    with open(xyz_file, 'w') as f:
        f.write(str(len(coords[frame]))+'\n')
        f.write('Type name\t\t\tPosition X\t\t\tPosition Y\t\t\tPosition Z\n')
        for line in coords[frame]:
            f.write(line['atoms']+'\t\t\t\t'+'\t\t'.join(['{:.8e}'.format(float(x)) for x in line['coordinates']])+'\n')


def read_trajectory(dcd_file: str, xyz_file: str, frame: Optional[int] = None) -> np.ndarray:
    """
    Read coordinates from a DCD trajectory file and match them with atom labels from an XYZ file.

    Parameters:
        dcd_file (str): Path to the DCD file.
        xyz_file (str): Path to the XYZ file containing atom labels.
        frame (int, optional): The specific frame to return. If None, returns all frames.

    Returns:
        numpy.ndarray: A structured array of coordinates for the specified frame with fields 'atoms' and 'coordinates'.
    """
    
    import mdtraj.formats as md
    
    with md.DCDTrajectoryFile(dcd_file) as f:
        xyz, cell_lengths, cell_angles = f.read()
        #coords = (xyz)

    atoms = read_xyz(xyz_file)[0]['atoms']

    coords = create_structured_array(atoms, (xyz))

    if frame:
        return coords[frame]
    else:
        return coords
        
        
def xyz_to_pdb(xyz_file: str, pdb_file: str) -> None:
    """
    Update coordinates in a PDB file using new coordinates from an XYZ file.

    Parameters:
        xyz_file (str): Path to the XYZ file containing new coordinates.
        pdb_file (str): Path to the PDB file with old coordinates.

    Returns:
        None: Writes a new PDB file with updated coordinates prefixed by 'new_'.
    """
    
    with open(xyz_file) as xyz:
        lines = xyz.readlines()
        xyz_list = [line.strip().split() for line in lines[2:]]

    with open(pdb_file) as pdb:
        lines = pdb.readlines()
        head = []
        atoms = []
        foot = []
        atoms_done = False
        for line in lines:
            if line.startswith('ATOM'):
                atoms.append(line.strip().split())
            elif line.startswith('END'):
                atoms_done = True
                foot.append(line)
            elif atoms_done:
                foot.append(line)
            else:
                head.append(line)

    with open('new_' + pdb_file, 'w') as new_pdb:
        for line in head:
            new_pdb.write(line)
        for idx, atom in enumerate(atoms):
            updated_atom = "{ATOM}{atom_num} {atom_name}{alt_loc_ind}{res_name} {chain_id}{res_seq_num}{res_code}   {x_coord}{y_coord}{z_coord}{occ}{temp}      {seg_id}{element}{charge}\n".format( #The values on the lines below can be edited to comply with other PDB types
                ATOM=atom[0].ljust(6),
                atom_num=atom[1].rjust(5),
                atom_name=atom[2].ljust(4),
                alt_loc_ind=' ',
                res_name=atom[3].rjust(3),
                chain_id=' ',
                res_seq_num=atom[4].rjust(4),
                res_code=' ',
                x_coord=str('%3.3f' % (float(xyz_list[idx][1]))).rjust(8),
                y_coord=str('%3.3f' % (float(xyz_list[idx][2]))).rjust(8),
                z_coord=str('%3.3f' % (float(xyz_list[idx][3]))).rjust(8),
                occ=str('%1.2f' % (float(atom[8]))).rjust(6),
                temp=str('%1.2f' % (float(atom[9]))).rjust(6),
                seg_id=atom[10].ljust(4),
                element=atom[11].rjust(2),
                charge='  '
            )

            new_pdb.write(updated_atom)
        for line in foot:
            new_pdb.write(line)
#TODO: Switch to new strucutred array method


def xyz_to_input(xyz: Union[str, np.ndarray], input_file: str, charge: Optional[str] = None, fixed: Optional[Union[int, List[int]]] = None, frame: int = 0) -> None:
    """
    Convert an XYZ file to an MBN Explorer input file.

    Parameters:
        xyz (str or numpy.ndarray): Path to an XYZ file or a structured array of coordinates.
        input_file (str): Path to the output MBN Explorer input file.
        charge (str, optional): Charge of the atoms, if specified. Default is None.
        fixed (int or list of int, optional): Index or indices of fixed blocks. Default is None.
        frame (int, optional): The simulation frame to use if multiple are specified. Defaults to the first frame.

    Returns:
        None: Writes the input file for MBN Explorer.
    """
    
    if type(xyz) == str:
        coords = read_xyz(xyz)[0]
    elif type(xyz) == np.ndarray:
        coords = xyz[0]
    
    with open(input_file, 'w') as f:
        for i, line in enumerate(coords):
            if i==fixed:
                f.write('<*\n')
            f.write(line['atoms']+(':'+charge if charge else '')+'\t\t\t'+'\t\t'.join(['{:.8e}'.format(float(x)) for x in line['coordinates']])+'\n')
        if fixed:
            f.write('>')
#TODO: Add dictionary of charges
#TODO: Change fixed to be lists. First item is start index, second is end index. Also list of lists. Single length list corresponds to fixed block after given index to end of file
