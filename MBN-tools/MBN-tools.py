'''
MBN Tools

A series of functions for use with the MBN Explorer software (https://mbnresearch.com/get-mbn-explorer-software).

Functions include running MBN Explorer simulations, data analysis and file manipulation.
'''

import numpy as np
from scipy.spatial import KDTree
import warnings


def create_structured_array(atoms, frames):
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


def run_MBN(task_file, MBN_path, show_output=False):
    '''This function will run an MBN Explorer simulation using a specified Task file. Function returns the
stanard output and standard error. The show_output optional argument will print the simulation output to
the screen.'''
    import subprocess
    result = subprocess.run(MBN_path + ' -t ' + task_file, capture_output=True, text=True, creationflags=0x08000000)

    if show_output:
        print('Output:')
        if len(result.stdout) > 0:
            print(result.stdout)
        else:
            print(result.stderr)

    return result.stdout, result.stderr


def read_task(task_file, flatten=False):
    '''This function splits a given Task file into a a dictionary of task file options and their respective
parameters. Useful if you want to get a specific set of task file options for use in calculations for example.
by default it is split into the task file subsections as a dictionary of dictionaries. Use flatten=True to flatten
into a single dictionary.'''
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
    

def write_task(task_file, file_options):
    '''This function writes the contents nested dictionary containing task file options to a task file.
The input format is the same as a non-flattened read_task output dictionary.'''
    with open(task_file, 'w') as f:
        f.write(';\n')
        for key in file_options:
            f.write('\n; '+key+'\n')
            for sub_key in file_options[key]:
                if sub_key == 'Random':
                    f.write('\n'+'{:<31}'.format(sub_key)+'= '+file_options[key][sub_key]+'\n')
                else:
                    f.write('{:<31}'.format(sub_key)+'= '+file_options[key][sub_key]+'\n')


def read_xyz(xyz_file):
    '''This function returns the contents of an xyz file as a strucutred array.'''
    with open(xyz_file, 'r') as f:
        data = [i.split() for i in f.read().split('\n')[2:-1]]
        atoms = [i[0] for i in data]
        coordinates = [[[float(i[1]), float(i[2]), float(i[3])] for i in data]]

    xyz_data = create_structured_array(atoms, coordinates)

    return xyz_data


def write_xyz(coords, xyz_file, frame=0):
    '''This function writes an xyz file from a strucutred array of coordinates. Frame argument defines
the simulation frame to use of multiple are specified. Defaults to the first frame.'''
    with open(xyz_file, 'w') as f:
        f.write(str(len(coords[frame]))+'\n')
        f.write('Type name\t\t\tPosition X\t\t\tPosition Y\t\t\tPosition Z\n')
        for line in coords[frame]:
            f.write(line['atoms']+'\t\t\t\t'+'\t\t'.join(['{:.8e}'.format(float(x)) for x in line['coordinates']])+'\n')


def read_trajectory(dcd_file, xyz_file, frame=None):
    '''This function uses the mdtraj module to return the coordinates of a particular DCD file.
A frame number can be specified. Whether to include the atoms labels can be specified.
Note that this function is only compatible on unix systems due to the mdtraj module.'''
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
        
        
def xyz_to_pdb(xyz_file, pdb_file):
    '''This function is meant to update coordinates of a particular system from an XYZ file
to a corresponding PDB file for the same system. It receives an XYZ file (xyz_file)
containing the new coordinates and a PDB file (pdb_file) containing old coordinates to
be updated, and returns a PDB file with the name of the one supplied prefixed by 'new_'.
It assumes one is using only the record types 'ATOM' or 'HETATM' on the PDB file, but can
be adjusted for other usages provided the string updated_atom is correctly modified.'''
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


def xyz_to_input(xyz, input_file, charge=None, fixed=None, frame=0):
    '''This function converts an xyz file to an MBN explorer input file. xyz can either be a filename or a list
containing the atom type, and x, y, and z coordinates. Charge of the atom can be specified. This will be converted
to a dictionary in the future.'''
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


def melting_temperature_calculation():
    pass


def calculate_rdf(coordinates, step, r_max, frame=0, box_size=None, select_atoms=None):
    coordinates = coordinates[frame]
    if box_size is None:
        # Calculate box dimensions from given coors
        mins = np.min(coordinates['coordinates'], axis=0)
        maxs = np.max(coordinates['coordinates'], axis=0)
        box_size = maxs - mins
        warnings.warn('Using calculated box size. This is unadvisable and may produce incorrect results. For best results define a box size manually.')

    if select_atoms:
        # Reduce coordinates to selection of atom types
        coordinates = coordinates[np.where(coordinates['atoms'] == select_atoms)]
        if len(coordinates) == 0:
            raise ValueError('No atoms of the selected type.')

    num_atoms = len(coordinates)
    num_bins = int(r_max / step)
    
    rdf = np.zeros(num_bins)
    
    # Apply PBC for coordinates
    pbc_coordinates = coordinates['coordinates'] % box_size
    
    # Create KD-tree
    tree = KDTree(pbc_coordinates)
    
    # Average density rho for cuboidal box
    box_volume = np.prod(box_size)
    density = num_atoms / box_volume
    
    # Find pairs within max_distance using KD-tree
    pairs = tree.query_pairs(r=r_max, output_type='ndarray')
    
    # Calculate distances and bin them
    distances = np.linalg.norm(pbc_coordinates[pairs[:, 0]] - pbc_coordinates[pairs[:, 1]], axis=1)
    hist, bin_edges = np.histogram(distances, bins=num_bins, range=(0, r_max))
    
    r_values = (bin_edges[:-1] + bin_edges[1:]) / 2  # Mid points of bins
    
    # Convert histogram to RDF
    for i_bin in range(num_bins):
        r = r_values[i_bin]
        shell_volume = 4 * np.pi * r**2 * step
        rdf[i_bin] = hist[i_bin] / (shell_volume * density * num_atoms)

    return r_values, rdf


def rmsd_analysis():
    pass


def diffusion_analysis():
    pass


def kinetic_energy_distribution():
    pass


def temperatur_fluctuation_calculation():
    pass


def heat_capacity_calculation():
    pass


def spectral_statistics_analyser():
    pass