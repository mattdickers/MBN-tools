'''
MBN Tools

A series of functions for use with the MBN Explorer software (https://mbnresearch.com/get-mbn-explorer-software).

Functions include running MBN Explorer simulations, data analysis and file manipulation.
'''

import numpy as np

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


def read_xyz(xyz_file, skip_atoms=False, flatten=False):
    '''This function returns the contents of an xyz file as a list. skip_atoms excludes the atom labels 
from the list. Useful if you just want the raw coordinates but do not care about atom type. Flatten will
flatten the nested list into a single list of atom and coordinates of the form [['Atom', 'x', 'y', 'z'], ...].'''
    with open(xyz_file, 'r') as f:
        coords = [i.split() for i in f.read().split('\n')[2:-1]]
    
    if skip_atoms:
        for i in range(len(coords)):
            del coords[i][0]
            coords[i][0] = float(coords[i][0])
            coords[i][1] = float(coords[i][1])
            coords[i][2] = float(coords[i][2])
    else:
        new_coords = []
        for i in range(len(coords)):
            if flatten:
                new_coords.append([coords[i][0], float(coords[i][1]), float(coords[i][2]), float(coords[i][3])])
            else:
                new_coords.append([coords[i][0], [float(coords[i][1]), float(coords[i][2]), float(coords[i][3])]])
        coords = new_coords
    
    return coords


def write_xyz(coords, xyz_file, atoms=None):
    '''This function writes an xyz file from a list of coordinates. Optional atoms argument
defines a seperate list of atoms if coords containes only xyz coordinates and no atom labels.'''
    if atoms:
        try:
            coords = coords.tolist()
        except AttributeError: pass
        coords_new = []
        for i, atom in enumerate(atoms):
            combine = [atom, coords[i]]
            coords_new.append(combine)
        coords = coords_new
    
    with open(xyz_file, 'w') as f:
        f.write(str(len(coords))+'\n')
        f.write('Type name\t\t\tPosition X\t\t\tPosition Y\t\t\tPosition Z\n')
        for coord in coords:
            if len(coord)>2: #Check if the list has been flattened and if so, unflatten.
                coord = [coord[0], [coord[1], coord[2], coord[3]]]
            atom = coord[0]
            f.write(atom+'\t\t\t\t'+'\t\t'.join(['{:.8e}'.format(float(x)) for x in coord[1:][0]])+'\n')


def read_trajectory(dcd_file, frame=None):#, include_atoms=False):
    '''This function uses the mdtraj module to return the coordinates of a particular DCD file.
A frame number can be specified. Whether to include the atoms labels can be specified.
Note that this function is only compatible on unix systems due to the mdtraj module.'''
    loaded = False
    try:
        import mdtraj.formats as md
        loaded = True
    except:
        assert loaded, 'Cannot import the mdtraj module; please run in unix'

    with md.DCDTrajectoryFile(dcd_file) as f:
        xyz, cell_lengths, cell_angles = f.read()
        coords = (xyz)
    
##    if include_atoms:
##        atoms = [i[0] for i in read_xyz(dcd_file.replace('.dcd', '_dcd.xyz'))]
##        
##        for frame_coords in coords:
##            for i, coord in enumerate(frame_coords):
##                np.insert(coord, 0, atoms[i])
##            #print(i, frame_coords)

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


def xyz_to_input(xyz, input_file, charge=None, fixed=None):
    '''This function converts an xyz file to an MBN explorer input file. xyz can either be a filename or a list
containing the atom type, and x, y, and z coordinates. Charge of the atom can be specified. This will be converted
to a dictionary in the future.'''
    if type(xyz) == str:
        coords = read_xyz(xyz)
    elif type(xyz) == list:
        coords = xyz
    
    with open(input_file, 'w') as f:
        for i, coord in enumerate(coords):
            if i==fixed:
                f.write('<*\n')
            f.write(str(coord[0])+(':'+charge if charge else '')+'\t\t\t'+'{:.8e}'.format(coord[1][0])+'\t\t'+'{:.8e}'.format(coord[1][1])+'\t\t'+'{:.8e}'.format(coord[1][2])+'\n')
        if fixed:
            f.write('>')
#TODO: Add dictionary of charges
#TODO: Change fixed to be lists. First item is start index, second is end index. Also list of lists. Single length list corresponds to fixed block after given index to end of file


def melting_temperature_calculation():
    pass


def rdf_calculaton():
    pass


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