def runMBN(taskFile, MBN_path, show_output=False):
    import subprocess
    result = subprocess.run(MBN_path + ' -t ' + taskFile, capture_output=True, text=True, creationflags=0x08000000)

    if show_output:
        print('Output:')
        if len(result.stdout) > 0:
            print(result.stdout)
        else:
            print(result.stderr)

    return result.stdout, result.stderr





def readTaskFile(taskFile):
    fileOptions = {}
    with open(taskFile) as f:
        data = f.read().split('\n')

    for i in data:
        if '=' in i:
            i=i.split('=')
            try:
                fileOptions[i[0].strip()] = eval(i[1].strip())
            except (NameError, SyntaxError):
                fileOptions[i[0].strip()] = i[1].strip()

    return fileOptions






def readTrajectory(dcdFile, frame=None):
    loaded = False
    try:
        import mdtraj.formats as md
        loaded = True
    except:
        assert loaded, 'Cannot import the mdtraj module; please run in linux'

    with md.DCDTrajectoryFile(dcdFile) as f:
        xyz, cell_lengths, cell_angles = f.read()
        coords = (xyz)

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

