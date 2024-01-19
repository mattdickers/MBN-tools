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
