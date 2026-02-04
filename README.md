
# MBN Tools

**MBN Tools** provides a set of utilities for running, processing, and analysing simulations with [MBN Explorer](https://mbnresearch.com/get-mbn-explorer-software).

## Features

- **Run MBN Explorer simulations**  
  Run MBN Explorer simulation tasks directly from Python.

- **File handling utilities**  
  Create, read, modify, and convert simulation input/output files including `.task`, `.xyz`, `.dcd`, and`.pdb`.

- **Trajectory analysis**  
  Analyze simulation outputs, including radial distribution functions (RDF), mean squared displacement (MSD), RMSD, and more.

- **Crystallography utilities**  
  Group atoms by crystallographic planes, translate coordinates to planes, and remove atoms near simulation boundaries.

- **3D visualisation**  
  Render atoms, planes, and simulation boxes in 3D using `VisPy` with customisable colours, sizes, and markers.

## Installation

```bash
pip install MBN-tools
```

## Usage
```python
import MBN_tools as MBN
from MBN_tools import analysis, crystallography, visualise

# Run an MBN Explorer simulation
stdout, stderr = MBN.run_task("example.task", "/path/to/MBN_Explorer")

# Read an XYZ file
xyz_data = MBN.read_xyz("simulation.xyz")

# Compute RMSD
rmsd_values = analysis.calculate_rmsd(xyz_data, box_size=[50, 50, 100])

# Group atoms by a crystalline plane
planes = crystallography.group_atoms_by_plane(xyz_data, hkl=(1,1,0), tolerance=0.5)

# Visualize atoms in 3D
canvas = visualise.create_canvas()
visualise.draw_atoms(xyz_data['coordinates'], parent=canvas.central_widget.add_view())
```

## Dependencies

- `numpy >= 1.26`
- `scipy >= 1.12`
- `mdtraj`
- `vispy`