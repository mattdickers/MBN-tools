"""
MBN Tools

A collection of functions for use with the MBN Explorer software (https://mbnresearch.com/get-mbn-explorer-software).

This module includes functions for:
- Running MBN Explorer simulations.
- File manipulation.
- Data analysis.
- Visualisation.

Dependencies:
- numpy
- scipy
- mdtraj
- vispy

Example:
    import MBN_tools as MBN

    # Run MBN Explorer simulation
    stdout, stderr = MBN.run_MBN('task_file.task', '/path/to/MBN_Explorer')

    # Read XYZ file
    xyz_data = MBN.read_xyz('xyz_file.xyz')
"""


from .core import (
    create_structured_array,
    run_task,
    read_task,
    write_task,
    read_xyz,
    write_xyz,
    read_trajectory,
    xyz_to_pdb,
    xyz_to_input
)

from . import analysis
from . import crystallography
from . import visualise


__all__ = [
    # Core functions
    'create_structured_array',
    'run_task',
    'read_task',
    'write_task',
    'read_xyz',
    'write_xyz',
    'read_trajectory',
    'xyz_to_pdb',
    'xyz_to_input',

    # Submodules
    'analysis',
    'crystallography',
    'visualise',
]


def __dir__():
    return __all__