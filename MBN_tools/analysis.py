"""
Analysis submodule of MBN_tools.
Provides utility functions for MBN Explorer trajectory file analysis.

Functions:
- calculate_rdf
- calculate_msd
- calculate_rmsd

Example:
    import MBN_tools as MBN
    from MBN_tools import analysis

    xyz_data = MBN.read_xyz('xyz_file.xyz')
    RMSD = analysis.calculate_rmsd(xyz_data, box_size=[50, 50, 100], directions='xyz')
"""

import warnings
import numpy as np
from scipy.spatial import KDTree
from typing import Optional, Tuple

def calculate_rdf(coordinates: np.ndarray, step: float, r_max: float, frame: int = 0, box_size: Optional[np.ndarray] = None, select_atoms: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the radial distribution function (RDF) from coordinates.

    Parameters:
        coordinates (numpy.ndarray): A structured array of coordinates.
        step (float): The width of the bins for RDF calculation.
        r_max (float): The maximum distance for RDF calculation.
        frame (int, optional): The simulation frame to use if specified. Defaults to the first frame.
        box_size (numpy.ndarray, optional): Dimensions of the simulation box. If None, it will be calculated from coordinates.
        select_atoms (str, optional): Atom type to select. Default is None.

    Returns:
        tuple: (r_values, rdf) where r_values is an array of bin centers and rdf is the radial distribution function.
    """

    coordinates = coordinates[frame]
    if box_size is None:
        # Calculate box dimensions from given coors
        mins = np.min(coordinates['coordinates'], axis=0)
        maxs = np.max(coordinates['coordinates'], axis=0)
        box_size = maxs - mins
        warnings.warn('Using calculated box size. This is unadvisable and may produce incorrect results. For best results define a box size manually.')

    if select_atoms:
        # Reduce coordinates to selection of atom types
        coordinates = coordinates[np.isin(coordinates['atoms'], select_atoms)]
        if len(coordinates) == 0:
            raise ValueError('No atoms of the selected type(s).')

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


def calculate_msd(structured_array, box_size, directions='xyz'):
    """
    Calculate the Mean Squared Displacement (MSD).
    
    Parameters:
        structured_array (np.ndarray): A structured array of coordinates.
        box_size (float or list): Size of the simulation box. A single float assumes a cubic box, 
                                  and a list of three values specifies [Lx, Ly, Lz] for non-cubic boxes.
        directions (str): Directions to include in the MSD calculation ('x', 'y', 'z', or combinations like 'xyz', 'xy', etc.).
    
    Returns:
        np.ndarray: MSD values as a function of time.
    """
    # Extract coordinates
    coordinates = structured_array['coordinates']  # Shape: (n_timesteps, n_atoms, 3)
    n_timesteps, n_atoms, _ = coordinates.shape

    # Handle box size
    if isinstance(box_size, (int, float)):
        box_size = np.array([box_size] * 3)  # Convert to cubic box
    else:
        box_size = np.array(box_size)  # Ensure array format

    if box_size.shape != (3,):
        raise ValueError("Box size must be a single value for cubic boxes or a list of three values for non-cubic boxes.")

    # Calculate displacements
    displacements = np.zeros_like(coordinates)
    for t in range(1, n_timesteps):
        step_displacement = coordinates[t] - coordinates[t - 1]
        # Apply minimum image convention for PBC
        step_displacement -= box_size * np.round(step_displacement / box_size)
        displacements[t] = displacements[t - 1] + step_displacement

    # Initial positions
    initial_positions = displacements[0]  # Shape: (n_atoms, 3)

    # Total displacements
    total_displacements = displacements - initial_positions  # Shape: (n_timesteps, n_atoms, 3)

    # Filter for specified directions
    direction_indices = {'x': 0, 'y': 1, 'z': 2}
    selected_indices = [direction_indices[dim] for dim in directions]
    selected_displacements = total_displacements[..., selected_indices]  # Select desired dimensions

    # Squared displacements
    squared_displacements = np.sum(selected_displacements**2, axis=-1)  # Sum over selected dimensions

    # Mean squared displacement (MSD) over all atoms
    msd = np.mean(squared_displacements, axis=1)  # Average over particles

    return msd


def calculate_rmsd(structured_array, box_size, directions='xyz'):
    """
    Calculate the Root Mean Squared Displacement (RMSD),
    handling periodic boundary conditions (PBC).

    Parameters:
        structured_array (np.ndarray): A structured array of coordinates.
        box_size (float or list): Size of the simulation box. A single float assumes a cubic box,
                                  and a list of three values specifies [Lx, Ly, Lz] for non-cubic boxes.
        directions (str): Directions to include in the displacement calculation ('x', 'y', 'z', or combinations like 'xyz', 'xy', etc.).

    Returns:
        np.ndarray: RMSD values as a function of time.
    """
    # Compute MSD first
    msd = calculate_msd(structured_array, box_size, directions=directions)

    # Compute RMSD
    rmsd = np.sqrt(msd)

    return rmsd


def melting_temperature_calculation() -> None:
    """
    Placeholder function for melting temperature calculation.

    Returns:
        None
    """
    
    pass


def diffusion_analysis() -> None:
    """
    Placeholder function for diffusion analysis.

    Returns:
        None
    """
    
    pass


def kinetic_energy_distribution() -> None:
    """
    Placeholder function for kinetic energy distribution analysis.

    Returns:
        None
    """
    
    pass


def temperature_fluctuation_calculation() -> None:
    """
    Placeholder function for temperature fluctuation calculation.

    Returns:
        None
    """
    
    pass


def heat_capacity_calculation() -> None:
    """
    Placeholder function for heat capacity calculation.

    Returns:
        None
    """
    
    pass


def spectral_statistics_analyser() -> None:
    """
    Placeholder function for spectral statistics analysis.

    Returns:
        None
    """
    
    pass

