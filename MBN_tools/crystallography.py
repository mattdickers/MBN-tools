"""
MBN_tools.crystallography

A collection of crystallography-related utility functions for analysing
atomic structures produced by MBN Explorer simulations.

This module provides tools for:
- Grouping atoms relative to crystalline planes based on Miller indices.
- Translating atomic coordinates onto specified planes.
- Removing atoms within an exclusion region at the boundaries of a simulation box.

Dependencies:
- numpy

Example:
    import MBN_tools as MBN
    from MBN_tools import crystallography

    # Read coordinates from .xyz file
    xyz_data = MBN.read_xyz('xyz_file.xyz')
    atom_coords = xyz_data['coordinates']

    # Remove all atoms within a certain distance from simulation box edge
    filtered_atoms, removed_atoms = crystallography.remove_atoms_exclusion_region(atom_coords, box_size=30.0, shell_size=5.0)

    # Finds centre coordinate of all 110 crystalline planes
    planes = crystallography.group_atoms_by_plane(filtered_atoms, (1, 1, 0), tolerance=0.5)

    # Translate the cebtre point of planes to a the 001 normal
    for atoms, centre in planes:
        coords.append(centre)

    translated_coords = crystallography.translate_to_plane(coords, np.array([0, 0, 1]))
"""

import numpy as np
from typing import List, Tuple

def group_atoms_by_plane(structured_array: np.ndarray,
                         hkl: Tuple[int, int, int],
                         tolerance: float = 0.6) -> List[Tuple[np.ndarray, np.ndarray]]:

    """
    Group atoms by their position relative to a particular crystalline plane and calculate the centre of each plane,
    using a tolerance-based approach to more accurately group atoms to planes.
    
    Parameters:
        structured_array (numpy.ndarray): The structured array containing atom types and coordinates.
        hkl (tuple): The Miller indices (h, k, l) defining the crystalline plane.
        tolerance (float): Tolerance value for grouping atoms to planes (in the same units as lattice constant and coordinates).
        
    Returns:
        list: A list of tuples where each tuple contains an array of atoms on that plane and their centre coordinates.
    """
    
    # Unpack the Miller indices
    h, k, l = hkl
    
    # Define the plane normal
    plane_normal = np.array([h, k, l])
    
    # Normalize the plane normal
    norm_plane_normal = np.linalg.norm(plane_normal)
    plane_normal = plane_normal / norm_plane_normal
    
    # Project the atomic coordinates onto the plane normal to find their positions along this direction
    projected_positions = np.dot(structured_array['coordinates'], plane_normal)
    
    # Initialize a list to store groups of atoms for each plane
    planes = []
    
    # Sort the projected positions to process them in order
    sorted_indices = np.argsort(projected_positions)
    sorted_positions = projected_positions[sorted_indices]
    sorted_atoms = structured_array[sorted_indices]
    
    # Initialize the first plane
    current_plane_d = sorted_positions[0]
    current_plane_atoms = []
    
    # Group atoms into planes based on their position
    for pos, atom in zip(sorted_positions, sorted_atoms):
        # If the atom is within the tolerance of the current plane, add it to the plane
        if abs(pos - current_plane_d) <= tolerance:
            current_plane_atoms.append(atom)
        else:
            # Finalize the current plane and calculate its centre
            if current_plane_atoms:
                plane_atoms_array = np.array(current_plane_atoms, dtype=structured_array.dtype)
                centre_coordinates = np.nanmean(plane_atoms_array['coordinates'], axis=0)
                planes.append((plane_atoms_array, centre_coordinates))
            
            # Start a new plane
            current_plane_d = pos
            current_plane_atoms = [atom]
    
    # Finalize the last plane
    if current_plane_atoms:
        plane_atoms_array = np.array(current_plane_atoms, dtype=structured_array.dtype)
        centre_coordinates = np.nanmean(plane_atoms_array['coordinates'], axis=0)
        planes.append((plane_atoms_array, centre_coordinates))
    
    return planes


def translate_to_plane(coordinates: np.ndarray,
                       normal_vector: np.ndarray,
                       d: float = 0) -> np.ndarray:

    """
    Translate coordinates so that they lie on the specified plane.

    Parameters:
    - coordinates (numpy.ndarray): The atomic coordinates array, shape (N, 3).
    - normal_vector (numpy.ndarray): The normal vector of the plane.
    - d (float): The constant d in the plane equation ax + by + cz + d = 0. Default is 0.

    Returns:
    - translated_coordinates (numpy.ndarray): Coordinates translated to lie on the plane.
    """
    # Normalize the normal vector
    normal_vector = normal_vector / np.linalg.norm(normal_vector)
    
    # Calculate the distance from each coordinate to the plane
    distances = np.dot(coordinates, normal_vector) + d
    
    # Translate each coordinate by the distance along the normal vector
    translated_coordinates = coordinates - np.outer(distances, normal_vector)
    
    return translated_coordinates


def remove_atoms_exclusion_region(structured_array: np.ndarray,
                                  box_size: float,
                                  shell_size: float) -> Tuple[np.ndarray, np.ndarray]:

    """
    Remove atoms that are within a specified boundary distance from the edges of a simulation box.

    Parameters:
        structured_array (numpy structured array): A structured array containing atom data, including coordinates.
        box_size (float): The size of the simulation box (assumed cubic).
        shell_size (float): The distance from the edge of the box within which atoms should be removed.
    
    Returns:
        numpy structured array: A filtered structured array with atoms near the edges removed.
    """
    
        # Calculate the center of the crystal (assuming cubic box)
    centre = box_size / 2.0
    
    # Get the coordinates of all atoms
    coordinates = structured_array['coordinates']

    condition_remove = np.where(np.any(np.abs(coordinates)>centre - shell_size, axis=1))[0]
    
    # Atoms outside the boundary
    removed_atoms = structured_array[condition_remove]
    
    # Atoms inside the boundary
    filtered_atoms = structured_array[~np.any(np.abs(coordinates)>centre - shell_size, axis=1)]
    
    return filtered_atoms, removed_atoms