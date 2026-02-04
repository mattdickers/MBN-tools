"""
MBN_tools.visualise

A collection of 3D visualisation utilities for MBN Explorer simulations using VisPy.

This submodule provides tools to:
- Draw transparent simulation boxes with wireframe edges.
- Render arbitrary 3D planes with configurable orientation and border.
- Plot atom positions as markers with customisable appearance.

Dependencies:
- numpy
- vispy

Example:
    import MBN_tools as MBN
    from MBN_tools import visualise
    import vispy.scene

    # Create a canvas and a 3D view
    canvas = vispy.scene.SceneCanvas(keys='interactive', show=True, bgcolor='white')
    view = canvas.central_widget.add_view()

    # Read coordinates from .xyz file
    xyz_data = MBN.read_xyz('xyz_file.xyz')
    atom_coords = xyz_data['coordinates']

    # Draw a simulation box, a plane, and some atoms
    visualisation.draw_simulation_box(30, 30, 30, parent=view.scene)
    visualisation.draw_plane(centre=(0,0,0), size_x=10, size_y=10, normal=(0,0,1), parent=view.scene)
    visualisation.draw_atoms(atom_coords, face_color=(1,0,0,1), symbol='o', size=12, parent=view.scene)
"""

import numpy as np
import vispy.scene
from vispy.scene import visuals
from typing import Tuple, Optional, Union, Any

ColorType = Union[Tuple[float, float, float, float], str]

def draw_simulation_box(size_x: float, size_y: float, size_z: float,
                        edge_color: ColorType = (0, 0, 0, 1),
                        edge_width: int = 2,
                        parent: Optional[Any] = None) -> vispy.scene.visuals.Line:
    """
    Draw a transparent 3D box (rectangular cuboid) with visible wireframe edges.

    Parameters:
    - size_x (float): Size of the box in the x-direction.
    - size_y (float): Size of the box in the y-direction.
    - size_z (float): Size of the box in the z-direction.
    - edge_color (tuple): RGBA color of the box edges. Default is black.
    - edge_width (int): Width of the box edges. Default is 2.
    - parent (Node, optional): Parent VisPy node to which this box will be attached.

    Returns:
    - vispy.scene.visuals.Line: The wireframe visual representing the box.
    """
    # Define the 8 vertices of the box
    vertices = np.array([
        [-size_x / 2, -size_y / 2, -size_z / 2],
        [ size_x / 2, -size_y / 2, -size_z / 2],
        [ size_x / 2,  size_y / 2, -size_z / 2],
        [-size_x / 2,  size_y / 2, -size_z / 2],
        [-size_x / 2, -size_y / 2,  size_z / 2],
        [ size_x / 2, -size_y / 2,  size_z / 2],
        [ size_x / 2,  size_y / 2,  size_z / 2],
        [-size_x / 2,  size_y / 2,  size_z / 2],
    ])

    # Define the 12 edges of the box
    edges = np.array([
        [0, 1], [1, 2], [2, 3], [3, 0], # bottom face
        [4, 5], [5, 6], [6, 7], [7, 4], # top face
        [0, 4], [1, 5], [2, 6], [3, 7]  # vertical edges
    ])

    # Create the box edges (wireframe only)
    box_edges = vispy.scene.visuals.Line(pos=vertices, connect=edges, color=edge_color, 
                                   parent=parent, width=edge_width)

    # Set the box to be transparent by not adding a face visual
    # If you were using a filled mesh, you would set `face_color` to be transparent

    return box_edges


def draw_plane(centre: Tuple[float, float, float],
               size_x: float, size_y: float,
               normal: Tuple[float, float, float],
               color: ColorType = (1, 1, 1, 1),
               edge_color: ColorType = (0, 0, 0, 1),
               parent: Optional[Any] = None) -> Tuple[vispy.scene.visuals.Mesh, vispy.scene.visuals.Line]:

    """
    Draw a plane in 3D space with adjustable size in each direction, a specified normal vector, and a border.

    Parameters:
    - centre: A tuple (x, y, z) specifying the centre of the plane.
    - size_x: The size of the plane along the first basis vector direction.
    - size_y: The size of the plane along the second basis vector direction.
    - normal: A 3-element tuple specifying the normal vector of the plane.
    - color: The color of the plane with RGBA format (default is white).
    - edge_color: The color of the edges with RGBA format (default is black).
    - parent: The parent node to which this plane will be added (default is None).

    Returns:
    - A tuple containing the Mesh visual of the plane and the Line visual of the border.
    """

    # Normalize the normal vector
    normal = np.array(normal) / np.linalg.norm(normal)

    # Create an arbitrary vector that is not parallel to the normal
    if np.allclose(normal, [0, 0, 1]):
        basis1 = np.array([1, 0, 0])
    else:
        basis1 = np.cross(normal, [0, 0, 1])
        basis1 /= np.linalg.norm(basis1)  # Normalize basis1

    # Create the second basis vector
    basis2 = np.cross(normal, basis1)

    # Create the vertices of the plane using the local basis vectors
    half_size_x = size_x / 2
    half_size_y = size_y / 2
    vertices = np.array([
        -half_size_x * basis1 - half_size_y * basis2,
         half_size_x * basis1 - half_size_y * basis2,
         half_size_x * basis1 + half_size_y * basis2,
        -half_size_x * basis1 + half_size_y * basis2,
    ]) + np.array(centre)

    # Define the indices that connect the vertices to form two triangles
    indices = np.array([0, 1, 2, 0, 2, 3])

    # Create the plane as a Mesh visual
    plane = vispy.scene.visuals.Mesh(vertices=vertices, faces=indices.reshape((-1, 3)),
                               color=color, parent=parent)

    # Create the edges of the plane as a Line visual
    edges = vispy.scene.visuals.Line(pos=np.append(vertices, [vertices[0]], axis=0),  # Close the loop
                               color=edge_color, width=2, parent=parent, connect='strip')

    return plane, edges


def draw_atoms(atom_coordinates: np.ndarray,
               face_color: ColorType = (1, 1, 1, 1),
               edge_color: ColorType = (0, 0, 0, 1),
               symbol: str = 'o',
               edge_width: float = 1,
               size: float = 10,
               parent: Optional[Any] = None) -> vispy.scene.visuals.Markers:

    """
    Draw a collection of atoms as markers in 3D space.

    Parameters:
    - atom_coordinates (numpy.ndarray): Nx3 array of atomic coordinates.
    - face_color (tuple): RGBA color of marker faces. Default is white.
    - edge_color (tuple): RGBA color of marker edges. Default is black.
    - symbol (str): Marker shape (e.g. 'o' for circle, 's' for square).
    - edge_width (float): Width of marker edges.
    - size (float): Size of the markers.
    - parent (Node, optional): Parent VisPy node to attach the markers to.

    Returns:
    - vispy.scene.visuals.Markers: The visual representing the atoms.
    """
    atoms = visuals.Markers(parent=parent)
    atoms.set_data(atom_coordinates, face_color=face_color, edge_color=edge_color,
                   edge_width=edge_width, size=size, symbol=symbol)

    return atoms
