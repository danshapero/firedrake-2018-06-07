
import numpy as np
from numpy.linalg import norm
import firedrake

def minimum_angle(mesh):
    mesh.init()
    coords = mesh.coordinates.dat.data_ro
    cells = mesh.coordinates.cell_node_map().values

    if not ((mesh.cell_dimension() == 2) and
            mesh.is_piecewise_linear_simplex_domain()):
        raise ValueError("Only works on 2D triangular mesh!")

    angle = np.inf
    for cell in cells:
        for k in range(3):
            x, y, z = coords[np.roll(cell, k)]
            u, v = y - x, z - x
            angle = min(angle, np.arccos(np.inner(u, v) / (norm(u) * norm(v))))

    return angle
