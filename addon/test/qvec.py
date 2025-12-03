import itertools
import numpy as np
from numpy.typing import NDArray
from dynasor.logging_tools import logger


def get_spherical_qpoints(
        cell: NDArray[float],
        q_value: float,
) -> NDArray[float]:
    r"""Generates all q-points on the reciprocal lattice inside a given radius
    :attr:`q_max`.  This approach is suitable if an isotropic sampling of
    q-space is desired.  The function returns the resulting q-points in
    Cartesian coordinates as an ``Nx3`` array.

    If the number of generated q-points are large, points can be removed by
    specifying the :attr:`max_points`. The q-points will be randomly removed in
    such a way that the q-points inside are roughly uniformly distributed with
    respect to :math:`|q|`. If the number of q-points are binned w.r.t. their
    norm the function would increase quadratically up until some distance P
    from which point the distribution would be constant.

    Parameters
    ----------
    cell
        real cell with cell vectors as rows
    q_max
        maximum norm of generated q-points (in units of rad/Å, i.e. including factor of 2pi)
    max_points
        Optionally limit the set to __approximately__ :attr:`max_points` points
        by randomly removing points from a "fully populated mesh". The points
        are removed in such a way that for :math:`q > q_\mathrm{prune}`, the
        points will be radially uniformly distributed. The value of
        :math:`q_\mathrm{prune}` is calculated from :attr:`max_q`,
        :attr:`max_points`, and the shape of the cell.
    seed
        Seed used for stochastic pruning

    """
    q_min = 0.95 * q_value
    q_max = 1.05 * q_value

    # inv(A.T) == inv(A).T
    # The physicists reciprocal cell
    rec_cell = np.linalg.inv(cell.T) * 2 * np.pi

    # We want to find all points on the lattice defined by the reciprocal cell
    # such that all points within max_q are in this set
    inv_rec_cell = np.linalg.inv(rec_cell.T)  # cell / 2pi

    # h is the height of the rec_cell perpendicular to the other two vectors
    h = 1 / np.linalg.norm(inv_rec_cell, axis=1)

    # If a q_point has a coordinate larger than this number it must be further away than q_max
    N = np.ceil(q_max / h).astype(int)

    # Create all q-points within a sphere
    lattice_points = list(itertools.product(*[range(-n, n+1) for n in N]))
    q_points = lattice_points @ rec_cell

    # Calculate distances for pruning
    q_distances = np.linalg.norm(q_points, axis=1)  # Find distances

    # Prune based on distances
    q_points = q_points[(q_distances >= q_min) & (q_distances <= q_max)]
    q_distances = q_distances[(q_distances >= q_min) & (q_distances <= q_max)]

    # keep only 50 central q-points if too many
    if q_points.shape[0] > 50:
        q_points = q_points[q_points.shape[0]//2-25:q_points.shape[0]//2+25]
        
    return q_points

