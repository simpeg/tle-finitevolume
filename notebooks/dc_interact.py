# Import numpy, python's n-dimensional array package,
# the mesh class with differential operators from SimPEG
# matplotlib, the basic python plotting package
import numpy as np
import discretize
from SimPEG import utils, Solver
import matplotlib.pyplot as plt


def dc_resistivity(
    log_sigma_background=1.,  # Conductivity of the background, S/m
    log_sigma_block=2,        # Conductivity of the block, S/m
    plot_type='potential'     # "conductivity", "potential", or "current"
):
    from pylab import rcParams
    rcParams['figure.figsize'] = 10, 10

    # Define a unit-cell mesh
    mesh = discretize.TensorMesh([100, 100])  # setup a mesh on which to solve

    # model parameters
    sigma_background = 10**log_sigma_background
    sigma_block = 10**log_sigma_block

    # add a block to our model
    x_block = np.r_[0.4, 0.6]
    y_block = np.r_[0.4, 0.6]

    # assign them on the mesh
    # create a physical property model
    sigma = sigma_background * np.ones(mesh.nC)

    block_indices = ((mesh.gridCC[:, 0] >= x_block[0]) &  # left boundary
                     (mesh.gridCC[:, 0] <= x_block[1]) &  # right boundary
                     (mesh.gridCC[:, 1] >= y_block[0]) &  # bottom boundary
                     (mesh.gridCC[:, 1] <= y_block[1]))   # top boundary

    # add the block to the physical property model
    sigma[block_indices] = sigma_block

    # Define a source
    a_loc, b_loc = np.r_[0.2, 0.5], np.r_[0.8, 0.5]
    source_locs = [a_loc, b_loc]

    # locate it on the mesh
    source_loc_inds = mesh.closest_points_index(source_locs)
    a_loc_mesh = mesh.gridCC[source_loc_inds[0], :]
    b_loc_mesh = mesh.gridCC[source_loc_inds[1], :]

    if plot_type == 'conductivity':
        plt.colorbar(mesh.plot_image(sigma)[0])
        plt.plot(a_loc_mesh[0], a_loc_mesh[1], 'wv', markersize=8)
        plt.plot(b_loc_mesh[0], b_loc_mesh[1], 'w^', markersize=8)
        plt.title('electrical conductivity, $\sigma$')
        return

    # Assemble and solve the DC resistivity problem
    Div = mesh.face_divergence
    Sigma = mesh.get_face_inner_product(sigma, invert_model=True, invert_matrix=True)
    Vol = utils.sdiag(mesh.cell_volumes)

    # assemble the system matrix
    A = Vol * Div * Sigma * Div.T * Vol

    # right hand side
    q = np.zeros(mesh.nC)
    q[source_loc_inds] = np.r_[+1, -1]

    # solve the DC resistivity problem
    Ainv = Solver(A)  # create a matrix that behaves like A inverse
    phi = Ainv * q

    if plot_type == 'potential':
        plt.colorbar(mesh.plot_image(phi)[0])
        plt.title('Electric Potential, $\phi$')
        return

    if plot_type == 'current':
        j = Sigma * mesh.face_divergence.T * utils.sdiag(mesh.cell_volumes) * phi
        plt.colorbar(mesh.plot_image(
            j,
            v_type='F',
            view='vec',
            stream_opts={'color': 'w'}
        )[0])
        plt.title('Current, $j$')
        return
