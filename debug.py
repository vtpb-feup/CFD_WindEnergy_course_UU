import numpy as np

def _interpolate_u(p,variable='v'):
    """
    Interpolates the input data to the U grid for staggered arrangement 
    calculations.

    Takes a 3D array `p`. Assumes a structured, uniform mesh.

    Parameters:
    ----------
    p : numpy.ndarray
        A 3D array with one of the u, v, w or p fields.

    variable : str, optional
        A string that specifies the variable for interpolation. It can take 
        the following values:
        - 'v': Interpolates for the v-component
        - 'w': Interpolates for the w-component
        - 'p': Interpolates for the pressure component by averaging in the 
               x direction only.
        Default is 'v'.

    Returns:
    -------
    numpy.ndarray
        A 3D array of interpolated values on the U grid, with the same shape 
        as `p` but only for inner nodes (reduced in size based on the 
        averaging operations).

    Notes:
    -----
    - The function uses numpy's `roll` function to shift the array for 
      averaging.
    - The output array will have dimensions reduced by 2 in the second 
      dimension due to the slicing operation `p[:,:, 1:-2]` and the averaging 
      process.
    """
    # interpolate to U grid (v_se, w_be, p_e)
    q=p[:,:, 1:-2] + p[:,:, 2:-1]       # averaging in x direction
    if variable=='v':
        q += np.roll(q,-1,axis=1)       # averaging in y direction
        q/=4
    elif variable=='w':
        q += np.roll(q,-1,axis=0)       # averaging in z direction
        q/=4
    elif variable=='p':
        q/=2

    return q

def init_monitor_points(xturb, yturb, zturb, rturb,
                        x_u, y_u, z_u, 
                        x_v, y_v, z_v,
                        x_w, y_w, z_w,
                        x_p, y_p, z_p):
    """
    Get coordinates of nearest neighbour nodes to monitor points.
    Monitor points are defined at the rotor centre, 2D downstream of each 
    rotor, 2D upstream of the first row and near the outlet.
    Does not allow for variable hub heights.

    Parameters:
    ----------
    xturb : list[float]
        A 1D list with the hub x-coordinates.

    yturb : list[float]
        A 1D list with the hub y-coordinates.

    zturb : list[float]
        A 1D list with the hub z-coordinates.

    rturb : list[float]
        A 1D list with the rotor radii.

    x_u, y_u, z_u : numpy.ndarray
        1D arrays with the x, y and z-coordinates of the u mesh nodes.

    x_v, y_v, z_v : numpy.ndarray
        1D arrays with the x, y and z-coordinates of the v mesh nodes.

    x_w, y_w, z_w : numpy.ndarray
        1D arrays with the x, y and z-coordinates of the w mesh nodes.

    Returns:
    -------
    coords : numpy.ndarray
        A 2D array with the coordinates of the monitoring points.

    nn_indx : numpy.ndarray
        A 3D array with the indices of the nearest neighbours to the 
        monitoring points for each variable's mesh (u, v, w, p). Dimension 0 
        is for each of the monitoring points, dimension 1 refers to the 3 
        coordinates and dimension 2 to the four variables.

    filepaths : list[str]
        List with filepaths for monitoring points.
    """

    ndim = 3
    nvar = 4
    downstream_dist = 4. * np.array(rturb)

    nturb = np.size(xturb)
    column_y = np.unique(np.array(yturb))
    nmonitor = int(2*nturb + np.size(column_y))
    coords = np.zeros((nmonitor, ndim))

    # assign hub locations
    coords[0:nturb,0] = np.array(xturb)
    coords[0:nturb,1] = np.array(yturb)
    coords[0:nturb,2] = np.array(zturb)

    # assign downstream locations
    coords[nturb:-np.size(column_y),0] = np.array(xturb) + downstream_dist
    coords[nturb:-np.size(column_y),1] = np.array(yturb)
    coords[nturb:-np.size(column_y),2] = np.array(zturb)
    # print(coords[:,0])
    # print(coords[:,1])
    # print(coords[:,2])

    # append also near outlet
    # (-5 due to p array size after interpolation to u mesh)
    coords[-np.size(column_y):,0] = np.array(x_p[-4])
    coords[-np.size(column_y):,1] = column_y
    coords[-np.size(column_y):,2] = np.array(zturb[-1])

    # Get nearest neighbour nodes for each variable
    nn_indx = np.zeros((nmonitor,ndim,nvar), dtype=int)
    for k in range(nmonitor):
        nn_indx[k,0,0] = int(np.argmin(np.abs(x_u - coords[k,0])))
        nn_indx[k,1,0] = int(np.argmin(np.abs(y_u - coords[k,1])))
        nn_indx[k,2,0] = int(np.argmin(np.abs(z_u - coords[k,2])))
        nn_indx[k,0,1] = int(np.argmin(np.abs(x_v - coords[k,0])))
        nn_indx[k,1,1] = int(np.argmin(np.abs(y_v - coords[k,1])))
        nn_indx[k,2,1] = int(np.argmin(np.abs(z_v - coords[k,2])))
        nn_indx[k,0,2] = int(np.argmin(np.abs(x_w - coords[k,0])))
        nn_indx[k,1,2] = int(np.argmin(np.abs(y_w - coords[k,1])))
        nn_indx[k,2,2] = int(np.argmin(np.abs(z_w - coords[k,2])))
        nn_indx[k,0,3] = int(np.argmin(np.abs(x_p - coords[k,0])))
        nn_indx[k,1,3] = int(np.argmin(np.abs(y_p - coords[k,1])))
        nn_indx[k,2,3] = int(np.argmin(np.abs(z_p - coords[k,2])))

    # nn_x = [int(np.argmin(np.abs(sim.x_u - xi))) for xi in x_conv]
    # nn_y = [int(np.argmin(np.abs(sim.y_u - yi))) for yi in y_conv]
    # nn_z = [int(np.argmin(np.abs(sim.z_u - zi))) for zi in z_conv]
    # print(nn_x)
    # print(nn_y)
    # print(nn_z)

    filepaths = init_convergence_file(coords)

    return coords, nn_indx, filepaths

def init_convergence_file(monitor_coords):
    """
    Initializes convergence files for monitoring points.

    This function creates a series of text files to log the outputs of a 
    simulation at set monitoring locations. Each file is named according to 
    the monitoring points' index and contains the coordinates along with 
    headers for the data that will be recorded.

    Parameters:
    ----------
    monitor_coords : numpy.ndarray
        A 2D array of shape (nmonitor, 3) where each row represents the 
        (x, y, z) coordinates of a monitoring point.

    Returns:
    -------
    list
        A list of file paths for the created convergence files. Each file is 
        named in the format 'conv_k.plt', where k is the index of the monitor 
        coordinate.
    """

    nmonitor, _ = np.shape(monitor_coords)
    filepaths = []
    for k in range(nmonitor):
        filepath = f'conv_{k}.plt'
        with open(filepath, 'w') as f:
            f.write('Convergence monitoring at:\n')
            f.write(f'{monitor_coords[k,0]}, '
                    f'{monitor_coords[k,1]}, '
                    f'{monitor_coords[k,2]}\n')
            f.write(f'i, time, dt, u, v, w, p, MAX_DIV\n')

        filepaths.append(filepath)

    return filepaths

def output_monitor(filepaths, nn_indx, 
                              iteration, time, dtime, u, v, w, p, max_div,
                              interpolate_to_u=True):
    """
    Outputs monitoring data to specified files.

    Parameters:
    ----------
    filepaths : list
        A list of file paths where the convergence data will be written. 
        Each file corresponds to a monitor point.

    nn_indx : numpy.ndarray
        A 3D array of shape (nmonitor, 3, 4) containing the indices for 
        accessing the variables (u, v, w, p) at each monitor point.

    iteration : int
        The current iteration number in the simulation.

    time : float
        The current simulation time.

    dtime : float
        The current time step used in the simulation.

    u : numpy.ndarray
        A 3D array representing the u-component values.

    v : numpy.ndarray
        A 3D array representing the v-component values.

    w : numpy.ndarray
        A 3D array representing the w-component values.

    p : numpy.ndarray
        A 3D array representing the pressure values.

    max_div : float
        The maximum divergence value to be recorded.

    interpolate_to_u : bool, optional
        A flag indicating whether to interpolate the v, w, and p values to 
        the U grid. 
        Default is True.

    Returns:
    -------
    None
    """

    nmonitor, _, _ = np.shape(nn_indx)

    if interpolate_to_u:
        for k in range(nmonitor):
            with open(filepaths[k], 'a') as f:
                line = f'{iteration}, {time}, {dtime}, ' \
                       f'{u[nn_indx[k,2,0],nn_indx[k,1,0],nn_indx[k,0,0]]}, ' \
                       f'{_interpolate_u(v, variable='v')[nn_indx[k,2,1],nn_indx[k,1,1],nn_indx[k,0,1]]}, ' \
                       f'{_interpolate_u(w, variable='w')[nn_indx[k,2,2],nn_indx[k,1,2],nn_indx[k,0,2]]}, ' \
                       f'{_interpolate_u(p, variable='p')[nn_indx[k,2,3],nn_indx[k,1,3],nn_indx[k,0,3]]}, ' \
                       f'{max_div}\n'
                f.write(line)
    else:
        for k in range(nmonitor):
            with open(filepaths[k], 'a') as f:
                line = f'{iteration}, {time}, {dtime}, ' \
                       f'{u[nn_indx[k,2,0],nn_indx[k,1,0],nn_indx[k,0,0]]}, ' \
                       f'{v[nn_indx[k,2,1],nn_indx[k,1,1],nn_indx[k,0,1]]}, ' \
                       f'{w[nn_indx[k,2,2],nn_indx[k,1,2],nn_indx[k,0,2]]}, ' \
                       f'{p[nn_indx[k,2,3],nn_indx[k,1,3],nn_indx[k,0,3]]}, ' \
                       f'{max_div}\n'
                f.write(line)

