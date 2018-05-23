
# add check for which direction need reducing
def convert_to_com(array):
    """
    return an array which has converted a vector of points into a new vector of less points based on mean ring_positions
    Does this using KMeans algorithm from sklearn
    e.g. ||  ||  ||  ||   -> |  |  |  |  etc
    """

    if array.size < =2:
        return array
    else:
        cluster_number = math.ceil(array.siz e /2)

        kmeans = KMeans(n_clusters = cluster_number,  random_state=0).fit(array)  # the number of clusters here is half
        # the no. atoms along the length dimension
        converted = kmeans.cluster_centers_  # runs kmeans and gets the centers

        return converted


def simplify_2dsheet(array):  # , site, add_on_ends = True):
    """

    """

    new_positions = np.array([[0 ,0]])
    h_planes = np.unique(array[: ,0]) # unique values in z direction. i.e. the sites to iterate over

    for j in h_planes:
        row = np.asarray([i[1] for i in array if i[0] == j]) # output the row of x values for that z value
        row = row.reshape(-1, 1) # change format so works with Kmeans in convert to centre of mass function
        row = convert_to_com(row)  # reduce the points down from pairs to singles

        # build a new array with these rows
        row = np.insert(row, 0, j, axis=1) # add the z value to each x coordinate
        new_positions = np.append(new_positions, row, axis=0) # append the coordinates to the master list


    new_positions = np.delete(new_positions, 0, 0)  # delete the initialised start value

    return new_positions # this returns the simplified 2d sheet of atoms


def ring_centers(array, tol = 0.001, zw=1):

    x_edge1, x_edge2 = np.amin(array[: ,1]), np.amax(array[: ,1]) # extreme values. i.e. edges of cell
    print("edges: ")
    print(x_edge1, x_edge2)

    z_edge1, z_edge2 = np.amin(array[: ,0]), np.amax(array[: ,0])
    h_planes = np.unique(array[: ,0]) # unique values in z direction. i.e. the sites to iterate over

    positions = np.array([[0 ,0]]) # initialise positions. Delete this element before returning!


    for j in h_planes:

        row = np.asarray([i[1] for i in array if i[0] == j]) # get the row for this z value
        print("row")
        print(row)

        if not math.isclose(max(row), x_edge2, abs_tol = tol) and (z w -1):
            positions = np.append(positions, [[j, x_edge2]], axis=0) # if edge positions are empty add the coordinates of this point to positions
            print("edge. max: ")
            print(max(row))

        if not math.isclose(min(row), x_edge1, abs_tol = tol):
            positions = np.append(positions, [[j, 0]], axis=0)
            print("edge. min: ")
            print(min(row))

        else:
            print("no edge. row:")
            print(row)
            midpoints = np.asarray([(row[i ] -row[ i -1] ) / 2 +row[ i -1] for i in range(1 ,len(row))]) # gets all the interpolated midpoints
            print("midpoints")
            print(midpoints)
            midpoints = midpoints.reshape(-1, 1)

            midpoints = np.insert(midpoints, 0, j, axis=1)  # add the z value to each x coordinate


            positions = np.append(positions, midpoints, axis=0)  # append the coordinates to the master list


    positions = np.delete(positions, 0, 0) # delete the initialised start value

    return positions



def intercalate(atoms, layers, xwidth=1, zwidth=1):
    """
    finds the bottom of the layers, then iterates up and finds the sites to put alkali metal atoms into
    returns a np array with all atom positions, the end ones of which are the alkali metal positions
    """
    tol = 0.01 # tolerance for isclose

    int_atoms = np.array([[0 ,0 ,0]])
    temp_atoms = np.copy(atoms)

    for n in range(layer s -1): # iterate over all atoms
        print("n " +str(n))

        lower_atom_plane = np.amin(temp_atoms[: ,1]) # number which is y height in cell
        temp_atoms = temp_atoms[np.invert
            (np.isclose(temp_atoms[: ,1], lower_atom_plane, tol))]  # delete lowest plane so you can get second lowest

        upper_atom_plane = np.amin(temp_atoms[: ,1])

        half_plane_thickness = (upper_atom_plane - lower_atom_plane ) /2
        plane_centre = half_plane_thickness + lower_atom_plane

        temp_atoms = temp_atoms[np.invert
            (np.isclose(temp_atoms[: ,1], upper_atom_plane, tol))]  # delete lowest plane so you can get second lowest

        second_lower_plane = np.amin(temp_atoms[: ,1])

        nth_gap = second_lower_plane - upper_atom_plane

        # Don't forget to use this!
        y_value_for_int_atoms = upper_atom_plane + nth_ga p /2 # put halfway between layers


        temp_atoms2 = np.copy(atoms)
        mask = np.isclose(temp_atoms2[: ,1], plane_centre, atol=( half_plane_thickness + tol)) # make a vector of ones and zeros for rows corresponding to the layer of interest
        twod_plane_monolayer_atoms = atoms[mask]
        twod_plane_monolayer_atoms = twod_plane_monolayer_atoms[: ,[0 ,2]] # select the a and b axes

        # problem is in simplify function!!
        simplified_monolayer = simplify_2dsheet(twod_plane_monolayer_atoms)
        intercalant_atom_positions = ring_centers(simplified_monolayer)
        intercalant_atom_positions = np.insert(intercalant_atom_positions, 1, y_value_for_int_atoms, axis=1)  # add the z value to each x coordinate

        int_atoms = np.append(int_atoms, intercalant_atom_positions, axis=0)  # append the coordinates to the master list



    int_atoms = np.delete(int_atoms, 0, 0)  # delete the initialised start value


    return int_atoms



