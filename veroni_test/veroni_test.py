import pymatgen as mg
from pymatgen.analysis.defects.point_defects import Interstitial, ValenceIonicRadiusEvaluator
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
import operator
import numpy as np
import macrodensity as md
import math
import subprocess
import sys
import scipy.spatial.distance as spd



########################################################################################################################
# Author: Mitch Watts
# Date: 4-6-18
# References: some aspects of this code draw from the public github package 'twod_materials' intercalation code,
#             by author Michael Ashton (ashtonmv)
########################################################################################################################


########################################################################################################################
#TO DO
########################################################################################################################
# before generating combinations remove any sites in bad positions (within layers, outside the structure).
# remove any sites within layers
# remove any sites with wildly different potentials??
# remove any sites too close to the edges?

# add Ewald sum functionality
# sort by Ewald summation and remove any with wildly different values

# tidy up
# docstrings
# add potential of setup to info when structure is produced
# add user inputs/argparse so can select density of intercalant and the species, as well as filenames
########################################################################################################################

class Intercalate(object):


    def __init__(self, intercalant_type, intercalant_density, poscar_path = '', locpot_path = '', layered = True, only_sym_dist = False, coord_file_type = 'poscar', output_filepath = '',
                 output_filename = 'POSCAR_intercalated' , sens_fact = 0.97, vac_thresh = 6, default_vac = 3, voronoi_prox_crit = 1.5, voronoi_rad_cutoff_ratio = 0.3):
        self.int_atom = intercalant_type
        self.int_dens = intercalant_density
        self.layered = layered
        self.sym_dist = only_sym_dist
        self.coord_file_type = coord_file_type
        self.output_filepath = output_filepath
        self.atom_rad = Element(self.int_atom).atomic_radius
        self.sens_factor = sens_fact
        self.output_filename = output_filename
        self.vac_thresh = vac_thresh
        self.default_vac = default_vac
        self.voronoi_prox_crit = voronoi_prox_crit
        self.voronoi_rad_cutoff_ratio = voronoi_rad_cutoff_ratio

        self.poscar_path = poscar_path
        try:
            if self.poscar_path:
                self.structure = Structure.from_file(self.poscar_path)
            else:
                self.structure = Structure.from_file('POSCAR')
        except FileNotFoundError as e:
            print(e)
            print('The given filepath for the poscar did not find a poscar file, ensure POSCAR is at the end '
                  'e.g. /POSCAR. if no filepath was given then no poscar was found in the current directory.')

        self.no_ions = int(self.int_dens * self.structure.num_sites)

        self.locpot_path = locpot_path
        try:
            if self.locpot_path:
                self.vasp_pot, self.NGX, self.NGY, self.NGZ, self.lattice = md.read_vasp_density(self.locpot_path)
            else:
                self.vasp_pot, self.NGX, self.NGY, self.NGZ, self.lattice = md.read_vasp_density('LOCPOT')
            self.grid_pot, self.electrons = md.density_2_grid(self.vasp_pot, self.NGX, self.NGY, self.NGZ)  # generate the grid of potentials

        except FileNotFoundError as e:
            print(e)
            print('The given filepath for the locpot did not find a locpot file, ensure LOCPOT is at the end '
                  'e.g. /LOCPOT. if no filepath was given then no locpot was found in the current directory.')



        # then run as normal, once voronoi done and cells found, convert back and carry on
        if not self.check_orthog_lat():
            print('exiting prog')
            sys.exit()
        self.find_vacuums() #this changes structure to one without vacuum. copies the old structure as self.orig_struct



        self.evaluator = ValenceIonicRadiusEvaluator(self.structure) #computes site valence and ionic radii using bond valence analyser.
        # note this uses the sites, periodic table, bond valence, composition and local env packages! (as well as structure)
        self.radii = self.evaluator.radii
        self.valences = self.evaluator.valences


        self.interstitial = Interstitial(self.structure, radii = self.radii, valences = self.valences, symmetry_flag = self.sym_dist)
        # this evaluates the structure and uses radii and valence to generate voronoi sites depending on whether vertex,
        # facecenter or edge center is selected.
        # note: not sur eif hsould set oxi_state = False. By default it is true and it then uses ionic radii in the calculation,
        #  shouldnt really change centres? Returns and interstitial object. ._defect_sites returns the coordinates and labels of
        #  all intersittial sites which were found from the voronoi nodes/edges/faces/all depending on settings.
        # sym if False so we get all sites, including symmetry inequivalent sites!

        print('\n')
        self.interstitial_sites = [
            ([site._fcoords[0], site._fcoords[1], site._fcoords[2]], site.properties.get('voronoi_radius', None))
            for site in self.interstitial._defect_sites]  # this creates a list of tuples: [ (array of site location
        # fractional coordinates, Voronoi radius') ]
        #  replaced  with self.get_nearest_neighbour_dist(site)


        # shift vornoi up if shifted!
        self.interstitial_sites = [([i[0][0]*self.structure.lattice.a - self.downwards_shifts[0][0],
                                     i[0][1]*self.structure.lattice.b - self.downwards_shifts[1][1],
                                     i[0][2]*self.structure.lattice.c - self.downwards_shifts[2][2]],
                                    i[1]) for i in self.interstitial_sites] #shift it and convert to cart coords at the same time


        # go back to the original, unshifted structure
        self.shifted_structure = self.structure.copy()
        self.structure = self.orig_struct.copy()
        self.interstitial_sites = [([i[0][0]/self.structure.lattice.a, i[0][1]/self.structure.lattice.b,
                                     i[0][2]/self.structure.lattice.c], i[1]) for i in self.interstitial_sites] # convert it back to fractional coords

        print('removing any voronoi radii which are too small. the cutoff is anything less than: ', self.voronoi_rad_cutoff_ratio, 'times smaller than the ion size of: ', self.atom_rad)
        self.exclude_small_voronoi() # gets rid of any stupidly small spaces
        print('\n')

        print('removing and voronoi sites which are outside the bounds of the overall structure, i.e. theyre in the vacuum/on the surface')
        self.remove_surface_atoms() # get rid of any that aren't intercalated, and are just on the surfaces
        print('\n')

        # remove any sites within layers here

        # remove any sites with wildly different potentials??

        # remove any sites too close to the edges?



        print('sorting sites by potential/by voronoi radius')
        self.add_potentials()  # add the potentials and variance to the interstitial_sites list
        self.sort_interstitials() # sort them by voronoi radii, then by potentials if possible, else by ewald summation
        print('\n')


        print("adding the sites and all possible combinations of these sites to the structures")
        self.find_all_combos() # find all combinations of the interstitials found so far
        len_before_removal = len(self.structure_list) # how many elements before you remove th eoverlapping ones?
        self.remove_sites_too_close() # remove the overlapping points
        self.structs_removed = len_before_removal - len(self.structure_list)  # how many were removed

        # this needs to be redone to choose similar energy combinations of sites, rather than individual sites!
        #print('choosing the most energetic voronoi points which have the same potentials, within a factor of: ',
        #      self.sens_factor)
        # note you should only do this if we have the potentials!!!!
        #self.equivalent_interstitials()  # finds equivalent points based on their potentials
        #print('\n')
        self.calc_no_combos()

        self.add_interstitials() # adds the lowest energy sites to structure.
        # any potentials within e.g. 0.98 i.e. 98% of the maximum are included, and variations of all these structures are produced
        # the final structures are called: self.intercalated_structs

        # sort by Ewald summation here and remove any with widly different values

        self.report_expected_sites() # tell the user how many structures were made
        print('\n')

        print("writing all structures to files in a new folder called: ", self.output_filepath + 'intercalated_structures')
        self.write_to() # write it to a POSCAR

        self.write_report()

    #this might need to have atm_rad taken off. Also does distance to closest atom take account for periodicity?
    def get_nearest_neighbour_dist(self, site):
        distances = []
        atom_rad = 0
        for i in self.structure._sites:
            distances.append(i.distance_and_image_from_frac_coords(site._fcoords)[0])
            atom_rad = i.specie.atomic_radius

        return (min(distances) - atom_rad)


    def exclude_small_voronoi(self):
        self.interstitial_sites = [i for i in self.interstitial_sites if i[1] > self.atom_rad * self.voronoi_rad_cutoff_ratio] #this filters out anything smaller than factor* the ion radius


    def add_potentials(self):
        try:
            interstitials = []
            for i in self.interstitial_sites:
                origin = i[0]  # fractional coordinates of centre of sphere
                int_rad = i[1]  # radius of voronoi cell
                a = round(self.NGX * (int_rad / self.lattice[0, 0]))  # get radius of voronoi cell as fractional coordinates
                # in each axis, then get the no.grid points this corresponds to. Should be the same number for all three?
                # there will be a small roundign error here!
                b = round(self.NGY * (int_rad / self.lattice[1, 1]))
                c = round(self.NGZ * (int_rad / self.lattice[2, 2]))
                sphere_rad = [a, b, c]
                av_sphere, var_sphere = spherical_potential(self.grid_pot, sphere_rad, origin, self.NGX, self.NGY, self.NGZ)
                interstitials.append([i[0], i[1], av_sphere, var_sphere])
            self.interstitial_sites = interstitials

        except NameError as e:
            print('using Ewald summation method as no LOCPOT file has been loaded')
            return self.ewald_summation(self.structure)


    def sort_interstitials(self):

        # sort the sites by size of voronoi radius
        self.interstitial_sites.sort(key=operator.itemgetter(1))  # sort sites using the key operator.itemgetter(1).
        # This means: key(site) returns site[1], i.e. the voronoi radius and sorts by this in increasing order.
        self.interstitial_sites.reverse()  # reverse it to get in decreasing order.
        if self.vasp_pot.any():
            # now sort by the potentials
            self.interstitial_sites.sort(key=operator.itemgetter(2))  # sort sites using the key operator.itemgetter(2).
            # This means: key(site) returns site[2], i.e. the potential and sorts by this in increasing order.
            self.interstitial_sites.reverse()  # reverse it to get in decreasing order.
        else:
            print('using Ewald summation method as LOCPOT file has not been read correctly')
            return self.ewald_summation(self.structure)


    def ewald_summation(self, structure):
        print('Working on Ewald summation to sort sites. Not implimented yet so sites are sorted just by Voronoi radii')


    def equivalent_interstitials(self):

        self.equivalent_interstitials_lst = []
        print(self.interstitial_sites)
        for i in self.interstitial_sites:
            print(i)
            if i[2] >= self.interstitial_sites[0][2] * self.sens_factor:
                self.equivalent_interstitials_lst.append(i)


    def calc_no_combos(self):
        # no. of combinations of no_ions added is: len(equivalent_interstitials_lst) choose no_ions. i.e. ( n )
        #                                                                                                 ( k )
        # where n is len(equivalent_interstitials) and k is no_ions. this is n!/(k!(n-k)!)
        self.n_choose_k = self.n_choose(len(self.interstitial_sites),
                                            self.no_ions)  # no structures to make


    def n_choose(self, n, k):
        n_choose = (math.factorial(n)) / (math.factorial(k) * math.factorial(n - k))
        return n_choose


    def find_all_combos(self):
        ## use recurssion to find every combination of sites into subsets of size no_ions
        # put in external function
        def ion_combos(site_list, no_ions):
            site_combo_lst = []
            for i in range(len(site_list)):
                if no_ions == 1:
                    site_combo_lst.append(site_list[i])
                else:
                    for c in ion_combos(site_list[i + 1:], no_ions - 1):
                        site_combo_lst.append(site_list[i] + c)
            return site_combo_lst

        structure_list = ion_combos([[i[0]] for i in self.interstitial_sites], self.no_ions)

        self.structure_list = structure_list
        return structure_list

    def remove_sites_too_close(self):
        a = self.structure.lattice.a
        b = self.structure.lattice.b
        c = self.structure.lattice.c
        lat_vecs = [a, b, c]

        # remove overlapping combinations here!
        for i in range(len(self.structure_list)):
            sites = [[j[k] * lat_vecs[k] for k in range(3)] for j in
                     self.structure_list[i]]  # convert from fractional into cartesian

            dist_mat = spd.cdist(sites, sites,
                                 'euclidean')  # euclidean distances between all coodinate points. diagonals are zeros
            elements_too_close = dist_mat[np.where(
                dist_mat < self.voronoi_prox_crit)]  # create a matrix of values which are below the threshold closeness
            elements_too_close = elements_too_close[np.where(
                elements_too_close > 0)]  # if any of these aren't diagonals (where d = 0) which correspond to distance to themselves. then keep them in this matrix

            if elements_too_close.any():
                self.structure_list[i] = 'delete'  # mark for deletion

        self.structure_list = [i for i in self.structure_list if i != 'delete']  # remove elements marked for deletion


    def add_interstitials(self):
        # iterate through combinations and add them to pymatgen structures, so we end up with a list of all structures with combinations of sites
        #could put this in an external function to keep it tidy
        final_structs = []
        for i in self.structure_list:
            new_struct = self.structure.copy()
            for j in i:
                new_struct.append(self.int_atom, j)
            final_structs.append(new_struct)
        self.intercalated_structs = final_structs


    def report_expected_sites(self):
        if len(self.intercalated_structs) == int(self.n_choose_k):
            print('the number of structures was as expected, and is: ', int(self.n_choose_k))
        elif self.structs_removed:
            print('some structures were removed because the combinations of sites were not suitable. The number removed was: ', self.structs_removed)
            if len(self.intercalated_structs) + self.structs_removed == int(self.n_choose_k):
                print('the number of structures before removal was as expected, and was: ', str(len(self.intercalated_structs) + self.structs_removed))
                print('the number generated after removal was: ', len(self.intercalated_structs))
            else:
                print('the number generated after removal was: ', len(self.intercalated_structs))
                print('this differs from what was expected, which was: ', str(self.n_choose_k - self.structs_removed))
        else:
            print('there is a mismatch between the number of expected structures and the actual number')
            print('the number is expected was: ', self.n_choose_k)
            print('the number produced was: ', len(self.intercalated_structs))
        print('the number of structures was decided by the sensitivity factor, which was set at: ', self.sens_factor)


    def write_to(self):
        # write it to a poscar
        a = subprocess.Popen(['mkdir', self.output_filepath + 'intercalated_structures/'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)  # makes the directory, if its already there it is silent and error is stored in a
        a = a.communicate()

        for i in range(len(self.intercalated_structs)):
            filepath = self.output_filepath + 'intercalated_structures/' + self.output_filename + '_structure_no_' + str(i) + '.vasp'
            self.intercalated_structs[i].to(self.coord_file_type, filepath)


    def find_vacuums(self):

        a = self.structure.lattice.a
        b = self.structure.lattice.b
        c = self.structure.lattice.c
        lat_vecs = [a, b, c]
        vector = [0, 0, 0]
        sites_list = [i for i in range(self.structure.num_sites)]


        frac_coords = []
        for site in self.structure._sites:
            frac_coords.append(site._fcoords)

        self.orig_struct = self.structure.copy()
        self.downwards_shifts = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        self.lowest_atoms = [0, 0, 0]
        self.highest_atoms = [a, b, c]

        for i in range(3):
            min_frac_coord = min([j[i] for j in frac_coords])
            self.lowest_atoms[i] = min_frac_coord * lat_vecs[i]
            bot_vac = min_frac_coord * lat_vecs[i]
            max_frac_coord = max([j[i] for j in frac_coords])
            self.highest_atoms[i] = max_frac_coord * lat_vecs[i]
            top_vac = lat_vecs[i] - max_frac_coord*lat_vecs[i]

            if bot_vac >= self.vac_thresh:
                print('vacuum detected at the bottom of structure in lattice vector: ', str(i))
                print('vacuum size is: ', bot_vac)
                shift_vec = vector
                shift_vec[i] = 1
                shift_val = float(self.default_vac - bot_vac)
                shift_vec = [ i * shift_val for i in shift_vec]
                self.structure.translate_sites(sites_list, shift_vec, frac_coords = False)
                self.downwards_shifts[i] = shift_vec  # record the shift downward so you can move the voronoi vertices up

                lat_vecs[i] = lat_vecs[i] - bot_vac + self.default_vac

            if top_vac >= self.vac_thresh:
                print('vacuum detected at the top of structure in lattice vector: ' + str(i))
                print('vacuum size is: ', top_vac)

                lat_vecs[i] = lat_vecs[i] - top_vac + self.default_vac


        coords = [i._coords for i in self.structure._sites]
        #lat_vecs = [lat_vecs[0], lat_vecs[1], lat_vecs[2]]
        lat_vecs = [[lat_vecs[0], 0, 0],[0,lat_vecs[1], 0],[0, 0, lat_vecs[2]]]
        self.structure = Structure(lat_vecs, self.structure.species, coords, coords_are_cartesian = True)



    def remove_surface_atoms(self):
        a = self.structure.lattice.a
        b = self.structure.lattice.b
        c = self.structure.lattice.c
        lat_vecs = [a, b, c]
        self.interstitial_sites = [([i[0][j] for j in range(3) if ((i[0][j]*lat_vecs[j] < self.highest_atoms[j]) and (i[0][j]*lat_vecs[j] > self.lowest_atoms[j]))], i[1]) for i in self.interstitial_sites] #this will make a list only from coordinates that are between the bindings of atoms
        self.interstitial_sites = [i for i in self.interstitial_sites if len(i[0]) == 3] # the above list conc doesn't delete sites, just coordinates, so here delete any sites which have had any coordinates removed


    def remove_sites_near_edges(self):
        a = self.structure.lattice.a
        b = self.structure.lattice.b
        c = self.structure.lattice.c
        lat_vecs = [a, b, c]
        # this needs to be done layer by layer??



    def remove_sites_within_layers(self):
        return

    def find_layers(self):
        # return coordinates of all points in layer, plus bounding limits in each axis
        return


    def check_orthog_lat(self):
        alpha = int(self.structure.lattice.alpha)
        beta = int(self.structure.lattice.beta)
        gamma = int(self.structure.lattice.gamma)
        if alpha == 90 and beta ==90 and gamma == 90:
            return True
        else:
            print('structure is not orthogonal, dont have support for this yet')
            return False

    def write_report(self):
        return False

def spherical_potential(grid_pot, sphere_rad, origin, NGX, NGY, NGZ):

    #sphere_rad = [a, b, c]

    # Recalc the origin as grid point coordinates
    n_origin = np.zeros(shape=(3))
    n_origin[0] = int(origin[0]*NGX)
    n_origin[1] = int(origin[1]*NGY)
    n_origin[2] = int(origin[2]*NGZ)

    sphere_rad_a = int(sphere_rad[0])
    sphere_rad_b = int(sphere_rad[1])
    sphere_rad_c = int(sphere_rad[2])

    potential_sphere = np.zeros(shape=(sphere_rad_a*2, sphere_rad_b*2, sphere_rad_c*2))

    for x in range(-sphere_rad_a, sphere_rad_a):
        for y in range(-sphere_rad_b, sphere_rad_b):
            for z in range(-sphere_rad_c, sphere_rad_c):
                if (x**2)/(sphere_rad_a**2) + (y**2)/(sphere_rad_b**2) + (z**2)/(sphere_rad_c**2) <= 1:
                    # x^2/a^2 + y^2/b^2 + z^2/c^2 =1 is the sufrace of an ellipsoid (the radii are slightly different
                    #  because of VASP I think  -sampling NGX etc differently in each direction)
                    # Assign the values of coordinates in the original grid
                    xv = int(n_origin[0] + x)
                    yv = int(n_origin[1] + y)
                    zv = int(n_origin[2] + z)
                    # Minimum image convention
                    xv = int(xv - NGX * round(xv / NGX))
                    yv = int(yv - NGY * round(yv / NGY))
                    zv = int(zv - NGZ * round(zv / NGZ))

                    xind = sphere_rad_a + x
                    yind = sphere_rad_b + y
                    zind = sphere_rad_c + z

                    potential_sphere[xind, yind, zind] = grid_pot[int(xv), int(yv), int(zv)]

    return potential_sphere.mean(), np.var(potential_sphere) #how does having zeros outside of the sphere affect the mean??


    #return md.cuboid_average(grid_pot, cube, origin, travelled, NGX, NGY, NGZ, magnitude)








intercalated = Intercalate('Li', 0.2)