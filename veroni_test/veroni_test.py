import pymatgen as mg
import argparse
from pymatgen.analysis.defects.point_defects import Interstitial, ValenceIonicRadiusEvaluator
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.analysis.local_env import CrystalNN
from scipy.spatial import ConvexHull
from scipy.optimize import linprog
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
########################################################################################################################


########################################################################################################################
#TO DO
########################################################################################################################
# remove any sites with wildly different potentials??
# remove any sites too close to the edges?
# remove any structures with low potential sum

#write details of each structure energy to csv

# make some of the methods more general so they'll work for any material

# tidy up
## create extra functions where needed - esp inside __init__
## remove some functions from class
## add all printing tags somewhere else? i.e. one reporting function
## try if calling a function, try not to rebind an instance attribute inside that function, instead just return a value and bind the value of the function in the call -> makes functions pure where possible
## make the main class inherit all of Structure class, then don't need to call structure everytime, can just call self!

# docstrings
# add variance to the dictionary for each structure produced
# write energy info to a file

# improve argparse so that it can actually take some user input for other settings!
########################################################################################################################

class Settings():
     def __init__(self, intercalant_type, intercalant_density, poscar_path = '', locpot_path = '', layered = True, only_sym_dist = False, coord_file_type = 'poscar', output_filepath = '',
                 output_filename = 'POSCAR_intercalated' , vac_thresh = 6, default_vac = 3, voronoi_prox_crit = 1.8, voronoi_rad_cutoff_ratio = 0.2, delta = 1):

         self.int_atom = intercalant_type
         self.int_dens = intercalant_density
         self.layered = layered
         self.sym_dist = only_sym_dist
         self.coord_file_type = coord_file_type
         self.output_filepath = output_filepath
         self.atom_rad = Element(self.int_atom).atomic_radius
         self.output_filename = output_filename
         self.vac_thresh = vac_thresh
         self.default_vac = default_vac
         self.voronoi_prox_crit = voronoi_prox_crit
         self.voronoi_rad_cutoff_ratio = voronoi_rad_cutoff_ratio
         self.delta = delta  # for vdw bonding criteria

         self.poscar_path = poscar_path
         self.locpot_path = locpot_path


class Intercalate(object):


    def __init__(self, settings_obj,):
        self.settings_obj = settings_obj

        try:
            if self.settings_obj.poscar_path:
                self.structure = Structure.from_file(self.settings_obj.poscar_path)
            else:
                self.structure = Structure.from_file('POSCAR')
        except FileNotFoundError as e:
            print(e)
            print('The given filepath for the poscar did not find a poscar file, ensure POSCAR is at the end '
                  'e.g. /POSCAR. if no filepath was given then no poscar was found in the current directory.')

        self.no_ions = int(self.settings_obj.int_dens * self.structure.num_sites)

        try:
            if self.settings_obj.locpot_path:
                self.vasp_pot, self.NGX, self.NGY, self.NGZ, self.lattice = md.read_vasp_density(self.settings_obj.locpot_path)
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


        self.interstitial = Interstitial(self.structure, radii = self.radii, valences = self.valences, symmetry_flag = self.settings_obj.sym_dist)
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




    #this might need to have atm_rad taken off. Also does distance to closest atom take account for periodicity?
    def get_nearest_neighbour_dist(self, site):
        distances = []
        atom_rad = 0
        for i in self.structure._sites:
            distances.append(i.distance_and_image_from_frac_coords(site._fcoords)[0])
            atom_rad = i.specie.atomic_radius

        return (min(distances) - atom_rad)


    def exclude_small_voronoi(self):
        self.interstitial_sites = [i for i in self.interstitial_sites if i[1] > self.settings_obj.atom_rad * self.settings_obj.voronoi_rad_cutoff_ratio] #this filters out anything smaller than factor* the ion radius


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

    def sort_sites_by_pot_sum(self):

        # get the sum for each config
        for i in range(len(self.structure_list)):
            pot_sum = 0
            for j in self.structure_list[i]:
                pot_sum += j[2]
            self.structure_list[i].insert(0, {'sum_of_pots_from_locpot':pot_sum}) # insert it to the 0th element of the list of sites for that config

        # then sort by potential
        self.structure_list.sort(key = lambda k: k[0]['sum_of_pots_from_locpot'])
        self.structure_list.reverse()



    def ewald_summation(self, structure):
        ewald_energy = EwaldSummation(structure).total_energy
        return ewald_energy


    def oxidation_decorate(self, structure, oxidation_states):
        # decorates with charges to make suitable for ewald summation by changing oxidation states
        structure.add_oxidation_state_by_element(oxidation_states)
        return structure

    def ewald_all_structs(self):
        #self.intercalated_structs #this is the list of structures at index 0, dictionary at 1
        for i in range(len(self.intercalated_structs)):
            # self.intercalated_structs[i][0] is the structure.
            struct_to_sum = self.intercalated_structs[i][0].copy()
            frac_charge = - self.no_ions / (struct_to_sum.num_sites)
            oxidation_states = {self.settings_obj.int_atom : 1, str(struct_to_sum.species[0]):frac_charge} #this assumes the structure you're adding is monoelemental
            struct_to_sum = self.oxidation_decorate(struct_to_sum, oxidation_states)
            ewald_energy = self.ewald_summation(struct_to_sum)
            self.intercalated_structs[i][1]['Ewald_energy'] = ewald_energy # put this energy into the dictionary for safe-keeping


    def ewald_sorter(self):
        self.ewald_all_structs() #add an ewald energy to the dictionary for each structure
        #now sort by this ewald sum
        self.intercalated_structs.sort(key=lambda k: k[1]['Ewald_energy']) #sort in ascending order
        #self.intercalated_structs.reverse() #reverse so highest first



    ## this is no longer used, might be revived so keeping here for now
    def equivalent_interstitials(self):

        self.equivalent_interstitials_lst = []
        for i in self.interstitial_sites:
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
        def ion_combos(site_list, no_ions):
            site_combo_lst = []
            for i in range(len(site_list)):
                if no_ions == 1:
                    site_combo_lst.append(site_list[i])
                else:
                    for c in ion_combos(site_list[i + 1:], no_ions - 1):
                        site_combo_lst.append(site_list[i] + c)
            return site_combo_lst

        structure_list = ion_combos([[i] for i in self.interstitial_sites], self.no_ions)
        self.structure_list = structure_list

    def remove_sites_helper(self, sites1, sites2):

        # remove ions too close together
        dist_mat = spd.cdist(sites1, sites2,
                             'euclidean')  # euclidean distances between all coodinate points. diagonals are zeros
        elements_too_close = dist_mat[np.where(
            dist_mat < self.settings_obj.voronoi_prox_crit)]  # create a matrix of values which are below the threshold closeness
        return elements_too_close



    def remove_sites_too_close(self):
        a = self.structure.lattice.a
        b = self.structure.lattice.b
        c = self.structure.lattice.c
        lat_vecs = [a, b, c]

        for i in range(len(self.structure_list)):
            sites = [[j[0][k] * lat_vecs[k] for k in range(3)] for j in
                     self.structure_list[i]]  # convert from fractional into cartesian

            elements_too_close = self.remove_sites_helper(sites, sites)
            elements_too_close = elements_too_close[np.where(
                elements_too_close > 0)]  # if any of these aren't diagonals (where d = 0)
            # which correspond to distance to themselves. then keep them in this matrix

            if elements_too_close.any():
                self.structure_list[i] = 'delete'  # mark for deletion

            # do same again for ion - structure sites
            lattice_sites = []
            for site in self.structure._sites:
                lattice_sites.append(site._coords)
            elements_too_close = self.remove_sites_helper(sites, lattice_sites)
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
                if type(j) is list:
                    new_struct.append(self.settings_obj.int_atom, j[0])
            final_structs.append([new_struct, i[0]]) # now the details of the site are included next to each structure
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


    def write_to(self):
        # write it to a poscar
        a = subprocess.Popen(['mkdir', self.settings_obj.output_filepath + 'intercalated_structures/'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)  # makes the directory, if its already there it is silent and error is stored in a
        a = a.communicate()

        for i in range(len(self.intercalated_structs)):
            filepath = self.settings_obj.output_filepath + 'intercalated_structures/' + self.settings_obj.output_filename + '_structure_no_' + str(i) + '.vasp'
            self.intercalated_structs[i][0].to(self.settings_obj.coord_file_type, filepath)


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

            if bot_vac >= self.settings_obj.vac_thresh:
                print('vacuum detected at the bottom of structure in lattice vector: ', str(i))
                print('vacuum size is: ', bot_vac)
                shift_vec = vector
                shift_vec[i] = 1
                shift_val = float(self.settings_obj.default_vac - bot_vac)
                shift_vec = [ i * shift_val for i in shift_vec]
                self.structure.translate_sites(sites_list, shift_vec, frac_coords = False)
                self.downwards_shifts[i] = shift_vec  # record the shift downward so you can move the voronoi vertices up

                lat_vecs[i] = lat_vecs[i] - bot_vac + self.settings_obj.default_vac

            if top_vac >= self.settings_obj.vac_thresh:
                print('vacuum detected at the top of structure in lattice vector: ' + str(i))
                print('vacuum size is: ', top_vac)

                lat_vecs[i] = lat_vecs[i] - top_vac + self.settings_obj.default_vac


        coords = [i._coords for i in self.structure._sites]
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




    def remove_intralayer_atoms(self, structure, interstitial_sites):
        a = structure.lattice.a
        b = structure.lattice.b
        c = structure.lattice.c

        layers_dict = self.find_layers(structure) #returns a dictionary with layer numbers as keys, values as lists of coordinates

        #now if any sites in interstitial_sites lie within any layers, then we remove that site
        #make boxes from layer coordinates
        for i in layers_dict.keys():
            layers = layers_dict[i]
            #find bounding coordinates
            min_a = min(layers, key=lambda x: x[0])[0]
            max_a = max(layers, key=lambda x: x[0])[0]
            min_b = min(layers, key=lambda x: x[1])[1]
            max_b = max(layers, key=lambda x: x[1])[1]
            min_c = min(layers, key=lambda x: x[2])[2]
            max_c = max(layers, key=lambda x: x[2])[2]

            for j in range(len(interstitial_sites)):
                int_site = (
                interstitial_sites[j][0][0] * a, interstitial_sites[j][0][1] * b, interstitial_sites[j][0][2] * c)

                if min_a < int_site[0] < max_a:
                    if min_b < int_site[1] < max_b:
                        if min_c < int_site[2] < max_c:
                            #then its within the layer bounding box!
                            print('site found inside layer', i)
                            interstitial_sites[j] = 'delete'
            interstitial_sites = [i for i in interstitial_sites if i != 'delete']

        return interstitial_sites

    def find_layers(self, structure):
        # this only works if there is one species, otherwise it will assume every atom has vdw radius equal to the first sites'

        vdw_rad = Element(structure.species[0]).van_der_waals_radius
        cutoff = vdw_rad * 2 - self.settings_obj.delta
        struct_cpy = structure.copy()
        sites_list = [i._coords for i in struct_cpy._sites]
        layer_dict = self.atom_connectivity(sites_list, cutoff) # a dictionary where keys are the layer number and values are lists of coordinates of sites
        return layer_dict


    def atom_connectivity(self, site_list, cutoff):
        layers = []
        dist_mat = spd.cdist(site_list, site_list)

        # make a dictionary of which sites are connected to which sites
        for i in range(len(dist_mat)):
            connected = [i]
            for j in range(len(dist_mat[i])):
                if dist_mat[i][j] < cutoff and j != i:
                    connected.append(j)
            layers.append(tuple(connected))

        # now iterate through this and make a dictionary of layer number as key, value as a list of site indices
        layer_dict = {}
        for i in range(len(layers)): #maximum iterations is that there are no layers and every single atom is not bonded
            try:
                layer = list(layers[0])

                #return a list of all the values in tuples which have k in them
                def search(k, lst):
                    rtn_lst = []
                    for i in range(len(lst)):
                        try:
                            if k in lst[i]:
                                for j in lst[i]:
                                    rtn_lst.append(j)
                                lst[i] = 'deleted'
                        except TypeError:
                            pass
                    return (rtn_lst, lst)

                # iterate over every tuple, iterate over its values and pull out all values from any other tuple containing one of the numbers uou're currently iterating over. only if its in the list 'layer'
                # i.e. layer contains all the atoms which are in the same layer as the seed position which is layers[0] for every iteration
                for j in range(len(layers)):
                    for k in layers[j]:
                        if k in layer:
                            layer_add, layers = search(k, layers)
                            layer += layer_add
                    layer = list(set(layer))  # to get unique values only

                layers = [i for i in layers if i != 'deleted']  # deletes tuples which were linked
                layer_dict[i] = layer #key is the layer number, i.e. the iteration are in


            except IndexError as e:
                break

        #now convert site indices into coordinates
        for i in layer_dict.keys():
            for j in range(len(layer_dict[i])):
                k = layer_dict[i][j]
                layer_dict[i][j] = site_list[k]
        return layer_dict

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



def intercalation_routine(metal, density):

    settings = Settings(metal, density)

    intercalated = Intercalate(settings)


    print(
    'removing any voronoi radii which are too small. the cutoff is anything less than: ', settings.voronoi_rad_cutoff_ratio,
    'times the ion size of: ', settings.atom_rad)
    intercalated.exclude_small_voronoi()  # gets rid of any stupidly small spaces
    print('\n')

    print(
        'removing any voronoi sites which are outside the bounds of the overall structure, i.e. theyre in the vacuum/on the surface')
    intercalated.remove_surface_atoms()  # get rid of any that aren't intercalated, and are just on the surfaces
    print('\n')

    print('removing any voronoi sites which are within layers of the structure')
    intercalated.interstitial_sites = intercalated.remove_intralayer_atoms(intercalated.structure,
                                                                           intercalated.interstitial_sites)  # get rid of any that aren't intercalated, and are just on the surfaces
    print('\n')

    print('sorting sites by potential/by voronoi radius')
    intercalated.add_potentials()  # add the potentials and variance to the interstitial_sites list
    intercalated.sort_interstitials()  # sort them by voronoi radii, then by potentials if possible, else by ewald summation
    print('\n')

    print("adding the sites and all possible combinations of these sites to the structures")
    intercalated.find_all_combos()  # find all combinations of the interstitials found so far
    len_before_removal = len(intercalated.structure_list)  # how many elements before you remove th eoverlapping ones?
    intercalated.remove_sites_too_close()  # remove the overlapping points and points to close to structure atoms
    intercalated.structs_removed = len_before_removal - len(intercalated.structure_list)  # how many were removed

    # this needs to be redone to choose similar energy combinations of sites, rather than individual sites!
    # print('choosing the most energetic voronoi points which have the same potentials, within a factor of: ',
    #      self.sens_factor)
    # note you should only do this if we have the potentials!!!!
    # self.equivalent_interstitials()  # finds equivalent points based on their potentials
    print('\n')
    intercalated.calc_no_combos()  # get the number of combinations

    intercalated.sort_sites_by_pot_sum()  # sort the found interstitial combos by their potential sum, largest to smallest
    # self.structure_list = self.structure_list[0:10] # pick the top in terms of potential points

    intercalated.add_interstitials()  # adds the sites to structure, most energetically favourable first
    # the final structures are called: self.intercalated_structs

    # sort by Ewald summation here and remove any with widly different values
    # note: self.intercalated_structs elements 0 are dictionaries containing energies etc.  'sum_of_pots_from_locpot' for sum of pot energy for structure

    ##### this doesn't take into accoutn periodic images!!
    intercalated.ewald_sorter()  # this will add the ewald sum energy for this configuration to the dictionary in position 1 in each list under 'Ewald_energy' key

    intercalated.report_expected_sites()  # tell the user how many structures were made
    print('\n')

    print("writing all structures to files in a new folder called: ", settings.output_filepath + 'intercalated_structures')
    intercalated.write_to()  # write it to a POSCAR

    intercalated.write_report()







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Intercalate a vasp structure')
    parser.add_argument('details', nargs = '*', help = 'enter first arg as metal e.g Li, second as a number which is density, e.g. 0.1')
    ns = parser.parse_args()

    metal = str(ns.details[0]) #first arg is the metal e.g. Li
    density = float(ns.details[1]) #second is density e.g. 0.1

    intercalation_routine(metal, density) #call routine to produce intercalated structures
