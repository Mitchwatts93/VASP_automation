import subprocess       # used for cp, rm, mkdir
import os               # only used for cwd
import datetime         # used in log files
import re               # used for easier reading in of data for analysis sake
import csv              # used for analysis output file
import argparse         # used for cmd line input




##########################################################################################
# to do
##########################################################################################
# - add docstrings and doctests
# - clean it up and add some abstraction barriers
# - add more info to the output of analysis document e.g. run time and scf steps
# - swap trees to class implementation
# - fix re-intialising trees etc in new classes, just refer to the originals! instead initialise
# e.g. the settings object
# - test running from path, do directory paths etc work or need editing with cwd?
# - The file directories to create on the fly have to be named either E_cutoff or Kpts. Fix this.
###########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
#Format of conv_params file
##########################################################################################
##########################################################################################

# Ensure that INCAR has ENCUT in it, and KPOINTS has Gamma in it, just before where KPOINTS should be added.

#note that this file should be names conv_params and be located in your INPUT folder
#Your INPUT folder should contian all files you wish to copy. This script will make directories as specified by
# the format below, it will edit releavnt files but won't touch others, e.g. your job files won't be edited

#E_low_range E_multiple E_higg_range    e.g. 200 50 300 will make 200, 250 and 300 values
#list of space separated kpoints which are dash delimited. e.g. 4-4-4 3-3-3 2-2-2
#default energy cutoff as a number
#default Kpts as k-k-k

##########################################################################################
##########################################################################################


##########################################################################################
# note on subprocess. This may not be implemented well. Have not fully read documentation.
# using shell = True can be a security hazard due to shell injection esp if command string is constructed from external
# input (not really an issue here but not used for good practice). Only placed used is for cp input files. This
# may be an issue if you can somehow get a string in there but I think its okay from external input

# further note: uisng pipe to capture output may cause issues if pipe buffer is filled, in which case it will cause the program to hang.
# remove stderr/stdout from subprocess calls if this happens - this will cause bash outputs to be seen during running and will not allow details to be saved to log

# documentation suggests using Popen with communicate() method when pipes are needed -> this seems to be after already using pipe
# e.g. a.communicate() is the same as a.stderr.read(). not sure if this will cause an issue though...
##########################################################################################




##########################################################################################
##########################################################################################
##########################################################################################
# program START #
##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
# abstraction functions
##########################################################################################

# This can all be put into a class, maybe comnbined with the current class for analysis reading files

def file_reader(filename = 'conv_params'):
    """
    Reads a conv_params.txt file parameters and returns them as a list.


    Args:
        filename: a string object which is the name of the file with settings in it. 'INPUT/conv_params' by default

    Returns:
          a list of settings, each element is a line in file:
            0: the energies which will be tested as strings. first is lower energy bound, second is the step size
                third is the upper range, e.g. [['450'],['50'],['600']]
            1: the kpoints to be tested as a list, each part of the string is an element in a sublist. e.g.
                [['7','1','1'],['8','1','1']].
            2: the default energy as a string e.g. ['400']
            3: the default kpoints as a list e.g. ['7','1','1']

    Raises:
        Currently nothing.

    To do:
        Add checks for format and raise errors here if structure is incorrect:
            - first line must be three space seperated numbers (integers?)
            - second line must be any amount of space seperated dash delimited kpoints
            - a single number
            - a single set of dash delimited kpoints
    """

    with open(filename) as f:
        info = f.readlines()
    info = [i.strip() for i in info] # each element in info is a string of a line from the file
    info = [i.split() for i in info] # split each whitespace delimited element into a list of lists
    info = [[i.split('-') for i in j] for j in info] # note info is 3 layers deep

    info[2] = info[2][0] # makes default E just a single string of the number
    info[3] = info[3][0]

    return info

def get_no_atoms(filename = 'INPUT/POSCAR'):

    with open(filename) as f:
        lines = f.readlines()
    lines = [i.strip() for i in lines] # each element in lines is a string of a line from the file
    lines = lines[5:7] #only interested in lines 6 and 7 in the file -  (python 0 indexing)

    no_atoms = 0
    for i in lines:
        i = ''.join([j for j in i if j.isdigit()])
        if i: #i.e. if there are any numbers in the line!
            no_atoms = int(i)

    if not no_atoms:
        no_atoms = 1
        print('COULDNT FIND THE NUMBER OF ATOMS!!')

    return no_atoms


def edit_kpts(param_label, i, dir, line_key = 'Gamma', file = 'KPOINTS'):
    """
    Edits a KPOINTS file to have kpoints equal to those passed in.

    Args:
        param_label: the name of the parent directory as a string. e.g. 'Kpts'
        i: a list of kpoints of length 3 where each element is a kpoint to be written
        dir: the name of the parent directory, currently equal to the kpoints to be written e.g. '7-1-1'

    Returns:
        False: using this convention as when add assertion errors then it will return a value -> True if anything goes
            wrong

    Raises:
        Currently nothing

    To do:
        Assertion errors. Note currently this only works if 'Gamma' is in the kpoints file, and the kpoints to be edited
            are directly below this. Note this will currently work if you have 'Gamma' commented out above the kpoints
            to be edited so its not an issue.
        Can probably get rid of dir argument as you can generate it by: i.join('-'), but allows flexibility elsewhere?
    """

    replacement_line = "  " + i[0] + "  " + i[1] + "  " + i[2]
    gen_file_editor(param_label, dir, file, replacement_line, line_key)

    return False


def edit_incar(param_label, i, dir, line_key = 'ENCUT', file = 'INCAR'):
    """
    Edits an INCAR file according to inputs to have the correct ENCUT.

    Args:
        param_label: the name of the parent directory as a string. e.g. 'E_cutoff'
        i: a list of lists of energies length 1, e.g. [['400'],['500']]
        dir: the name of the parent directory, currently equal to the energy to be written e.g. '400'

    Returns:
        False: using this convention as when add assertion errors then it will return a value -> True if anything goes
            wrong

    Raises:
        Currently nothing

    To do:
        Assertion errors. Note currently this only works if 'Gamma' is in the kpoints file, and the kpoints to be edited
            are directly below this. Note this will currently work if you have 'Gamma' commented out above the kpoints
            to be edited so its not an issue.
        Can probably get rid of dir argument as you can generate it by: i.join('-'), but allows flexibility elsewhere?

    """

    replacement_line = "  ENCUT   = " + i[0] + "       ! Plane-wave cutoff"
    gen_file_editor(param_label, dir, file, replacement_line, line_key)

    return False


def gen_file_editor(param_label, i, file, replacement_line, line_key):
    """
    A general file editor that uses inputs passed from either the edit_kpts or edit_incar functions.

    Args:
        param_label: the name of the parent directory as a string. e.g. 'E_cutoff'
        i: the name of the parent directory, currently equal to the thing to be written e.g '400' or '7-7-1'
        file: the name of the file to edit e.g. 'KPOINTS' or 'INCAR'
        replacement_line: the line to be inserted
        line_key: the key for where to input the replacement line. This should be changed to be a line no. and the
            key used in either edit_kpts or edit_incar decides the number...


    Returns:
        False: using this convention as when add assertion errors then it will return a value -> True if anything goes
            wrong

    Raises:
        Currently nothing

    To do:
        by this point everything passed in should have been checked by previous functions. If it hasn't been passed in
            you may have issues so do a check on that somehow with some dummy argument? Then if it is called from
            somewhere else do checks on the format of everything passed in.
    """

    file = param_label + '/' + "-".join(root(i)) + '/' + file  # path to file
    lines = open(file).read().splitlines()  # read lines into elements in a list where each element is a new line in the
        # file

    # this could be moved into another function?
    # find where to replace the line based on what we are finding.
    line_no = False
    for i in lines:
        if line_key in i:
            line_no = lines.index(i)

    # if the line key isn't there
    if not line_no:
        return 'Line key not found in files. Ensure that INCAR has ENCUT in it, and KPOINTS has Gamma in it, just before' \
               ' where KPOINTS should be added'

    # for kpoints we find the 'Gamma' string and replace the line below with the kpoints we want to write
    if line_key == 'Gamma':
        line_no += 1

    # replace the
    lines[line_no] = replacement_line  # replace the kpoints details in the list
    open(file, 'w').write('\n'.join(lines))  # write list to file

    return False


def get_params(path = 'INPUT/conv_params'):
    """
    A function which calls other functions in order to get the user parameters, clean them and make an object from them

    Args:
        path: path to the conv_params user settings. default is 'INPUT/conv_params'

    Returns:
        a Settings object with the settings of the users choosing
    """

    # cd to Input, read the conv_params file in and pass each line to file reader
    list = file_reader(path)

    Ecuts =  list[0] # first element returned from filereader is the energies
    start = int(Ecuts[0][0]) # the first element of this is the lower energy to start from. convert to integer for maths
    multiplier = int(Ecuts[1][0]) # middle element is the step size
    end = int(Ecuts[2][0]) # last element is upper bound on energy
    E_range = (end - start)//multiplier +1 # the number of energies you will create
    Es = [i*multiplier for i in range(E_range)] # take steps in the E_range of step size multiplier
    Ecuts = [[str(i+start)] for i in Es] # add the start energy to all these steps to shift them to correct energies
                                            # convert the numbers to strings for ease of file writing later

    kpts = list[1] # kpoints list is first element returned
    def_E = list[2] # default energy
    def_k = list[3] # default kpoints
    params = Settings(Ecuts, kpts, def_E, def_k) # create the settings object

    return params # return the object



def remove_fly_bodge(path):
    path_fly = path +'fly/'
    subprocess.call("mv "+path_fly+"* "+path, shell=True)
    subprocess.call("rm -r "+path_fly, shell=True)

    return False



##########################################################################################
# tree stuff - put this into a class and change all calls to trees appropriately
##########################################################################################

# Class Tree:
#    """
    #You need to implement this!
#   """


def tree(root, branches = []):
    """
    A tree function which builds a tree object from a root label and a list of branches

    Args:
        root: in this application either 'Ecuts' or 'Kpts', but is supposed to be the parent label
        branches: A list of the elements to be made, e.g. the list of energies you are testing

    Returns:
        A tree, which is a list with first element a label for the tree, second element is a list of branches which are
            themselves trees whose roots are the sub directories to be created

    Raises:
        If the branches aren't themselves trees then it raises an error that they must be trees themselves.
    """

    for branch in branches:
        assert is_tree(branch), 'branches must be trees'
    return [root] + list(branches)

def root(tree):
    """
    Returns the root of a tree.

    Args:
        the tree

    Returns:
        the root label of a tree, i.e. the first element of a list.

    Raises:
        Nothing

    """

    return tree[0]

def branches(tree):
    """
    Returns the branches of a tree.

    Args:
        the tree

    Returns:
        the branches of a tree, i.e. the rest of the elements of the list
    """

    return tree[1:]

def is_tree(tree):
    """
    Tests if an argument is a tree.

    Args:
        the tree to be tested

    Returns:
        True if its a tree, false otherwise

    """

    if type(tree) != list or len(tree) < 1:
        return False
    for branch in branches(tree):
        if not is_tree(branch):
            return False
    return True

def is_leaf(tree):
    """
    Tests if a tree is a leaf. i.e. it has no branches

    Args:
        tree to be tested

    Returns:
        False if it has branches, True if the list of branhces is empty
    """

    return not branches(tree)




def flatten_list(lst):
    """
    Flattens a nested list.

    Args:
        list to be flattened

    Returns:
        a list that has depth 1

    Raises:
        assertion error if a list isn't passed in
    """
    assert isinstance(lst, list), "you didn't pass a list!"

    if isinstance(lst[0], list):
        if len(lst[0])>1:
            return ['-'.join(i) for i in lst] # then its a kpoints list
        return flatten_list([i[0] for i in lst])
    else:
        return [i for i in lst]



##########################################################################################
##########################################################################################




##########################################################################################
# Classes
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


class Settings:
    """
    A settings object which stores the user settings, only method is to create a tree from these settings.
    """

    def __init__(self, Ecuts, kpts, def_E, def_k):
        """
        Initialises the object with the user settings. This could be changed to instead pass just a path and
            call the necessary functions to pull them out. This current form is more general though...

        Args:
            Ecuts: list of cutoff energies in nested list e.g. [['400'],['450'],['500']]
            kpts: list of kpoints e.g. [['7','1','1'],['8','1','1']]
            def_E: default energy to use when testing kpoints as a list e.g. ['400']
            def_k: default kpoints to use when testing energies e.g. ['7','1','1']
        """

        self.ecuts = Ecuts
        self.kpts = kpts
        self.def_e = def_E
        self.def_k = def_k

        self.e_tree = [] # initialise these as empty until the tree maker methods is called
        self.k_tree = []


    def make_file_tree(self, e_label = 'E_cutoff' , k_label = 'Kpts'):
        """
        Makes tree structures for energies and kpoints being tested.

        Args:
            e_label: the label for the tree for energies
            k_label: the label for the tree for kpoints

        Returns:
            a list of the mutated instance attributes e_tree and k_tree

        Raises:
            Nothing, but tree function call may raise errors
        """

        self.e_tree = tree(e_label, [tree(i) for i in self.ecuts]) # create the tree from the label and the list of
                                                                    # energies in tree format
        self.k_tree = tree(k_label, [tree(i) for i in self.kpts])

        return [self.e_tree, self.k_tree]



##########################################################################################
# for creating directories
##########################################################################################

class FileMaker:
    """
    A class with methods to make directories, add files to them and edit the files
    """

    def __init__(self, settings):
        """
        Initialise the settings to be used, as well as the log files recording what happens.

        Args:
            settings: the settings object which has everything in it
        """

        #self.e_tree = settings.e_tree # Don't do this!!!! change all calls to e_tree to self.settings.e_tree!
        #self.k_tree = settings.k_tree
        self.settings = settings
        self.made_parents = {}
        self.added_files = {}
        self.edit_files = {}


    def dir_maker(self, e_or_k=2):
        """

        """
        trees = [self.settings.e_tree, self.settings.k_tree]

        if e_or_k !=2:
            trees = [trees[e_or_k]]

        for i in trees:
            a = subprocess.Popen(['mkdir', root(i)], stdout = subprocess.PIPE, stderr = subprocess.PIPE) # makes either e.g. kpts or Ecuts
            a = a.communicate()
            b = [subprocess.Popen(["mkdir", root(i)+'/'+"-".join(root(j))], stdout = subprocess.PIPE, stderr = subprocess.PIPE) for j in branches(i)] #makes all the sub directories to be filled. join rejoins the kpoints into one string for this
            b = [i.communicate() for i in b]

            self.made_parents[root(i) + ' made output'] = a
            self.made_parents[root(i) + ' sub directories'] = dict(zip(flatten_list(branches(i)), b))

        return False


    def file_adder(self, e_or_k=2):
        """

        """
        trees = [self.settings.e_tree, self.settings.k_tree]

        if e_or_k !=2:
            trees = [trees[e_or_k]]

        for j in trees:


            a = [subprocess.Popen("cp INPUT/* "+ root(j)+'/'+"-".join(root(i)), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) for i in branches(j)]
            a = [i.communicate() for i in a]
            b = [subprocess.Popen(['rm', root(j) + '/' + "-".join(root(i))+'/conv_params'], stdout = subprocess.PIPE, stderr = subprocess.PIPE) for i in branches(j)] # lazy way to get rid of the accidental conv_params carried over
            b = [i.communicate() for i in b]

            self.added_files[root(j) + ' copied successfully'] = dict(zip(flatten_list(branches(j)), a))
            self.added_files['removed conv_params from all?'] = dict(zip(flatten_list(branches(j)), b))

        return False


    def file_editor(self, e_or_k=2):
        """

        """
        trees = [self.settings.e_tree, self.settings.k_tree]

        if e_or_k !=2:
            trees = [trees[e_or_k]]

        for j in trees:

            param_label = root(j)
            print(param_label)

            self.edit_files[param_label] = 'editing'
            status = []
            if 'Kpts' in param_label:
                for i in branches(j):
                    edit_kpts(param_label, i[0], i) # edit kpts file
                    edit_incar(param_label, self.settings.def_e, i) # edit incar file
                    status.append('Success')
            elif 'E_cutoff' in param_label:
                for i in branches(j):
                    edit_incar(param_label, i[0], i) # edit incar file
                    edit_kpts(param_label, self.settings.def_k, i) # edit kpts file
                    status.append('Success')
            else:
                return "enter a valid label for editing files"

            self.edit_files[param_label+' instances'] = dict(zip(flatten_list(branches(j)), status))
        return False








##########################################################################################
# for analysing
##########################################################################################

class FindFileStructure:
    """

    """

    def __init__(self):
        """

        """

        self.e_sub = []
        self.k_sub = []
        self.def_E = 0
        self.def_k = 0
        self.parent_dirs = {}

    def get_file_structure(self):
        """

        """

        ls = subprocess.check_output(['ls']).decode(
            'ascii')  # get ls and convert it into a normal string (subprocess gives it as repr format)
        ls = (ls.strip()).split(
            '\n')  # ls is newline delimited, split it into elements of a list for each directory or file



        # loop through the ls list and pick out the directory names similar to k or e and what their names are
        self.parent_dirs = self.test_directories(ls, self.kpts_or_ecut_dir)  # two element dictionary with names of parent directories found
        # key is 'e' or 'k' and value is the thing found

        #######################################################
        # figure out a way to do a count on the number of directories with this name to prevent errors
        #######################################################

        # get subdirectories as lists
        e_sub = subprocess.check_output(['ls', self.parent_dirs['e'] + '/']).decode('ascii')
        e_sub = (e_sub.strip()).split('\n')  # split any whitespace off and split into elements

        if 'fly' in e_sub:
            path = self.parent_dirs['e']+'/'
            remove_fly_bodge(path)
            e_sub = subprocess.check_output(['ls', self.parent_dirs['e'] + '/']).decode('ascii')
            e_sub = (e_sub.strip()).split('\n')  # split any whitespace off and split into elements


        k_sub = subprocess.check_output(['ls', self.parent_dirs['k'] + '/']).decode('ascii')
        k_sub = (k_sub.strip()).split('\n')

        if 'fly' in k_sub:
            path = self.parent_dirs['k']+'/'
            remove_fly_bodge(path)
            k_sub = subprocess.check_output(['ls', self.parent_dirs['k'] + '/']).decode('ascii')
            k_sub = (k_sub.strip()).split('\n')



        # check subdirectories are correct, throw away any that aren't the right format i.e. energies must be a integer
        # kpts must be a string containing two dashes
        e_sub = list(self.test_directories(e_sub, self.e_test).keys())  # dictionary where keys and values are equal. contains only directories that passed the test
        k_sub = list(self.test_directories(k_sub, self.k_test).keys())  # we pull out just the keys because we want the keys

        self.e_sub = e_sub
        self.k_sub = k_sub
        self.def_k = self.get_def(self.parent_dirs['e'], self.e_sub[0], 'KPOINTS', 'Gamma')
        self.def_E = self.get_def(self.parent_dirs['k'], self.k_sub[0], 'INCAR', 'ENCUT')

        return False


    def get_def(self, parent, dir, file, line_key):
        """

        """

        file = parent + '/' + dir + '/' + file
        lines = open(file).read().splitlines()  # read lines into elements in a list

        line_no = False
        for i in lines:
            if line_key in i:
                line_no = lines.index(i)
        if not line_no:
            return 'Line key not found in files. Ensure that INCAR has ENCUT in it, and KPOINTS has Gamma in it, just before where KPOINTS should be added'
        if line_key == 'Gamma':
            line_no += 1
        line_with_info = lines[line_no]
        default = self.strip_line(line_with_info, line_key)

        return default

    def strip_line(self, line, line_key):
        """

        """

        if line_key == 'Gamma':
            value = '-'.join((line.strip().split()))
        else:
            value = re.sub("[^0-9]", "", line) # this strips all non-numeric characters from the line

        return value

    def kpts_or_ecut_dir(self, directory):
        """

        """

        kpts = ['Kpts', 'kpts', 'KPOINTS', 'kpoints', 'KPTS', 'k', 'K', 'Ks', 'ks']
        ecuts = ['E_cutoff', 'E', 'Ecuts', 'Es', 'energy', 'Energy', 'ecuts']
        if directory in kpts:
            return 'k'
        elif directory in ecuts:
            return 'e'
        else:
            return False

    def test_directories(self, ls, test_fn):
        """

        """

        dir_names = {}
        for i in ls:
            test_for_dir = test_fn(i)
            if test_for_dir:  # i.e. if it satisfies the test
                dir_names[test_for_dir] = i

        return dir_names

    def e_test(self, directory):
        """

        """

        if type(eval(directory)) == int:
            return directory
        return False

    def k_test(self, directory):
        """

        """

        if type(directory) == str and directory.count('-') == 2:
            return directory
        return False



class FindEnergies:
    """

    """

    def __init__(self, e_tree, k_tree):
        """

        """

        self.e_tree = e_tree
        self.k_tree = k_tree
        self.energies_list = {}


    def parent_iterator(self):
        """

        """

        e_name = root(self.e_tree)
        self.file_iterator(e_name, self.e_tree)

        k_name = root(self.k_tree)
        self.file_iterator(k_name, self.k_tree)

        return False

    def file_iterator(self, parent_name, tree):
        """

        """

        self.energies_list[parent_name] = {}
        for b in branches(tree): #iterate over the directories listed in the tree
            file = parent_name + '/' + root(b) + '/' + 'OUTCAR'  # path to OUTCAR file
            lines = open(file).read().splitlines()  # read lines into elements in a list
            # iterate over the lines
            for i in lines:
                if 'y=' in i:
                    numbers_on_line = re.findall('\d+', i)
                    energy = -float(str(numbers_on_line[0]) + '.' + str(numbers_on_line[1])) # the first element of numbers is before decimal, second is after

                    self.energies_list[parent_name][root(b)] = energy
                    # add this energy to dictionary with the current dir as label

            if root(b) not in self.energies_list[parent_name]:
                self.energies_list[parent_name][root(b)] = "'y=' wasn't found in the OUTCAR file"
        return False

class WriteCSVFile:
    """

    """

    def __init__(self, energies_object, foldername, filename, settings):
        """

        """

        self.energies = energies_object.energies_list
        self.filename = filename+'.csv'
        self.foldername = foldername
        self.settings = settings
        self.e_name = root(energies_object.e_tree)
        self.k_name = root(energies_object.k_tree)
        self.accessed = datetime.datetime.now().strftime("%Y-%m-%d")
        self.cwd = os.getcwd()  # of the script!
        self.time = datetime.datetime.now().time().strftime("%H:%M:%S")

    def write(self):
        """

        """

        path = self.cwd + '/' + self.foldername + '/' + self.filename
        e_def = self.settings.def_e
        k_def = self.settings.def_k


        #write information at top
        lines = ["", "This file shows the convergence details for files checked in this folder.", "file was written (data recorded) on: "+ self.accessed + " at: " + self.time, ""]  # make lines
        open(path, 'a').write('\n'.join(lines))  # write list to file


        #title for ecuts, default kpts used
        lines = ["convergence details for cutoff energy testing follows", "note the default kpts used were: " + str(k_def), ""]  # make lines
        open(path, 'a').write('\n'.join(lines))  # write list to file



        #write ecuts to csv
        e_tuples = [(i, j) for i, j in self.energies[self.e_name].items()] # convert the dictionary into a list of tuples for sorting
        e_tuples = sorted(e_tuples, key = lambda x: int(x[0])) # this sorts the tuples - they came from a dict
        for i in range(len(e_tuples)-1):
            e_tuples[i+1] = (e_tuples[i+1][0], e_tuples[i+1][1], e_tuples[i+1][1] - e_tuples[i][1])
        with open(path, 'a') as f:
            w = csv.writer(f)
            w.writerow(['cutoff energy', ' energy (eV)', ' difference (eV)'])
            for row in e_tuples:
                w.writerow(row)


        # title for kpts, default ecuts used
        lines = ["","convergence details for kpts testing follows", "note the default energy used was: " + e_def +" ev",
                 ""]  # make lines
        open(path, 'a').write('\n'.join(lines))  # write list to file


        # write kpts to csv
        k_tuples = [(i, j) for i, j in self.energies[self.k_name].items()]
        k_tuples = sorted(k_tuples, key=lambda x: (int(x[0][0]), int(x[0][2]), int(x[0][4:]))) # sort them by the kpoints
        for i in range(len(k_tuples)-1):
            k_tuples[i+1] = (k_tuples[i+1][0], k_tuples[i+1][1], k_tuples[i+1][1] - k_tuples[i][1])
        with open(path, 'a') as f:
            w = csv.writer(f)
            w.writerow(['kpoints', 'energy (eV)'])
            for row in k_tuples:
                w.writerow(row)



        return False

    def add_POSCAR(self):
        """

        """

        a = [subprocess.Popen("cp INPUT/POSCAR " + str(self.foldername), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)]
        a = [i.communicate() for i in a]



        return False


##########################################################################################
# log file
##########################################################################################
class log_file:
    """

    """

    def __init__(self, settings, filemaker, title = 'log_file'):
        """

        """

        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.time = datetime.datetime.now().time().strftime("%H:%M:%S")
        self.cwd = os.getcwd() # of the script!
        self.settings = settings
        self.def_e = settings.def_e
        self.def_k = settings.def_k
        self.trees = [settings.e_tree, settings.k_tree]
        self.made_dirs = filemaker.made_parents
        self.added_files = filemaker.added_files
        self.edited_files = filemaker.edit_files
        self.title = title


    def create_log(self):
        """

        """

        path = self.cwd + '/' + self.title
        lines = ["Date: " + self.date, "Time: " + self.time, "Working Directory: " + self.cwd, "default E: " + str(self.def_e), "default Kpts: " + str(self.def_k) ,"","file trees: ", str(self.trees),
                 "","NOTE: the format of these logs is not great, if there are no obvious messages then everything is okay","","logs from making directories: ", str(self.made_dirs), "","logs from adding files: ", str(self.added_files), "",
                 "logs from editing files: ", str(self.edited_files)]  # make lines
        open(path, 'w').write('\n'.join(lines))  # write list to file


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



##########################################################################################
# functions to run
##########################################################################################
def make_conv_files():
    """

    """
    path = ''
    path += 'INPUT/conv_params'
    choices = get_params(path) # an object which has all the settings stored in it
    file_trees = choices.make_file_tree() # first element is the Ecut tree, second is kpts tree
    files = FileMaker(choices) # make files object with trees to be created and choices

    files.dir_maker() # make directories

    files.file_adder() # add stock files to directories

    files.file_editor() # edit the files as needed

    log = log_file(choices, files) # create log file
    log.create_log()



def analyse_conv_files():
    """

    """

    a = subprocess.Popen(['mkdir', 'convergence_details'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # makes a folder to put the convergence details in
    files = FindFileStructure()
    files.get_file_structure() # this sets objects values to the cutoff energies, kpts tested and the defaults that were used.
    # note that the defaults are only tested on the first listed directory that was deemed to be correct. e.g. 1000
    choices = Settings(files.e_sub, files.k_sub, files.def_E, files.def_k) # create settings object with the run conv testing settings

    file_trees = choices.make_file_tree(files.parent_dirs['e'], files.parent_dirs['k']) # list with first element as e_tree and second element as k_tree

    find_energies = FindEnergies(choices.e_tree, choices.k_tree)
    find_energies.parent_iterator()

    energy_dict = find_energies.energies_list # this returns a dictionary with two keys, E and K (or their names as directories found already)

    #scale the energy by the number of atoms in the cell
    no_atoms = get_no_atoms()
    for i in energy_dict:
        for j in energy_dict[i]:
            energy_dict[i][j] = energy_dict[i][j]/no_atoms

    # the values of these are dictionaries with keys as the names of the folders and values as the energy in the OUTCAR

    # you should edit that class to also save other important details from runs

    written = WriteCSVFile(find_energies, 'convergence_details', 'convergence data', choices)
    written.write() # writes all the convergence details analysed into a csv file.
    # this has details and list of energies and energy differences by sorted values
    written.add_POSCAR() # copies the poscar from input into this folder


#######################The file directories to create on the fly have to be named either E_cutoff or Kpts. Fix this.
def on_the_fly():
    """

    """

    create_or_edit = input("please enter either 'c' or 'e' to choose whether to create new directories on the fly, "
                           "or to edit all instances of a file in subdirectories of current directory. \n")
    if create_or_edit == 'e':
        print("Have only implemented create right now, sorry.")
        return False

    # create instructions
    k_or_e = input("Please enter either 'k' or 'e' to choose which convergence directories you would like to make: "
                   "a set of kpoints directories or a set of energy directories.\n")

    # read in the kpoints or energies to create
    if k_or_e =='k':
        kpts = input("please enter the kpoints you would like to create, in the format x-y-z,i-j-k, etc. i.e. comma "
                     "seperated the directories to create, and dash seperate the numbers. \n")
        def_e = [input("please enter the default energy you would like to run at, leave blank to get this from"
                      "the input folder \n")]

        kpts =  [i.split('-') for i in kpts.split(',')] # split kpoints into a list of lists of elements

        choices = Settings([], kpts, def_e, '0-0-0') # pass in the defaults for what you want to make
        k_tree = choices.make_file_tree()[1] # the ktree is the 1st element in the list returned
        k_tree[0] += '/fly' # put into a subdirectory fly for ease of running
        files = FileMaker(choices)#, tree(''), k_tree) # make the object for creating everything

        files.dir_maker(1) #k_tree) # make the directories

        files.file_adder(1) #files.k_tree) # add the default files from INPUT

        files.file_editor(1) #k_tree) # edit the files according to the user inputs

        date = datetime.datetime.now().strftime("%Y_%m_%d")
        time = datetime.datetime.now().time().strftime("%H_%M_%S")
        log = log_file(choices, files, 'on_the_fly' + date + "_" + time)
        log.create_log() # create log file with date and time


    # if they want to make some energies do this
    elif k_or_e =='e':
        energies = input("please enter the energies you would like to create, in the format xyz,ijk, etc. i.e. comma "
                     "seperated directories to create\n")
        def_k = input("please enter the default kpoints you would like to run at, format: x-y-z"
                      "leave blank to get this from the alphabetically first subdirectory here\n")

        energies = energies.split(',') # split the input on commas
        energies = [[i] for i in energies] # make sure the energies are themselves lists for objects to work properly
        def_k = def_k.split('-') #split the kpoints into a list of elements

        choices = Settings(energies, [], '0', def_k) # pass in the defaults for what you want to make
        e_tree = choices.make_file_tree()[0] # the etree is the 0th element in the list returned
        e_tree[0] += '/fly' # put into a subdirectory fly for ease of running
        files = FileMaker(choices) #, e_tree, tree('')) # make the object for creating everything

        files.dir_maker(0) #e_tree) # make the directories

        files.file_adder(0) #files.e_tree) # add the default files from INPUT

        files.file_editor(0) #e_tree) # edit the files according to the user inputs

        date = datetime.datetime.now().strftime("%Y_%m_%d")
        time = datetime.datetime.now().time().strftime("%H_%M_%S")
        log = log_file(choices, files, 'on_the_fly' + date + "_" + time)  # create log file
        log.create_log() # create log file with date and time

    # if neither k or e was entered
    else:
        print('You entered the wrong input!!')
        return False

    return False

##########################################################################################
##########################################################################################





##########################################################################################
# arg parser
##########################################################################################
parser = argparse.ArgumentParser(description = 'Make or analyse convergence docs for vasp.')
parser.add_argument('-f', '--function', help = 'Enter m for making files, a for analysing files', required = True)
args = vars(parser.parse_args())

print("please note this program is not fully tested and may contain bugs. I advise reading the readme and checking "
      "all files before running \n")

if args['function'] == 'm':
    make_conv_files()
elif args['function'] == 'a':
    analyse_conv_files()
elif args['function'] == 'fly':
    print('On the fly functionality not fully implemented yet! \n')
    on_the_fly()
else:
    print('enter either m or a to run the program. m will make files, a will analyse them')

##########################################################################################
#########################################################################################
