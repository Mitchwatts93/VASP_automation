import subprocess
import os
import datetime
import re
import csv
import argparse

# note on subprocess. This may not be implemented well. Have not fully read documentation.
# using shell = True can be a security hazard due to shell injection esp if command string is constructed from external
# input (not really an issue here but not used for good practice)

# further note: uisng pipe to capture output may cause issues if pipe buffer is filled, in which case it will cause the program to hang.
# remove stderr/stdout from subprocess calls if this happens - this will cause bash outputs to be seen during running and will not allow details to be saved to log

# documentation suggests using Popen with communicate() method when pipes are needed -> this seems to be after already using pipe
# e.g. a.communicate() is the same as a.stderr.read(). not sure if this will cause an issue though...


##########################################################################################
# to do
##########################################################################################
# - add docstrings and doctests
# - add more info to the output
# add functionality for creating kpoints or Ecuts on the fly
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
# abstraction functions
##########################################################################################

def file_reader(filename):
    with open(filename) as f:
        info = f.readlines()
    info = [i.strip() for i in info] # each element in info is a string of a line from the file
    info = [i.split() for i in info] #split each whitespace delimited element into a list of lists
    info = [[i.split('-') for i in j] for j in info] # note info is 3 layers deep

    #info[0] = [i[0] for i in info[0]] #removes some of the nested structure. now its just a
    info[2] = info[2][0] #makes default E just a number
    info[3] = info[3][0]

    return info

def edit_kpts(param_label, i, dir):
    line_key = 'Gamma'
    #line_no = 3
    file = 'KPOINTS'
    replacement_line = "  " + i[0] + "  " + i[1] + "  " + i[2]
    gen_file_editor(param_label, dir, file, replacement_line, line_key)
    return False

def edit_incar(param_label, i, dir):
    line_key = 'ENCUT'
    #line_no = 18
    file = 'INCAR'
    replacement_line = "  ENCUT   = " + i[0] + "       ! Plane-wave cutoff"
    gen_file_editor(param_label, dir, file, replacement_line, line_key)
    return False

def gen_file_editor(param_label, i, file, replacement_line, line_key):

    file = param_label + '/' + "-".join(root(i)) + '/' + file  # path to file
    lines = open(file).read().splitlines()  # read lines into elements in a list

    line_no = False

    for i in lines:
        if line_key in i:
            line_no = lines.index(i)


    if not line_no:
        return 'Line key not found in files. Ensure that INCAR has ENCUT in it, and KPOINTS has Gamma in it, just before where KPOINTS should be added'

    if line_key == 'Gamma':
        line_no += 1

    lines[line_no] = replacement_line  # replace the kpoints details in the list
    open(file, 'w').write('\n'.join(lines))  # write list to file

    return False

def get_params(path):

    # cd to Input, read the conv_params file in and pass each line to file reader
    list = file_reader(path)
    # return an object choices to be passed to settings
    Ecuts =  list[0]

    start = int(Ecuts[0][0])
    multiplier = int(Ecuts[1][0])
    end = int(Ecuts[2][0])
    E_range = (end - start)//multiplier +1

    Es = [i*multiplier for i in range(E_range)]
    Ecuts = [[str(i+start)] for i in Es]
    kpts = list[1]
    def_E = list[2]
    def_k = list[3]
    params = Settings(Ecuts, kpts, def_E, def_k)

    return params





##########################################################################################
#tree stuff
##########################################################################################
def tree(root, branches = []):
    """"give it elements of a two element  list. first is a string, either 'Ecuts' or 'kpts'. second is the list of elements to make"""
    for branch in branches:
        assert is_tree(branch), 'branches must be trees'
    return [root] + list(branches)

def root(tree):
    return tree[0]

def branches(tree):
    return tree[1:]

def is_tree(tree):
    if type(tree) != list or len(tree) < 1:
        return False
    for branch in branches(tree):
        if not is_tree(branch):
            return False
    return True

def is_leaf(tree):
    return not branches(tree)

def flatten_list(lst):

    if isinstance(lst[0], list):
        if len(lst[0])>1:
            return ['-'.join(i) for i in lst] # then its a kpoints list
        return flatten_list([i[0] for i in lst])
    else:
        return [i for i in lst]



##########################################################################################
##########################################################################################




##########################################################################################
##########################################################################################

##########################################################################################
# Classes
##########################################################################################
class Settings:

    def __init__(self, Ecuts, kpts, def_E, def_k):

        self.ecuts = Ecuts
        self.kpts = kpts
        self.def_e = def_E
        self.def_k = def_k

        self.e_tree = []
        self.k_tree = []


    def make_file_tree(self,e_label = 'E_cutoff' , k_label = 'Kpts'):

        self.e_tree = tree(e_label, [tree(i) for i in self.ecuts])
        self.k_tree = tree(k_label, [tree(i) for i in self.kpts])

        return [self.e_tree, self.k_tree]


class FileMaker:

    def __init__(self, settings, e_tree, k_tree):

        self.e_tree = e_tree
        self.k_tree = k_tree
        self.settings = settings
        self.made_parents = {}
        self.added_files = {}
        self.edit_files = {}


    def dir_maker(self, tree):

        a = subprocess.Popen(['mkdir', root(tree)], stdout = subprocess.PIPE, stderr = subprocess.PIPE) # makes either e.g. kpts or Ecuts
        a = a.communicate()
        b = [subprocess.Popen(["mkdir", root(tree)+'/'+"-".join(root(i))], stdout = subprocess.PIPE, stderr = subprocess.PIPE) for i in branches(tree)] #makes all the sub directories to be filled. join rejoins the kpoints into one string for this
        b = [i.communicate() for i in b]

        self.made_parents[root(tree) + ' made output'] = a
        self.made_parents[root(tree) + ' sub directories'] = dict(zip(flatten_list(branches(tree)), b))


        return False


    def file_editor(self, tree):
        e_def = self.settings.def_e  # default E and K to be used)
        k_def = self.settings.def_k
        param_label = root(tree)

        self.edit_files[param_label] = 'editing'
        status = []
        if param_label == 'Kpts':
            for i in branches(tree):
                edit_kpts(param_label, i[0], i) # edit kpts file
                edit_incar(param_label, e_def, i) # edit incar file
                status.append('Success')
        elif param_label =='E_cutoff':
            for i in branches(tree):
                edit_incar(param_label, i[0], i) # edit incar file
                edit_kpts(param_label, k_def, i) # edit kpts file
                status.append('Success')
        else:
            return "enter a valid label for editing files"

        self.edit_files[param_label+' instances'] = dict(zip(flatten_list(branches(tree)), status))
        return False


    def file_adder(self, tree):

        a = [subprocess.Popen("cp INPUT/* "+ root(tree)+'/'+"-".join(root(i)), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) for i in branches(tree)]
        a = [i.communicate() for i in a]
        b = [subprocess.Popen(['rm', root(tree) + '/' + "-".join(root(i))+'/conv_params'], stdout = subprocess.PIPE, stderr = subprocess.PIPE) for i in branches(tree)] # lazy way to get rid of the accidental conv_params carried over
        b = [i.communicate() for i in b]

        self.added_files[root(tree) + ' copied successfully'] = dict(zip(flatten_list(branches(tree)), a))
        self.added_files['removed conv_params from all?'] = dict(zip(flatten_list(branches(tree)), b))

        return False



class log_file:

    def __init__(self, settings, filemaker):
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.time = datetime.datetime.now().time().strftime("%H:%M:%S")
        self.cwd = os.getcwd() # of the script!
        self.settings = settings
        self.def_e = settings.def_e
        self.def_k = settings.def_k
        self.trees = [filemaker.e_tree, filemaker.k_tree]
        self.made_dirs = filemaker.made_parents
        self.added_files = filemaker.added_files
        self.edited_files = filemaker.edit_files


    def create_log(self):
        path = self.cwd + '/' + 'log_file'
        lines = ["Date: " + self.date, "Time: " + self.time, "Working Directory: " + self.cwd, "default E: " + str(self.def_e), "default Kpts: " + str(self.def_k) ,"","file trees: ", str(self.trees),
                 "","NOTE: the format of these logs is not great, if there are no obvious messages then everything is okay","","logs from making directories: ", str(self.made_dirs), "","logs from adding files: ", str(self.added_files), "",
                 "logs from editing files: ", str(self.edited_files)]  # make lines
        open(path, 'w').write('\n'.join(lines))  # write list to file



class FindFileStructure:

    def __init__(self):
        self.e_sub = []
        self.k_sub = []
        self.def_E = 0
        self.def_k = 0
        self.parent_dirs = {}

    def get_file_structure(self):

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

        if line_key == 'Gamma':
            value = '-'.join((line.strip().split()))
        else:
            value = re.sub("[^0-9]", "", line) # this strips all non-numeric characters from the line

        return value

    def kpts_or_ecut_dir(self, directory):

        kpts = ['Kpts', 'kpts', 'KPOINTS', 'kpoints', 'KPTS', 'k', 'K', 'Ks', 'ks']
        ecuts = ['E_cutoff', 'E', 'Ecuts', 'Es', 'energy', 'Energy', 'ecuts']
        if directory in kpts:
            return 'k'
        elif directory in ecuts:
            return 'e'
        else:
            return False

    def test_directories(self, ls, test_fn):
        dir_names = {}
        for i in ls:
            test_for_dir = test_fn(i)
            if test_for_dir:  # i.e. if it satisfies the test
                dir_names[test_for_dir] = i

        return dir_names

    def e_test(self, directory):
        if type(eval(directory)) == int:
            return directory
        return False

    def k_test(self, directory):
        if type(directory) == str and directory.count('-') == 2:
            return directory
        return False


class FindEnergies:

    def __init__(self, e_tree, k_tree):
        self.e_tree = e_tree
        self.k_tree = k_tree
        self.energies_list = {}


    def parent_iterator(self):

        e_name = root(self.e_tree)
        self.file_iterator(e_name, self.e_tree)

        k_name = root(self.k_tree)
        self.file_iterator(k_name, self.k_tree)

        return False

    def file_iterator(self, parent_name, tree):

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

    def __init__(self, energies_object, foldername, filename, settings):
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
            w.writerow(['cutoff energy', 'energy (eV)', 'difference (eV)'])
            for row in e_tuples:
                w.writerow(row)


        # title for kpts, default ecuts used
        lines = ["","convergence details for kpts testing follows", "note the default energy used was: " + e_def +" ev",
                 ""]  # make lines
        open(path, 'a').write('\n'.join(lines))  # write list to file


        # write kpts to csv
        k_tuples = [(i, j) for i,j in self.energies[self.k_name].items()]
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
        a = [subprocess.Popen("cp INPUT/POSCAR " + str(self.foldername), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)]
        a = [i.communicate() for i in a]



        return False

##########################################################################################
##########################################################################################

##########################################################################################
# functions to run
##########################################################################################
def make_conv_files():
    path = 'INPUT/conv_params'
    choices = get_params(path) # an object which has all the settings stored in it
    file_trees = choices.make_file_tree() # first element is the Ecut tree, second is kpts tree
    files = FileMaker(choices, file_trees[0], file_trees[1]) # make files object with trees to be created and choices

    files.dir_maker(files.e_tree) # make directories
    files.dir_maker(files.k_tree)

    files.file_adder(files.e_tree) # add stock files to directories
    files.file_adder(files.k_tree)

    files.file_editor(files.e_tree) # edit the files as needed
    files.file_editor(files.k_tree)

    log = log_file(choices, files) # create log file
    log.create_log()



def analyse_conv_files():
    a = subprocess.Popen(['mkdir', 'convergence_details'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # makes a folder to put the convergence details in
    files = FindFileStructure()
    files.get_file_structure() # this sets objects values to the cutoff energies, kpts tested and the defaults that were used.
    # note that the defaqults are only tested on the first listed directory that was deemed to be correct. e.g. 1000
    choices = Settings(files.e_sub, files.k_sub, files.def_E, files.def_k) # create settings object with the run conv testing settings

    file_trees = choices.make_file_tree(files.parent_dirs['e'], files.parent_dirs['k']) # list with first element as e_tree and second element as k_tree

    find_energies = FindEnergies(choices.e_tree, choices.k_tree)
    find_energies.parent_iterator()

    energy_dict = find_energies.energies_list # this returns a dictionary with two keys, E and K (or their names as directories found already)
    # the values of these are dictionaries with keys as the names of the folders and values as the energy in the OUTCAR

    # you should edit that class to also save other important details from runs

    written = WriteCSVFile(find_energies, 'convergence_details', 'convergence data', choices)
    written.write() # writes all the convergence details analysed into a csv file.
    # this has details and list of energies and energy differences by sorted values
    written.add_POSCAR() # copies the poscar from input into this folder

##########################################################################################
##########################################################################################


##########################################################################################
# arg parser
##########################################################################################
parser = argparse.ArgumentParser(description = 'Make or analyse convergence docs for vasp.')
parser.add_argument('-f', '--function', help = 'Enter m for making files, a for analysing files', required = True)
args = vars(parser.parse_args())

if args['function'] == 'm':
    make_conv_files()
elif args['function'] == 'a':
    analyse_conv_files()
else:
    print('enter either m or a to run the program. m will make files, a will analyse them')

##########################################################################################
##########################################################################################
