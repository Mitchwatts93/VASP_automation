import subprocess
import os
import datetime

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

def get_params():

    # cd to Input, read the conv_params file in and pass each line to file reader
    list = file_reader('INPUT/conv_params')
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


    def make_file_tree(self):

        e_tree = tree('E_cutoff', [tree(i) for i in self.ecuts])
        k_tree = tree('Kpts', [tree(i) for i in self.kpts])

        return [e_tree, k_tree]


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
        path = self.cwd + '/' + 'logfile'
        lines = ["Date: " + self.date, "Time: " + self.time, "Working Directory: " + self.cwd, "default E: " + str(self.def_e), "default Kpts: " + str(self.def_k) ,"","file trees: ", str(self.trees),
                 "","NOTE: the format of these logs is not great, if there are no obvious messages then everything is okay","","logs from making directories: ", str(self.made_dirs), "","logs from adding files: ", str(self.added_files), "",
                 "logs from editing files: ", str(self.edited_files)]  # make lines
        open(path, 'w').write('\n'.join(lines))  # write list to file



##########################################################################################
##########################################################################################

##########################################################################################
# main run
##########################################################################################


choices = get_params() # an object which has all the settings stored in it
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




# make short scripts for adding extra directories
# make program for pulling out key info into txt doc (maybe even excel?)
# make program for organising data for scp in format organised
