import subprocess

# open conv_params in INPUT, save them

# Edit the input files accordingly

# recursively create each directory, copy all the files from INPUT, edit them accordingly

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

def file_editer(filename, default):

    return True


def dir_maker(tree):





    return True


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
##########################################################################################
##########################################################################################

##########################################################################################
# funs and classes
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

def get_params():

    # cd to Input, read the conv_params file in and pass each line to file reader
    list = file_reader('INPUT/conv_params')
    # return an object choices to be passed to settings
    Ecuts =  list[0]
    kpts = list[1]
    def_E = list[2]
    def_k = list[3]
    params = Settings(Ecuts, kpts, def_E, def_k)

    return params


class FileMaker:

    def __init__(self, settings, e_tree, k_tree):

        self.e_tree = e_tree
        self.k_tree = k_tree
        self.settings = settings

    def dir_maker(self, tree):

        subprocess.call(['mkdir', root(tree)]) # makes either e.g. kpts or Ecuts
        a = [subprocess.call(['mkdir', root(tree)+'/'+"-".join(root(i))]) for i in branches(tree)] #makes all the sub directories to be filled. join rejoins the kpoints into one string for this
        return a

#next bit to do! edit each file in each directory as needed
    def file_editor(self, blah):
        e_def = self.settings.def_e  # default E and K to be used
        k_def = self.settings.def_k


        return file

    def file_adder(self, tree):

        a = [subprocess.call('cp INPUT/* '+root(tree)+'/'+"-".join(root(i)), shell=True) for i in branches(tree)]
        [subprocess.call('rm ' + root(tree) + '/' + "-".join(root(i))+'/conv_params', shell=True) for i in branches(tree)] # lazy way to get rid of the accidental conv_params carried over

        return a







##########################################################################################
##########################################################################################

##########################################################################################
# main run
##########################################################################################


choices = get_params() # an object which has all the settings stored in it
file_trees = choices.make_file_tree() # first element is the Ecut tree, second is kpts tree
files = FileMaker(choices, file_trees[0], file_trees[1])
rec = {} #record of making files etc
rec['E_dirs'] = files.dir_maker(files.e_tree) #make directories and add successful notes to rec
rec['K_dirs'] = files.dir_maker(files.k_tree)

rec['E_files'] = files.file_adder(files.e_tree)
rec['K_files'] = files.file_adder(files.k_tree)


# make short scripts for adding extra directories
# make program for pulling out key info into txt doc (maybe even excel?)
# make program for organising data for scp in format organised
