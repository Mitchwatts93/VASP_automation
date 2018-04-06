import subprocess

# open conv_params in INPUT, save them

# Edit the input files accordingly

# recursively create each directory, copy all the files from INPUT, edit them accordingly

##########################################################################################
# abstraction functions
##########################################################################################

def file_reader(filename):

    info = []

    return info

def file_editer(filename, default):

    return True


def dir_maker(list):
    subprocess.call(['mkdir', list[0]])

    return True
##########################################################################################
##########################################################################################

##########################################################################################
# funs and classes
##########################################################################################
class settings:

    def __init__(self, Ecuts, kpts, def_E, def_k):
        self.ecuts = Ecuts
        self.kpts = kpts
        self.def_e = def_E
        self.def_k = def_k


    def tree(list):
        return tree


def get_params():

    # cd to Input, read the conv_params file in and pass each line to file reader

    # return a list of lists to be passed to settings

    return list


##########################################################################################
##########################################################################################

##########################################################################################
# main run
##########################################################################################

def run():
     list = get_params()
     choices = settings(list[0], list[1], list[2], list[3])


# make short scripts for adding extra directories
# make program for pulling out key info into txt doc (maybe even excel?)
# make program for organising data for scp in format organised
