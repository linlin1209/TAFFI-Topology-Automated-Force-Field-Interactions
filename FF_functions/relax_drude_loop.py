#!/bin/env python                                                                                                                                                              

import sys,os,argparse

def main(argv):
    
    # This is not written with relax_drude.py since it's easier to do single relax test with independent ,py file

    parser = argparse.ArgumentParser(description='takes in the dimer xyz file in configs folder, do drude relaxation, this program assumes there is a FF folder that contains lin.prm and lin.rtf in vdw folder')

    parser.add_argument('-tag',dest='tag',default='',
                        help = 'tag for FF folder and relax drude, this assumes FF name is FF_tag and output drude: _drude_tag.xyz, default:None')

    
    args=parser.parse_args(argv)    

    fun(args.tag)

# make a function so it's easier for drude_fit to call
def fun(tag):

    current_dir = os.getcwd()

    folders = [ i for i in os.listdir(current_dir) if (os.path.isdir(i) and i != 'FF') ]

    for count_i,i in enumerate(folders):
         print('Working on {}, {} configurations left'.format(i,len(folders)-count_i),end="\r")
         os.chdir('{}/{}'.format(current_dir,i))
         if tag == '':
            relax_drude.main('-xyz {}.xyz --silent '.format(i).split())
         else:
            relax_drude.main('-xyz {}.xyz --silent -tag {}'.format(i,tag).split())
    return

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder, "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
