#!/bin/env python                                                                                                                                                              

# Description: This script combines msd data from several trajectories and plots the t vs msd curve
import sys,argparse,os,ast,re,fnmatch,matplotlib,subprocess,scandir
from numpy import *
from pylab import *

def main(argv):

    parser = argparse.ArgumentParser(description='Average over all files matching the -f string.')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-f', dest='Filename', default='*.in',
                        help = 'The program operates on all input files discovered during a directory walk from the working directory whose name matches this variable. '\
                               'For example, if the user supplies *.in, then all files ending in .in within all subfolders relative the working directory will be submitted. '\
                               'Note that the program uses python''s regular expression matching for finding files, so ''*'' can be used as a wildcard character '\
                               'or as a literal by using the ''\\'' as an escape character. When basic wildcards are used, the program will submit all input files matching '\
                               'the wildcard pattern as separate jobs.') 

    parser.add_argument('-p', dest='procs', default=1,
                        help = 'Specifies the number of processors for each job (default: 1; max 8)')
                        
    parser.add_argument('-o', dest='outputname', default="orca_mass",
                        help = 'Specifies the job name (default: orca_mass)')

    parser.add_argument('-t', dest='walltime', default=48,
                        help = 'Specifies the walltime for each job (default: 48, hours by default, if Xmin is used then the argument will be interpretted in minutes)')

    parser.add_argument('-q', dest='queue', default='ccm_queue',
                        help = 'Specifies the queue for the job (default: ccm_queue; see NERSC webpage for valid options)')

    parser.add_argument('-ppn', dest='ppn', default=24,
                        help = 'Specifies the number of processors per node on the cluster architecture. the -ppn %% -p should equal zero. (default: 24)') 

    parser.add_argument('-size', dest='size', default=100,
                        help = 'Specifies the number of calculations to bundle per job (default: 100)')

    parser.add_argument('--overwrite', dest='overwrite', default=0, action='store_const', const=1,
                        help = 'When set, if previous run data is discovered in the submission folders it is deleted. (default: off)')

    parser.add_argument('-path_to_exe',dest='path_to_exe', default="/global/homes/b/bsavoie/Orca/3_0_2/orca",
                        help = 'specifies the shell submission script to call.')

    args=parser.parse_args()
    if type(args.walltime) == str and "min" in args.walltime: args.walltime = int(args.walltime.split('min')[0]); min_flag = 1
    else: args.walltime = int(args.walltime); min_flag = 0
    args.procs = int(args.procs)
    args.ppn = int(args.ppn)
    args.size = int(args.size)
    Filename = args.Filename
    working_dir = os.getcwd()

    # Check that the number of processors per job divides into the number of processors per node.
    if args.ppn % args.procs != 0:
        print "ERROR: the -ppn % -p must be zero to ensure that jobs aren't split across nodes. Exiting..."
        quit()

    # Find the input files:
    # Use fnmatch to use wildcards and os.walk to perform a directory walk. Names holds the path and filename of each 
    # file to be processed by the program
    Names = [ _ for _ in set([ f for dp, dn, filenames in scandir.walk('.') for f in filenames if fnmatch.fnmatch(f,Filename) ]) ]

    # Create a dictionary from the filenames, where each dictionary key corresponds to a filename and each entry is a list
    # of subfiles to be processed as a batch. e.g., molecule.in might show up in twenty subfolders. molecule.in would end
    # up as a key, and the list of individual instances of molecule.in would constitute the entry.
    Files = {}
    Files[Filename] = [ os.path.join(dp, f) for dp, dn, filenames in scandir.walk('.') for f in filenames if fnmatch.fnmatch(f,Filename) ]

    #################################################################
    # Iterate over all discovered input files, cd to the containing #
    # directory and check if the job already has an output file     #
    #################################################################    
    input_files = []
    input_paths = []
    unsubmitted_files = []
    unsubmitted_paths = []
    incomplete_files = []
    incomplete_paths = []
    complete_files = []
    complete_paths = []

    for i in Files.keys():

        # Sort the filenames as string
        Current_files = sorted(Files[i])

        print "\nTotal number of jobs discovered: {}".format(len(Current_files))
#         print "{}".format("#"*80)
#         print "# {:^76s} #".format("PROCESSING THE FOLLOWING FILES")
#         print "{}".format("#"*80)
#         for j in Current_files:
#             print j

        for j in Current_files:

            # Change to the directory holding the current file
            path_to_file = '/'.join(j.split('/')[0:-1])
            if path_to_file != '':
                os.chdir(path_to_file)

            # If no output file is discovered in the run folder then add the file to the unsubmitted lists
            current_name = j.split('/')[-1]
            if current_name.split('.')[0]+".out" not in os.listdir('.'):
                unsubmitted_files += [current_name]
                unsubmitted_paths += [path_to_file]

            # If an output is present, a check for completeness is performed. 
            else:
                completed_flag = 0
                with open(current_name.split('.')[0]+".out",'r') as f:
                    for lines in f:
                        if "****ORCA TERMINATED NORMALLY****" in lines:
                            completed_flag = 1

                # Add incompleted files to the incomplete lists
                if completed_flag == 0:
                    incomplete_files += [current_name]
                    incomplete_paths += [path_to_file]

                # Add completed files to the complete lists
                else:
                    complete_files += [current_name]
                    complete_paths += [path_to_file]

#             elif args.overwrite == 1:
#                 files = [ k for k in os.listdir('.') if k != j.split('/')[-1] ]
#                 for k in files:
#                     os.remove(k)
#                 input_files += [current_name]
#                 input_paths += [path_to_file]
                                    
            os.chdir(working_dir)

    print "\nNumber of Unrun files: {}\n".format(len(unsubmitted_files))
    for count_i,i in enumerate(unsubmitted_files):
        print "\t{}/{}".format(unsubmitted_paths[count_i],i)

    print "\nNumber of Incomplete files: {}\n".format(len(incomplete_files))
    for count_i,i in enumerate(incomplete_files):
        print "\t{}/{}".format(incomplete_paths[count_i],i)

    print "\nNumber of Complete files: {}\n".format(len(Current_files)-len(unsubmitted_files)-len(incomplete_files))
    
    quit()
#     # If no viable jobs were discovered then quit
#     if len(input_files) == 0:
#         print "No jobs in need of running, exiting..."
#         quit()

#     # Calculate the number of separate jobs to be submitted
#     N_bundles = int(ceil(float(len(input_files))/float(args.size)))
    
#     # Bundle the jobs
#     bundled_files = [[] for i in range(N_bundles) ]
#     bundled_paths = [[] for i in range(N_bundles) ]
#     for i in range(N_bundles):
#         bundled_files[i] = input_files[i*args.size:(i+1)*args.size]
#         bundled_paths[i] = input_paths[i*args.size:(i+1)*args.size]

#     # Create input files and submit each bundle
#     for n in range(len(bundled_files)):
        
#         # Set input_files and input_paths to point towards the bundled_files and bundled_paths sublists
#         # NOTE: the reuse of variable names (input_files, input_paths) is just a convenience since the following loops weren't written for the bundled feature.
#         input_files = bundled_files[n]
#         input_paths = bundled_paths[n]

#         # Initialize working variable for the number of sub jobs being submitted, total cores, total nodes, and jobs per node
#         N_jobs = len(input_files)
#         N_cores = N_jobs*args.procs
#         N_nodes = int(ceil(float(N_cores)/float(args.ppn)))
#         N_jpn   = int(args.ppn/args.procs)

#         # Begin making input files for the bundled submission
#         with open('{}.{}.submit'.format(args.outputname,n),'w') as f:
#             f.write("#PBS -q {}\n".format(args.queue))
#             f.write("#PBS -l mppwidth={}\n".format(N_cores))
#             if min_flag == 0:
#                 f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
#             elif min_flag == 1:
#                 f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))
#             f.write("#PBS -N {}.{}\n".format(args.outputname,n))
#             f.write("#PBS -j oe\n\n")
#             f.write("# load ccm openmpi module\nmodule load ccm\n\n")
#             f.write("# cd into the submission directory\n")
#             f.write("cd {}\n\n".format(working_dir))
#             f.write("# {} script that spawns parallel jobs on each node\n".format(args.outputname+'.spawn.sh'))
#             f.write("ccmrun ./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

#         subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

#         with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
#             f.write('#!/bin/bash\n\n')
#             f.write("#get the compute nodes allocated to the CCM job\n")
#             f.write("ccm_nodes=`sort -u $PBS_NODEFILE`\n\n")
#             f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
#             f.write("j=0\n")
#             f.write("for node in $ccm_nodes; do\n")
#             f.write("    let j=$j+1\n")
#             f.write("    ssh -n -T $node $PBS_O_WORKDIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
#     #        f.write("    ssh -n -T $node {}/{}.node.${{j}}.sh &\n".format(working_dir,args.outputname))
#             f.write("done\n")
#             f.write("wait\n")

#         subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

#         for i in range(0,N_nodes):
#             with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
#                 f.write('\n# cd into the submission directory\n')
#                 f.write('cd {}\n\n'.format(working_dir))
#                 f.write('# Submit jobs to this node\n')
#                 for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
#                     f.write("cd {}\n".format(j))
#                     f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
#                     f.write("cd {}\n\n".format(working_dir))
#                 f.write("wait\n")
#             subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

#         subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)
        
if __name__ == "__main__":
   main(sys.argv[1:])
