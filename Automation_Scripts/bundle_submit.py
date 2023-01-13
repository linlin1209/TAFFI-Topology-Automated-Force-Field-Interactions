#!/bin/env python                                                                                                                                                              
import sys,argparse,os,ast,re,fnmatch,matplotlib,subprocess,shutil
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

    parser.add_argument('-d', dest='path', default=None,
                        help = 'The program operates on all files discovered during a directory walk that match the -f argument. Optionally, a directory name or any string can also be supplied via this '+\
                               'argument and only files matching -f whose directory string includes -d will be operated on. (default: "")')

    parser.add_argument('-d_exact', dest='path_exact', default=None,
                        help = 'The program operates on all files discovered during a directory walk that match the -f argument. Optionally, a directory name or any string can also be supplied via this '+\
                               'argument and only files matching -f whose directory string includes -d_exact, exactly (no partial match as for -d), will be operated on. (default: "")')

    parser.add_argument('-p', dest='procs', default=1,
                        help = 'Specifies the number of processors for each job (default: 1; max 8)')
                        
    parser.add_argument('-o', dest='outputname', default="orca_mass",
                        help = 'Specifies the job name (default: orca_mass)')

    parser.add_argument('-t', dest='walltime', default=4,
                        help = 'Specifies the walltime for each job (default: 48, hours by default, if Xmin is used then the argument will be interpretted in minutes)')

    parser.add_argument('-q', dest='queue', default='standby',
                        help = 'Specifies the queue for the job (default: ccm_queue; see NERSC webpage for valid options)')

    parser.add_argument('-ppn', dest='ppn', default=20,
                        help = 'Specifies the number of processors per node on the cluster architecture. the -ppn %% -p should equal zero. (default: 24)') 

    parser.add_argument('-size', dest='size', default=100,
                        help = 'Specifies the number of calculations to bundle per job (default: 100)')

    parser.add_argument('--overwrite', dest='overwrite', default=0, action='store_const', const=1,
                        help = 'When set, if previous run data is discovered in the submission folders it is deleted. (default: off)')

    parser.add_argument('--resubmit', dest='resubmit', default=0, action='store_const', const=1,
                        help = 'When set, if previous run data is discovered it is checked for completion. Any jobs without run data or incomplete run data are resubmitted (default: off)')

    parser.add_argument('--no_geo_swap', dest='no_geo_swap', default=0, action='store_const', const=1,
                        help = 'This flag only affects the behavior of the program when it is running in --resubmit mode. The default behavior is to replace the original input geometry '+\
                               'block with the intermediate job geometry, if it is available (for example, this is desirable for restarting geometry optimizations). When this flag is present '+\
                               'the program resubmits the job with the original run geometry. (default: off)')

    parser.add_argument('--no_cleanup', dest='no_cleanup', default=False, action='store_const', const=True,
                        help = 'This flag only affects the behavior of the program when it is running in --resubmit mode. The default behavior is to remove all old files except *.xyz\'s '+\
                               'When this flag is present all of the old files are left in place. (default: off)')

    parser.add_argument('--new_folder', dest='new_folder', default=0, action='store_const', const=1,
                        help = 'This flag only affects the behavior of the program when it is running in --resubmit mode. The default behavior is to replace the original input geometry '+\
                               'block in the original input file with the intermediate job geometry if it is available (for example, this is desirable for restarting geometry optimizations). When this flag is present '+\
                               'the program creates a new subfolder in the original run directory called "RESTART," copies the original input file, and replaces the geometry in the copy. This avoids overwriting the original '+\
                               'input file. (default: off)')

    parser.add_argument('-path_to_exe',dest='path_to_exe', default="/depot/bsavoie/apps/orca_4_1_2/orca",
                        help = 'specifies the shell submission script to call.')

    parser.add_argument('-path_to_mpi',dest='path_to_mpi', default=None,
                        help = 'specifies the openmpi path containing the bin and lib folders. Each sched option has a default value that can be replaced if this is assigned. (defualt: None)')

    parser.add_argument('-sched',dest='scheduler', default="slurm-halstead",
                        help = 'specifies the scheduler for the cluster ("torque" "slurm") (default: "slurm")')

    parser.add_argument('--silent',dest='silent', default=False, const=True, action='store_const',
                        help = 'Shortcircuits all script related print statements (scheduler will probably still echo the job id)')

    args=parser.parse_args()
    if type(args.walltime) == str and "min" in args.walltime: args.walltime = int(args.walltime.split('min')[0]); min_flag = 1
    else: args.walltime = int(args.walltime); min_flag = 0
    args.procs = int(args.procs)
    args.ppn = int(args.ppn)
    args.size = int(args.size)
    args.scheduler = args.scheduler.lower()
    Filename = args.Filename
    working_dir = os.getcwd()

    # Check that the number of processors per job divides into the number of processors per node.
    if args.ppn % args.procs != 0:
        print("ERROR: the -ppn % -p must be zero to ensure that jobs aren't split across nodes. Exiting...")
        quit()

    # Create a dictionary from the filenames, where each dictionary key corresponds to a filename and each entry is a list
    # of subfiles to be processed as a batch. e.g., molecule.in might show up in twenty subfolders. molecule.in would end
    # up as a key, and the list of individual instances of molecule.in would constitute the entry.
    Files = {}
    if args.path is None and args.path_exact is None:
        Files[Filename] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk('.') for f in filenames if fnmatch.fnmatch(f,Filename) ]
    elif args.path_exact is not None:
        args.path_exact = args.path_exact.split()
        Files[Filename] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk('.') for f in filenames if (fnmatch.fnmatch(f,Filename) and True in [ i in dp.split('/') for i in args.path_exact ]) ]
    elif args.path is not None:
        args.path = args.path.split()
        Files[Filename] = [ os.path.join(dp, f) for dp, dn, filenames in os.walk('.') for f in filenames if (fnmatch.fnmatch(f,Filename) and True in [ i in dp for i in args.path ]) ]

    #################################################################
    # Iterate over all discovered input files, cd to the containing #
    # directory and check if the job already has an output file     #
    #################################################################    
    input_files = []
    input_paths = []
    
    for i in list(Files.keys()):

        # Sort the filenames as string
        Current_files = sorted(Files[i])

        if args.silent is False:        
            print("{}".format("#"*80))
            print("# {:^76s} #".format("PROCESSING THE FOLLOWING FILES"))
            print("{}".format("#"*80))
            for j in Current_files:
                print(j)

        # Iterate over the discovered input files
        for j in Current_files:

            # Change to the directory holding the current file
            path_to_file = working_dir+'/'+'/'.join(j.split('/')[0:-1])
            if path_to_file != '':
                os.chdir(path_to_file)

            # Save the submission input files and paths to lists
            current_name = j.split('/')[-1]
            if current_name.split('.')[0]+".out" not in os.listdir('.'):
                input_files += [current_name]
                input_paths += [path_to_file]

            # Check for job completion if resubmit option is toggled
            elif args.resubmit == 1:
                completed_flag = 0
                with open(current_name.split('.')[0]+".out",'r') as f:
                    for lines in f:
                        if "****ORCA TERMINATED NORMALLY****" in lines:
                            completed_flag = 1

                # Add incompleted files to the submission lists and clean up previous run files
                if completed_flag == 0:

                    # --new_folder protocols
                    if args.new_folder == 1:

                        # If the RESTART folder already exists then avoid the resubmission
                        if os.path.isdir(path_to_file+'/RESTART'):
                            continue
                        # Else, copy the run files into the RESTART folder
                        else:
                            os.mkdir(path_to_file+'/RESTART')
                            shutil.copy(path_to_file+'/'+current_name,path_to_file+'/RESTART/'+current_name)
                            path_to_file = path_to_file
                            current_name = 'RESTART/'+current_name
                        
                    # Swap the geometry of the input file with the most up-to-date geometry from the previous calculation
                    # (when running the program with the --no_swap_geo flag, this action is avoided.
                    if args.no_geo_swap == 0:

                        # Parse the run name from the file
                        geo_name=[]
                        with open(current_name,'r') as f:
                            for lines in f:
                                if "%base" in lines:
                                    geo_name = lines.split()[1].split("\"")[1]+".xyz"
                                    break

                        # If an intermediate geometry exists, check for completeness and replace the input geometry specification
                        if geo_name != [] and os.path.isfile(geo_name):

                            # read in intermediate geometry
                            Geo,Elements = scrape_xyz(geo_name)

                            # Replace the geometry block in the original input file
                            # NOTE: scrape_xyz returns empty lists if a complete geometry wasn't discovered. This is a safeguard against corrupt files
                            if Geo != []:
                                if args.silent is False:                                
                                    print("replacing geometry in {}".format(path_to_file+'/'+current_name))
                                replace_geo(current_name,Geo,Elements)

                    # Add the incomplete filename and path to the submission lists
                    input_files += [current_name]
                    input_paths += [path_to_file]                    

                    # Clean up old run data (excluding .xyz files)
                    if args.new_folder == 0 and args.no_cleanup is False:
                        files = [ k for k in os.listdir('.') if k != j.split('/')[-1] and '.xyz' not in k ]
                        for k in files:
                            os.remove(k)
                        
            # When running in overwrite mode, the job is resubmitted regardless of the existence/completeness of output data
            elif args.overwrite == 1:

                # Clean up old run data
                files = [ k for k in os.listdir('.') if k != j.split('/')[-1] ]
                for k in files:
                    os.remove(k)

                # Add the input file and path to the submission lists
                input_files += [current_name]
                input_paths += [path_to_file]
                                    
            # Skip any that already have output files present
            else:
                if args.silent is False:
                    print("Skipped file {} because output was already found".format(current_name))
            os.chdir(working_dir)

    # If no viable jobs were discovered then quit
    if len(input_files) == 0:
        if args.silent is False:
            print("No jobs in need of running, exiting...")
        quit()

    # Insert escape characters for ( and )
    for i in [ "(", ")" ]:
        input_files = [ _.replace(i,"\{}".format(i)) for _ in input_files ]
        input_paths = [ _.replace(i,"\{}".format(i)) for _ in input_paths ]

    # Calculate the number of separate jobs to be submitted
    N_bundles = int(ceil(float(len(input_files))/float(args.size)))
    
    # Bundle the jobs
    bundled_files = [[] for i in range(N_bundles) ]
    bundled_paths = [[] for i in range(N_bundles) ]
    for i in range(N_bundles):
        bundled_files[i] = input_files[i*args.size:(i+1)*args.size]
        bundled_paths[i] = input_paths[i*args.size:(i+1)*args.size]

    # Create input files and submit each bundle
    for n in range(len(bundled_files)):
        
        # Set input_files and input_paths to point towards the bundled_files and bundled_paths sublists
        # NOTE: the reuse of variable names (input_files, input_paths) is just a convenience since the following loops weren't written for the bundled feature.
        input_files = bundled_files[n]
        input_paths = bundled_paths[n]

        # Initialize working variable for the number of sub jobs being submitted, total cores, total nodes, and jobs per node
        N_jobs = len(input_files)
        N_cores = N_jobs*args.procs
        N_nodes = int(ceil(float(N_cores)/float(args.ppn)))
        N_jpn   = int(args.ppn/args.procs)

        # Submission on ORNL-titan
        if args.scheduler == "torque-titan":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                f.write("#PBS -N {}.{}\n".format(args.outputname,n))
                f.write("#PBS -A chm114\n")
                f.write("#PBS -l nodes={}\n".format(N_nodes))
                if min_flag == 0:
                    f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))
                if args.queue == 'debug':
                    f.write("#PBS -q debug\n")
                f.write("#PBS -S /bin/sh\n")
                f.write("#PBS -o {}.out\n".format(args.outputname))
                f.write("#PBS -e {}.err\n\n".format(args.outputname))

                f.write("#Setting OPENMPI paths and variables here:\n")
                f.write("export PATH=$PATH:/ccs/proj/chm114/openmpi/3.0.0/bin\n")
                f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/proj/chm114/openmpi/3.0.0/lib\n")
#                f.write("export OMP_NUM_THREADS={}\n".format(args.procs))
                f.write("export OMPI_MCA_btl=self,tcp,sm\n")
#                f.write("export OMPI_UNIVERSE_SIZE={:<d}\n\n".format(len(input_files)))

                f.write("# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                # aprun based submission
                f.write('# Start submitting jobs\n')
                for count_i,i in enumerate(input_files):
                    f.write("cd {}\n".format(input_paths[count_i]))
                    f.write("aprun -n 1 -d {} {} {} >> {} &\n".format(args.procs,args.path_to_exe,i,'.'.join([ j for j in i.split('.')[:-1] ])+'.out'))
                    f.write("cd {}\n\n".format(working_dir))
                f.write("wait\n")

#                 # Wraprun based submission
#                 f.write('# Start submitting jobs\n')
#                 f.write('module load python wraprun\n')
#                 f.write('wraprun ')
#                 for count_i,i in enumerate(input_files):
#                     if count_i == len(input_files)-1:
#                         f.write('-n 1 -d {} --w-cd {} {} {} \n'.format(args.procs,input_paths[count_i],args.path_to_exe,i))
#                     else:
#                         f.write('-n 1 -d {} --w-cd {} {} {} : '.format(args.procs,input_paths[count_i],args.path_to_exe,i))
#                 f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)
            subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Submission on ITAP-halstead
        if args.scheduler == "torque-halstead":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                # Boilerplate
                f.write("#PBS -N {}.{}\n".format(args.outputname,n))
                f.write("#PBS -l nodes={}:ppn={}\n".format(N_nodes,int(args.ppn)))
                if min_flag == 0:
                    f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))                    
                f.write("#PBS -q {}\n".format(args.queue))
                f.write("#PBS -S /bin/sh\n")
                f.write("#PBS -o {}.out\n".format(args.outputname))
                f.write("#PBS -e {}.err\n\n".format(args.outputname))

                # Set up paths
                if args.path_to_mpi is None:
                    args.path_to_mpi = "/depot/bsavoie/apps/openmpi/3.0.0"
                f.write("#Setting OPENMPI paths and variables here:\n")
                f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                # Print diagnostics
                f.write("\n# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                # Execute jobs spawn script
                f.write("# {} is a script that spawns parallel jobs on each node\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            # Write job spawn script resonsible for spawning jobs across the allocated nodes
            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the job\n")
                f.write("job_nodes=`sort -u $PBS_NODEFILE`\n\n")
                f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                f.write("j=0\n")
                f.write("for node in $job_nodes; do\n")
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write("    let j=$j+1\n")
                f.write("    ssh -n -T $node $PBS_O_WORKDIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                f.write("done\n")
                f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            # Create the individual run scripts that are executed on each node
            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:

                    # Set up paths
                    f.write("#Setting OPENMPI paths and variables here:\n")
                    f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                    f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)
        # Submission on ITAP-halstead
        if args.scheduler == "slurm-halstead":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                # Boilerplate
                f.write("#!/bin/bash \n\n")
                f.write("#SBATCH --job-name {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH -N {}\n".format(N_nodes))
                f.write("#SBATCH -n {}\n".format(N_nodes*int(args.ppn))) # total cpu count, not per node
                if min_flag == 0:
                    f.write("#SBATCH -t {}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#SBATCH -t 00:{}:00\n".format(args.walltime))                    
                f.write("#SBATCH -A {}\n".format(args.queue))
                f.write("#SBATCH -o {}.out\n".format(args.outputname))
                f.write("#SBATCH -e {}.err\n\n".format(args.outputname))

                # Set up paths
                if args.path_to_mpi is None:
                    args.path_to_mpi = "/depot/bsavoie/apps/openmpi/3.0.0"
                f.write("#Setting OPENMPI paths and variables here:\n")
                f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                # Print diagnostics
                f.write("\n# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                # Execute jobs spawn script
                f.write("# {} is a script that spawns parallel jobs on each node\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            # Write job spawn script resonsible for spawning jobs across the allocated nodes
            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the job\n")
                f.write("echo $SLURM_JOB_NODELIST > nodefile\n")
                f.write("job_nodes=`sort -u nodefile`\n")
                f.write("rm nodefile\n\n")
                if N_nodes > 1:
                   f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                   f.write("j=0\n")
                   f.write("for node in $job_nodes; do\n")
                   f.write('    echo -e "{}"\n'.format("*"*100))
                   f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                   f.write('    echo -e "{}"\n'.format("*"*100))
                   f.write("    let j=$j+1\n")
                   f.write("    if [[ $node == $(hostname) ]]; then\n")
                   f.write("        $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                   f.write("    else\n") 
                   f.write("        ssh -n -T $node $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                   f.write("    fi\n")
                   f.write("done\n")
                   f.write("wait\n")
                # ssh is not required if only one node is being used. 
                else:
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write('echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${SLURM_JOB_NODELIST}"))
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write("$SLURM_SUBMIT_DIR/{}.{}.node.1.sh &\n".format(args.outputname,n))
                    f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            # Create the individual run scripts that are executed on each node
            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:

                    # Set up paths
                    f.write("#Setting OPENMPI paths and variables here:\n")
                    f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                    f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("sh {} &\n".format(input_files[input_paths.index(j)]))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("sbatch {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Submission on ITAP-halstead
        if args.scheduler == "torque-halstead-new":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                f.write("#PBS -N {}.{}\n".format(args.outputname,n))
                f.write("#PBS -l nodes={}\n".format(N_nodes))
                if min_flag == 0:
                    f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))                    
                f.write("#PBS -q {}\n".format(args.queue))
                f.write("#PBS -S /bin/sh\n")
                f.write("#PBS -o {}.out\n".format(args.outputname))
                f.write("#PBS -e {}.err\n\n".format(args.outputname))

                f.write("#Setting OPENMPI paths and variables here:\n")
                f.write("export PATH=$PATH:/depot/bsavoie/apps/openmpi/3.0.0/bin\n")
                f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/depot/bsavoie/apps/openmpi/3.0.0/lib\n\n")

                f.write("\n# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                f.write("# {} is a script that spawns parallel jobs on each node\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the job\n")
                f.write("job_nodes=`sort -u $PBS_NODEFILE`\n\n")
                f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                f.write("j=0\n")
                f.write("for node in $job_nodes; do\n")
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write("    let j=$j+1\n")
                f.write("    ssh -n -T $node $PBS_O_WORKDIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                f.write("done\n")
                f.write("wait\n")

            count = 0
            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
                    f.write("\n#Setting OPENMPI paths and variables here:\n")
                    f.write("export PATH=$PATH:/depot/bsavoie/apps/openmpi/3.0.0/bin\n")
                    f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/depot/bsavoie/apps/openmpi/3.0.0/lib\n\n")
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    tmp_folders    = []
                    move_back_locs = []
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("mkdir {}\n".format(count))
                        f.write("cd {}\n".format(count))
                        f.write("cp {}/{} .\n".format(j,input_files[input_paths.index(j)]))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd ..\n\n".format(working_dir))
                        move_back_locs += [j]
                        tmp_folders += [count]
                        count += 1
                    f.write("#Wait until the jobs complete\n")
                    f.write("wait\n\n")
                    f.write("# Copy the data back to the original folder(s)\n")
                    for count_j,j in enumerate(tmp_folders):
                        f.write("mv {}/* {}/. \n".format(j,move_back_locs[count_j]))
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Submission on xsede-comet version 2
        if args.scheduler == "slurm-comet2":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                # Boilerplate
                f.write("#!/bin/bash \n\n")
                f.write("#SBATCH --job-name {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH -N {}\n".format(N_nodes))
                f.write("#SBATCH -n {}\n".format(N_nodes*int(args.ppn))) # total cpu count, not per node
                if min_flag == 0:
                    f.write("#SBATCH -t {}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#SBATCH -t 00:{}:00\n".format(args.walltime))                    
                # XXX hard coded account
                f.write("#SBATCH -A prd180\n".format(args.queue))
                f.write("#SBATCH -o {}.out\n".format(args.outputname))
                f.write("#SBATCH -e {}.err\n\n".format(args.outputname))
                f.write("#SBATCH --export=ALL\n")
                f.write("#SBATCH --partition={}\n".format(args.queue))

                # Set up paths
                if args.path_to_mpi is None:
                    args.path_to_mpi = "/depot/bsavoie/apps/openmpi/3.0.0"
                f.write("#Setting OPENMPI paths and variables here:\n")
                #f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                #f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                # Print diagnostics
                f.write("\n# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                # Execute jobs spawn script
                f.write("# {} is a script that spawns parallel jobs on each node\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            # Write job spawn script resonsible for spawning jobs across the allocated nodes
            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the job\n")
                #f.write("echo $SLURM_JOB_NODELIST > nodefile\n")
                #f.write("job_nodes=`sort -u nodefile`\n")
                #f.write("rm nodefile\n\n")
                if N_nodes > 1:
                   # XXX need modification when multiples nodes are used
                   f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                   f.write("j=0\n")
                   f.write("for node in $job_nodes; do\n")
                   f.write('    echo -e "{}"\n'.format("*"*100))
                   f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                   f.write('    echo -e "{}"\n'.format("*"*100))
                   f.write("    let j=$j+1\n")
                   f.write("    if [[ $node == $(hostname) ]]; then\n")
                   f.write("        $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                   f.write("    else\n") 
                   f.write("        ssh -n -T $node $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                   f.write("    fi\n")
                   f.write("done\n")
                   f.write("wait\n")
                # ssh is not required if only one node is being used. 
                else:
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write('echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${SLURM_JOB_NODELIST}"))
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write("$SLURM_SUBMIT_DIR/{}.{}.node.1.sh &\n".format(args.outputname,n))
                    f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            # Create the individual run scripts that are executed on each node
            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:

                    # Set up paths
                    f.write("#Setting OPENMPI paths and variables here:\n")
                    #f.write("export PATH={}/bin:$PATH\n".format(args.path_to_mpi))
                    #f.write("export LD_LIBRARY_PATH={}/lib:$LD_LIBRARY_PATH\n\n".format(args.path_to_mpi))
                    f.write("module load openmpi_ib/3.1.4\n")
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("sh {} &\n".format(input_files[input_paths.index(j)]))
                        #f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("sbatch {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Submission on ORNL-rhea
        if args.scheduler == "torque-rhea":

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:

                f.write("#PBS -N {}.{}\n".format(args.outputname,n))
                f.write("#PBS -A chm114\n")
                f.write("#PBS -l nodes={}\n".format(N_nodes))
                if min_flag == 0:
                    f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))
                if args.queue == 'debug':
                    f.write("#PBS -q debug\n")
                f.write("#PBS -S /bin/sh\n")
                f.write("#PBS -o {}.out\n".format(args.outputname))
                f.write("#PBS -e {}.err\n\n".format(args.outputname))

                f.write("#Setting OPENMPI paths and variables here:\n")
                f.write("export PATH=$PATH:/ccs/proj/chm114/openmpi/3.0.0/bin\n")
                f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/proj/chm114/openmpi/3.0.0/lib\n")
                f.write("export OMPI_MCA_btl=self,tcp,sm\n")

                f.write("\n# cd into the submission directory\n")
                f.write("cd {}\n".format(working_dir))
                f.write("echo Working directory is ${}\n".format(working_dir))
                f.write("echo Running on host `hostname`\n")
                f.write("echo Time is `date`\n\n")

                f.write("# {} is a script that spawns parallel jobs on each node\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the job\n")
                f.write("job_nodes=`sort -u $PBS_NODEFILE`\n\n")
                f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                f.write("j=0\n")
                f.write("for node in $job_nodes; do\n")
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write("    let j=$j+1\n")
                f.write("    ssh -n -T $node $PBS_O_WORKDIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                f.write("done\n")
                f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
                    f.write("\n#Setting OPENMPI paths and variables here:\n")
                    f.write("export PATH=$PATH:/ccs/proj/chm114/openmpi/3.0.0/bin\n")
                    f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ccs/proj/chm114/openmpi/3.0.0/lib\n")
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Generate batch scripts and submission commands for torque clusters @ NERSC
        if args.scheduler == 'torque-NERSC':

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:
                f.write("#PBS -q {}\n".format(args.queue))
                f.write("#PBS -l mppwidth={}\n".format(N_cores))
                if min_flag == 0:
                    f.write("#PBS -l walltime={}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#PBS -l walltime=00:{}:00\n".format(args.walltime))
                f.write("#PBS -N {}.{}\n".format(args.outputname,n))
                f.write("#PBS -j oe\n\n")
                f.write("# load ccm openmpi module\nmodule load ccm\n\n")
                f.write("# cd into the submission directory\n")
                f.write("cd {}\n\n".format(working_dir))
                f.write("# {} script that spawns parallel jobs on each node\n".format(args.outputname+'.spawn.sh'))
                f.write("ccmrun ./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the CCM job\n")
                f.write("ccm_nodes=`sort -u $PBS_NODEFILE`\n\n")
                f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))
                f.write("j=0\n")
                f.write("for node in $ccm_nodes; do\n")
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                f.write('    echo -e "{}"\n'.format("*"*100))
                f.write("    let j=$j+1\n")
                f.write("    ssh -n -T $node $PBS_O_WORKDIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                f.write("done\n")
                f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)

            subprocess.call("qsub {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Generate batch scripts and submission commands for slurm clusters
        if args.scheduler == 'slurm-cori':

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:
                f.write("#!/bin/bash -l\n\n")
                f.write("#SBATCH --partition {}\n".format(args.queue))
                f.write("#SBATCH -N {}\n".format(N_nodes))       # Set total number of nodes                       
                if min_flag == 0:
                    f.write("#SBATCH -t {}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#SBATCH -t 00:{}:00\n".format(args.walltime))
                f.write("#SBATCH -J {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH -o {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH --ccm\n\n")
                f.write("module swap openmpi/1.6.5\n")
                f.write("export OMP_NUM_THREADS={}\n".format(args.procs))
                f.write("# cd into the submission directory\n")
                f.write("cd {}\n\n".format(working_dir))
                f.write("# {} script that spawns parallel jobs on each node\n".format(args.outputname+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the CCM job\n")
                f.write("job_nodes=$(/usr/common/usg/bin/gen_nodelist.sh)\n")
                f.write("echo -e ${job_nodes}\n\n")                
                f.write('echo -e "SLURM_NODELIST: ${SLURM_NODELIST}"\n\n')

                if N_nodes > 1:
                    f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))                
                    f.write("j=0\n")
                    f.write("for node in $(/usr/common/usg/bin/gen_nodelist.sh); do \n")
                    f.write('    echo -e "{}"\n'.format("*"*100))
                    f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                    f.write('    echo -e "{}"\n'.format("*"*100))
                    f.write("    let j=$j+1\n")
                    f.write("    if [[ $node == $(hostname) ]]; then\n")
                    f.write("        $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                    f.write("    else\n") 
                    f.write("        ssh -o StrictHostKeyChecking=no -n -T $node $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                    f.write("    fi\n")
                    f.write("done\n")
                    f.write("wait\n")

                # ssh is not required if only one node is being used. 
                else:
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write('echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${SLURM_NODELIST}"))
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write("$SLURM_SUBMIT_DIR/{}.{}.node.1.sh &\n".format(args.outputname,n))
                    f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)
            subprocess.call("sbatch {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Generate batch scripts and submission commands for slurm clusters
        if args.scheduler == 'slurm-edison':

            # Begin making input files for the bundled submission
            with open('{}.{}.submit'.format(args.outputname,n),'w') as f:
                f.write("#!/bin/bash -l\n\n")
                f.write("#SBATCH --partition {}\n".format(args.queue))
                f.write("#SBATCH -N {}\n".format(N_nodes))       # Set total number of nodes                       
                if min_flag == 0:
                    f.write("#SBATCH -t {}:00:00\n".format(args.walltime))
                elif min_flag == 1:
                    f.write("#SBATCH -t 00:{}:00\n".format(args.walltime))
                f.write("#SBATCH -J {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH -o {}.{}\n".format(args.outputname,n))
                f.write("#SBATCH --ccm\n\n")
                f.write("module swap openmpi-ccm/1.6.5\n")
                f.write("export OMP_NUM_THREADS={}\n".format(args.procs))
                f.write("# cd into the submission directory\n")
                f.write("cd {}\n\n".format(working_dir))
                f.write("# {} script that spawns parallel jobs on each node\n".format(args.outputname+'.spawn.sh'))
                f.write("./{}\nwait\n".format(args.outputname+'.'+str(n)+'.spawn.sh'))

            subprocess.call("chmod 777 {}.{}.submit".format(args.outputname,n), shell=True)

            with open('{}.{}.spawn.sh'.format(args.outputname,n),'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write("#get the compute nodes allocated to the CCM job\n")
                f.write("job_nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST)\n")
                f.write("echo -e ${job_nodes}\n\n")                
                f.write('echo -e "SLURM_NODELIST: ${SLURM_NODELIST}"\n\n')

                if N_nodes > 1:
                    f.write("#ssh to each compute node and execute {}.node.${{j}}.sh where j is the node.\n".format(args.outputname))                
                    f.write("j=0\n")
                    f.write("for node in $(/usr/common/usg/bin/gen_nodelist.sh); do \n")
                    f.write('    echo -e "{}"\n'.format("*"*100))
                    f.write('    echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${node}"))
                    f.write('    echo -e "{}"\n'.format("*"*100))
                    f.write("    let j=$j+1\n")
                    f.write("    if [[ $node == $(hostname) ]]; then\n")
                    f.write("        $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                    f.write("    else\n") 
                    f.write("        ssh -o StrictHostKeyChecking=no -n -T $node $SLURM_SUBMIT_DIR/{}.{}.node.${{j}}.sh &\n".format(args.outputname,n))
                    f.write("    fi\n")
                    f.write("done\n")
                    f.write("wait\n")

                # ssh is not required if only one node is being used. 
                else:
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write('echo -e "* {:^96s} *"\n'.format("RUNNING ON NODE ${SLURM_NODELIST}"))
                    f.write('echo -e "{}"\n'.format("*"*100))
                    f.write("$SLURM_SUBMIT_DIR/{}.{}.node.1.sh &\n".format(args.outputname,n))
                    f.write("wait\n")

            subprocess.call("chmod 777 {}.{}.spawn.sh".format(args.outputname,n), shell=True)

            for i in range(0,N_nodes):
                with open('{}.{}.node.{}.sh'.format(args.outputname,n,i+1),'w') as f:
                    f.write('\n# cd into the submission directory\n')
                    f.write('cd {}\n\n'.format(working_dir))
                    f.write('# Submit jobs to this node\n')
                    for j in input_paths[i*N_jpn:i*N_jpn+N_jpn]:
                        f.write("cd {}\n".format(j))
                        f.write("{} {} >> {} &\n".format(args.path_to_exe,input_files[input_paths.index(j)],'.'.join([ k for k in input_files[input_paths.index(j)].split('.')[:-1] ])+'.out'))
                        f.write("cd {}\n\n".format(working_dir))
                    f.write("wait\n")
                subprocess.call("chmod 777 {}.{}.node.{}.sh".format(args.outputname,n,i+1), shell=True)
            subprocess.call("sbatch {}".format(args.outputname+'.'+str(n)+'.submit'), shell=True)

        # Generate batch scripts and submission commands for tachus - Miller group cluster
        # Submission here consists of single jobs, with the batch scripts generated in the
        # same folder as the input file and submission also occuring in that folder
        if args.scheduler == 'slurm-tachus':
            
            current_dir=os.getcwd()
            for count_i,i in enumerate(input_paths):
                os.chdir(i)                                                
                base_name = '.'.join([ j for j in input_files[count_i].split('.')[:-1] ])

                # Write commands for the running the current job
                with open(base_name+'.submit','w') as f:
                    f.write("#!/bin/bash\n")
                    f.write("#\n")
                    f.write("#SBATCH --job-name={}\n".format(base_name))
                    f.write("#SBATCH --output={}.out\n".format(base_name))
                    f.write("#SBATCH --error={}.err\n".format(base_name))
                    f.write("#SBATCH --partition={}\n".format(args.queue))
                    f.write("#SBATCH --nodes=1\n")
                    f.write("#SBATCH --ntasks-per-node={}\n".format(args.procs))
                    if min_flag == 0:
                        f.write("#SBATCH -t {}:00:00\n".format(args.walltime))
                    elif min_flag == 1:
                        f.write("#SBATCH -t 00:{}:00\n".format(args.walltime))

                    f.write("\n# Load MPI module\n")
#                    f.write("module load mpi/openmpi/1.8.1/gcc44\n\n")
                    f.write("module load mpi/openmpi/1.8.3/intel15\n\n")

                    f.write("#Setting OPENMPI paths here:\n")
                    f.write('export PATH=$PATH:/share/apps/openmpi/1.8.3/intel15/bin\n')
                    f.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/openmpi/1.8.3/intel15/lib\n')

                    f.write("\n# Write out some information on the job\n")
                    f.write('echo Running on hosts: $SLURM_NODELIST\n')
                    f.write('echo Running on $SLURM_NNODES nodes.\n')
                    f.write('echo "Running on \$SLURM_NPROCS processors."\n')
                    f.write('echo "Current working directory is `pwd`"\n\n')
                    f.write('echo "Copying input file to scratch..."\n')
                    f.write('mkdir /scratch/$SLURM_JOB_NAME.$SLURM_JOB_ID\n')
                    f.write('cp {} /scratch/$SLURM_JOB_NAME.$SLURM_JOB_ID/.\n'.format(input_files[count_i]))
                    f.write('cd /scratch/$SLURM_JOB_NAME.$SLURM_JOB_ID/.\n\n')

                    f.write("# Run simulation\n")
                    f.write("sleep 2\n")
                    f.write("{} {}\n".format(args.path_to_exe,input_files[count_i].split('/')[-1]))
                    f.write("rm *.gbw\n")
                    f.write("mv * {}/.\n".format(i))
                    f.write("cd {}\n".format(i))
                    f.write("rm -r /scratch/$SLURM_JOB_NAME.$SLURM_JOB_ID\n")
                    f.write("sleep 2\n")

                subprocess.call("chmod 777 {}.submit".format(base_name), shell=True)
                subprocess.call("sbatch {}.submit".format(base_name), shell=True)
                os.chdir(current_dir)

# read in intermediate geometry
def scrape_xyz(xyz_file):

    # Initialize atom counter
    counter = 0

    # Search the xyz file for element and coordinate data
    with open(xyz_file,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            # grab the atom count from the first line and initialize Geo and Element array/list
            if lc == 0:
                atom_count = int(fields[0])
                Geo = zeros([atom_count,3])
                Elements = ["X"]*atom_count
                continue

            # parse Geo and Element info
            if len(fields) > 0 and lc > 1:
                Geo[lc-2,:] = array([float(fields[1]),float(fields[2]),float(fields[3])])
                Elements[lc-2] = fields[0]
                counter += 1

    # Only return a full geometry (safeguard against corrupt files)
    if counter == atom_count:
        return Geo,Elements
    else:
        return [],[]

# replaces the geometry block of the input file
def replace_geo(current_name,Geo,Elements):

    # Read the old input file into memory
    with open(current_name,'r') as f:
        content = f.readlines()

    # Find the start and end of the geometry specification
    geo_flag = 0
    for lc,lines in enumerate(content):
        fields = lines.split()
        if len(fields) > 1 and "*" == fields[0] and "xyz" == fields[1]:
            geo_start = lc+1
            geo_flag = 1
            continue
        if geo_flag == 1 and len(fields) > 0 and fields[0] == "*":
            geo_end = lc
            break

    # Write new input file
    with open(current_name,'w') as f:
        for lines in content[:geo_start]:
            f.write(lines)
        for count_i,i in enumerate(Geo):
            f.write('  {:<20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(Elements[count_i],i[0],i[1],i[2]))
        for lines in content[geo_end:]:
            f.write(lines)
            
    return

if __name__ == "__main__":
   main(sys.argv[1:])
