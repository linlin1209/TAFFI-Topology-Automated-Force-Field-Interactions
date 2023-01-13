
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob
import numpy as np
from matplotlib import pyplot as plt
# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
from monitor_jobs import *
from copy import deepcopy

def main(argv):

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')

    # parse configuration dictionary (c)
    print("parsing benchmark configurations...")
    c = parse_configuration(parser.parse_args())

    # run geometry optimizations and frequency analysis
    print("\nentering normal mode phase...")
    run_normal_analysis(c)

    # run md
    print("\nentering molecular dynamics phase...")
    #run_md(c)

    quit()

# Wrapper function for the normal mode analysis
def run_normal_analysis(c):
    
    # Make the Normal_Modes 
    if os.path.isdir("Normal_Modes") is False:
        os.mkdir("Normal_Modes")

    # Loop over the xyz files, copy them, create canonical conformer, and generate orca input files
    print("\tgenerating input files for  geometry optimizations...")
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        if os.path.isdir("Normal_Modes/{}_freq".format(name)) is False:
            
            # transify the geometry
            subprocess.check_output(subproc_string("python {}/Lib/transify.py {}".format(c["taffi_path"],i)), encoding='utf8')

            # use openbabel to minimize the geometry via the UFF force-field
            output = subprocess.Popen(subproc_string("obminimize -ff uff -sd -c 1e-20 -n 100000 straightened.xyz  > {}.opt.xyz 2> /dev/null".format(name)), stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
            with open("{}.opt.xyz".format(name),'w') as f:
                for lines in output:
                    fields = lines.split()
                    if len(fields) > 0 and fields[0] == "WARNING:":
                        continue
                    else:
                        f.write(lines)
                    
            if os.path.isfile("{}.opt.xyz".format(name)) is False:
                print("ERROR: file {} failed at openbabel optimization...".format(i))
                quit()

            output = subprocess.Popen(subproc_string("python {}/FF_functions/xyz_to_orca_FOM.py {}.opt.xyz -p {} -o {}_freq.in".format(c["taffi_path"],name,c["fom_geoopt_procs"],name)),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
            os.mkdir("Normal_Modes/{}_freq".format(name))
            shutil.move("{}.opt.xyz".format(name),"Normal_Modes/{}_freq/{}.opt.xyz".format(name,name))
            shutil.move("{}_freq.in".format(name),"Normal_Modes/{}_freq/.".format(name))
            os.remove("straightened.xyz")

    # Submit the geometry optimizations
    jobids = []
    print("\trunning geometry optimizations...")
    #substring = "python {}/Automation_Scripts/orca_submit.py -f *freq.in -p {} -t {} -ppn {} -q {} -sched {} -o frequencies -path_to_exe {} -size {} --silent"
    substring = "python {}/Automation_Scripts/orca_submit.py -d Normal_Modes/ -p {} -t {} -ppn {} -q {} -sched {} -o frequencies -path_to_exe {} -size {} --silent"
    substring = substring.format(c["taffi_path"],c["fom_geoopt_procs"],c["fom_geoopt_wt"],c["fom_geoopt_ppn"],c["fom_geoopt_q"],c["fom_geoopt_sched"],c["orca_exe"],c["fom_geoopt_size"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,'lin1209')

    # Remove submission files
    for f in glob.glob("frequencies.*"):
        os.remove(f)

    # Calculate the classical modes
    jobids=[]
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Check for incomplete
        if os.path.isdir("Normal_Modes/{}_minimized".format(name)) is True and os.path.isfile("Normal_Modes/{}_minimized/minimized.xyz".format(name)) is False:
            shutil.rmtree("Normal_Modes/{}_minimized".format(name))

        # Run the minimization and mode calculation
        if os.path.isdir("Normal_Modes/{}_minimized".format(name)) is False:

            # Check the requisite data is in place. 
            if os.path.isdir("Normal_Modes/{}_freq".format(name)) is False or os.path.isfile("Normal_Modes/{}_freq/{}_freq.out".format(name,name)) is False:
                print("An error occured when trying to process {}_freq. Exiting...".format(name))
                quit()

            # Assemble the submission string for the batch submission file
            substring="module load gcc &> /dev/null\n"
            substring+="python {}/FF_functions/lammps_minimize_molecule.py {}_freq/geo_opt.xyz {} -o {}_minimized -lammps_exe {} -gens {} > /dev/null \n"
            substring=substring.format(c["taffi_path"],name_safe,c["fom_ff"],name_safe,c["lammps_exe"],c['gens'])                

            substring+="if [ ! -f {}_minimized/minimized.xyz ]; then\n"
            substring+="    echo 'An error occured during the minimization of {}.xyz'\n"
            substring+="    exit\n"
            substring+="fi\n\n"
            substring = substring.format(name_safe,name)
            substring+="# Calculate the FF normal modes, and create an aligned geometry\n"
            substring+="python {}/Parsers/parse_normal_modes.py -geo {}_minimized/minimized.xyz -FF {} -o {}_FF -out {}_freq/{}_freq.out -gens {} > /dev/null \n"
            substring+="python {}/Parsers/align_points.py {}_freq/geo_opt.xyz {}_minimized/minimized.xyz -o {}_aligned > /dev/null  \n"
            substring = substring.format(c["taffi_path"],name_safe,c["fom_ff"],name_safe,name_safe,name_safe,c['gens'],c["taffi_path"],name_safe,name_safe,name_safe)

            # Write the batch script and submit the job
            shell_submit(substring,"{}.submit".format(name_sub),1,1,c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}/{}".format(os.getcwd(),"Normal_Modes"))
            output = subprocess.Popen([c["sub_cmd"],"{}.submit".format(name_sub)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
            jobids += [ output.split()[-1] ]

    # Wait until the jobs complete
    monitor_jobs(jobids,'lin1209')

    # Remove the submission scripts
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        if os.path.isfile("{}.submit".format(name_sub)) is True:
            for f in glob.glob("{}.submit*".format(name_sub)):
                os.remove(f)

    # Plot the normal mode comparisons
    if os.path.isfile("Normal_Modes/nm.txt") is False:
        current = os.getcwd()
        os.chdir("Normal_Modes")
        output = subprocess.Popen("python {}/Parsers/plot_normal_modes.py > nm.log".format(c["taffi_path"]),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True, encoding='utf8').communicate()[0].strip("\r\n")
        os.chdir(current)

    return

# Wrapper function for the md and figure of merit analysis
def run_md(c):

    # Make the md directory if it doesn't already exist
    if os.path.isdir("FOM") is False:
        os.mkdir("FOM")
        
    # Copy the xyz files if necessary
    for i in c["xyz"]:
        if os.path.isfile("FOM/{}".format(i)) is False:
            shutil.copy(i,"FOM/{}".format(i))

    # Generate the single molecule and condensed phase input files
    print("\tgenerating input files for gas-phase and condensed-phase simulations...")
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)
        
        # Grab the number of atoms in this molecule
        with open(i,'r') as f:
            for lc,lines in enumerate(f):
                N_atoms = int(lines.split()[0])
                break

        # Calculate the number of molecules to be run in the condensed phase simulation
        N_mol = int(np.ceil(float(5000)/float(N_atoms))) 

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)
 
        # If the gas-phase simulation folder doesn't exist
        if os.path.isdir("FOM/{}_gas".format(name)) is False:
            os.mkdir("FOM/{}_gas".format(name))
            
            # Generate five trajectories at each temperature
            for t in T:
                if os.path.isdir("FOM/{}_gas/{}_K".format(name,t)) is False:
                    os.mkdir("FOM/{}_gas/{}_K".format(name,t))
                for j in range(1,c["fom_n_traj"]+1):
                    if os.path.isdir("FOM/{}_gas/{}_K/run_{}".format(name,t,j)) is False:
                        substring="python {}/FF_functions/gen_mol_traj.py {} {} -T {} -T_A {} -o FOM/{}_gas/{}_K/run_{} -mixing_rule {} -gens {}"
                        substring=substring.format(c["taffi_path"],i,c["fom_ff"],t,t,name,t,j,'wh',c['gens'])
                        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
        
        # If the condensed-phase simulation folder doesn't exist
        if os.path.isdir("FOM/{}_cond".format(name)) is False:
            os.mkdir("FOM/{}_cond".format(name))
            
        # Generate five trajectories at each temperature
        for t in T:
            if os.path.isdir("FOM/{}_cond/{}_K".format(name,t)) is False:
                os.mkdir("FOM/{}_cond/{}_K".format(name,t))
            for j in range(1,c["fom_n_traj"]+1):
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j)) is False:
                    substring="python {}/FF_functions/gen_md_for_sampling.py {} -FF {} -N {} -T {} -T_A 10 -t_A {} -t 0 -t_ext {} -o FOM/{}_cond/{}_K/run_{} -mixing wh --tail -q 0 -gens {}" 
                    substring=substring.format(c["taffi_path"],i,c["fom_ff"],N_mol,t,c["fom_t_a"],c["fom_t_ext"],name,t,j,c['gens'])
                    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")

    # Run the MD initializations
    print("\trunning gas phase md and condensed phase initializations...")
    jobids = []
    substring="python {}/Automation_Scripts/lammps_submit.py -f run*.in.init -d gas -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} -path_to_exe {} -o {} --silent"
    substring = substring.format(c["taffi_path"],c["fom_md_gas_procs"],c["fom_md_gas_ppn"],c["fom_md_gas_size"],c["fom_md_gas_q"],c["fom_md_gas_sched"],c["fom_md_gas_wt"],c["acct"],c["lammps_exe"],"gas_init")
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    substring="python {}/Automation_Scripts/lammps_submit.py -f run*.in.init -d cond -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} -path_to_exe {} -o {} --silent"
    substring = substring.format(c["taffi_path"],c["fom_md_cond_procs"],c["fom_md_cond_ppn"],c["fom_md_cond_size"],c["fom_md_cond_q"],c["fom_md_cond_sched"],c["fom_md_cond_wt"],c["acct"],c["lammps_exe"],"cond_init")
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,'lin1209')

    # Check that all initializations completed and attempt "fom_restarts" number of resubmissions
    print("\tchecking if resubmissions are necessary...")
    for r in range(c["fom_restarts"]):
        jobids = []
        incomplete = []
        # Loop over jobs
        for i in c["xyz"]:

            # Save base name and versions to variables
            name=".".join(i.split('.')[:-1])
            name_sub=name.replace('(','-')
            name_sub=name_sub.replace(')','-')
            name_safe=esc_str(name)

            # Get the simulation temperatures
            T = get_temps(c["fom_temps"],name)

            # Check the presence of *.end.data files (only generated if the job completes)
            for t in T:
                for j in range(1,c["fom_n_traj"]+1):                
                    
                    # Only attempt condensation resubmission
                    if os.path.isfile("FOM/{}_cond/{}_K/run_{}/run_{}.end.data".format(name,t,j,j)) is False:
                        incomplete += [ "FOM/{}_cond/{}_K/run_{}/run_{}.in.init".format(name,t,j,j) ]

                        substring="python {}/Automation_Scripts/lammps_submit.py -f run*.in.init -d {}_cond^{}_K^run_{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} -path_to_exe {} -o {} --overwrite --silent"
                        substring = substring.format(c["taffi_path"],name,t,j,c["fom_md_cond_procs"],c["fom_md_cond_ppn"],c["fom_md_cond_size"],c["fom_md_cond_q"],c["fom_md_cond_sched"],c["fom_md_cond_wt"],c["acct"],\
                                                     c["lammps_exe"],"init_{}_{}_{}".format(name_sub,t,j))
                        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                        output = str(output,'utf-8')
                        jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
                        
        # If jobs were incomplete, print diagnostic and wait on submission
        if len(incomplete) > 0:
            print("\tresubmitting the following initializations:\n\n\t\t{}".format("\n\t\t".join(incomplete)))
        
        monitor_jobs(jobids,'lin1209')
            
    # Check that all initializations completed
    print("\tchecking initializations for completion...")
    incomplete = []
    diffuse = []

    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Check the presence of *.end.data files (only generated if the job completes)
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):                
                if os.path.isfile("FOM/{}_gas/{}_K/run_{}/run_{}.end.data".format(name,t,j,j)) is False:
                    incomplete += [ "FOM/{}_gas/{}_K/run_{}/run_{}.in.init".format(name,t,j,j) ]                
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/run_{}.end.data".format(name,t,j,j)) is False:
                    incomplete += [ "FOM/{}_cond/{}_K/run_{}/run_{}.in.init".format(name,t,j,j) ]

                # Calculate density in cond thermo to ensure condensation
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/thermo.avg".format(name,t,j)) is True:
                    d = np.mean(get_cols("FOM/{}_cond/{}_K/run_{}/thermo.avg".format(name,t,j),[2],head=20)) 
                    if d<0.5:
                        diffuse += [ "FOM/{}_cond/{}_K/run_{}/run_{}.in.init".format(name,t,j,j) ]

    # If jobs are incomplete, print diagnostic and exit
    if len(incomplete) > 0:
        print("\nERROR in md_runs: the following initializations didn't complete:\n\n\t{}".format("\n\t".join(incomplete)))
        quit()

    # If jobs didn't condense print diagnostic and exit
    if len(diffuse) > 0:
        print("WARNING in md_runs: the following initializations appear to be diffuse:\n\n\t{}".format("\n\t".join(diffuse)))
        quit()

    # Run condensed phase extensions
    N = int(np.ceil(float(c["fom_t_tot"])/float(c["fom_t_ext"])))
    restarts = { _:0 for _ in c["xyz"] }

    print("\trunning {} condensed phase extensions...".format(N))
    break_flag = False
    for z in range(N):
        jobids = []

        # Early escape if combined trajectory is discovered
        if break_flag is True:
            break

        # The job submission is nested inside of a restart loop
        for r in range(c["fom_restarts"]+1):
            for i in c["xyz"]:

                # Save base name and versions to variables
                name=".".join(i.split('.')[:-1])
                name_sub=name.replace('(','-')
                name_sub=name_sub.replace(')','-')
                name_safe=esc_str(name)

                # Get the simulation temperatures
                T = get_temps(c["fom_temps"],name)

                # Loop over temperatures and trajectories
                for t in T:
                    for j in range(1,c["fom_n_traj"]+1):

                        N_run = get_run_N("FOM/{}_cond/{}_K/run_{}/extend.in.init".format(name,t,j))

                        # If the extension doesn't exist then submit it                    
                        if N_run == z:

                            # Print warning if this is a resubmission
                            if os.path.isfile("FOM/{}_cond/{}_K/run_{}/{}.sys.lammpstrj".format(name,t,j,z)) is True:
                                print("WARNING: job run {} for FOM/{}_cond/{}_K/run_{}/extend.in.init is incomplete, attempting resubmission".format(z,name,t,j))

                            # Submit job
                            substring="python {}/Automation_Scripts/lammps_submit.py -f extend.in.init -d {}_cond^{}_K^run_{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} -path_to_exe {} -o {} --overwrite --silent"
                            substring = substring.format(c["taffi_path"],name,t,j,c["fom_md_cond_procs"],c["fom_md_cond_ppn"],c["fom_md_cond_size"],c["fom_md_cond_q"],c["fom_md_cond_sched"],c["fom_md_cond_wt"],c["acct"],\
                                                         c["lammps_exe"],"cond_ext_{}_{}_{}".format(name_sub,t,j))
                            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                            output = str(output,'utf-8')
                            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

                        # Check for combined trajectory
                        if os.path.isfile("FOM/{}_cond/{}_K/run_{}/combined.lammpstrj".format(name,t,j)) is True:
                            break_flag = True
                                          
            monitor_jobs(jobids,'lin1209')

    # Combine the output trajectories
    print("\tcombining trajectories...")
    N_start = int(np.ceil(float(c["fom_t_equil"])/float(c["fom_t_ext"])))   # N_start is based upon how long the user wants a pre-equilibration
    N_end = N - 1
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                # Combine the production trajectories 
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/combined.lammpstrj".format(name,t,j)) is False:
                    os.chdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j))
                    substring="{}/Automation_Scripts/combine_ext.sh *.sys.lammpstrj {} {}".format(c["taffi_path"],N_start,N_end)
                    output = subprocess.check_output(subproc_string(substring), encoding='utf8').strip("\r\n")
                    os.chdir(current_dir)

                # Combine the production thermo files
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/combined.thermo".format(name,t,j)) is False:
                    os.chdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j))
                    substring="{}/Automation_Scripts/combine_ext.sh *.thermo.avg {} {} --thermo".format(c["taffi_path"],N_start,N_end)
                    output = subprocess.check_output(subproc_string(substring), encoding='utf8').strip("\r\n")
                    os.chdir(current_dir)

    # Rerun trajectories to obtain unaveraged thermodynamic properties
    print("\trerunning combined trajectory to generate unaveraged thermodynamic data...")
    jobids=[]    
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):

                # # TMP REMOVE FOR PRODUCTION
                # if os.path.isfile("FOM/{}_cond/{}_K/run_{}/rerun.in.init".format(name,t,j)):
                #     os.remove("FOM/{}_cond/{}_K/run_{}/rerun.in.init".format(name,t,j))
                # if os.path.isfile("FOM/{}_cond/{}_K/run_{}/rerun.in.out".format(name,t,j)):
                #     os.remove("FOM/{}_cond/{}_K/run_{}/rerun.in.out".format(name,t,j))
                # if os.path.isfile("FOM/{}_cond/{}_K/run_{}/combined.thermo.rerun".format(name,t,j)):
                #     os.remove("FOM/{}_cond/{}_K/run_{}/combined.thermo.rerun".format(name,t,j))


                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/rerun.in.init".format(name,t,j)) is False:
                    write_rerun("FOM/{}_cond/{}_K/run_{}/extend.in.init".format(name,t,j),"FOM/{}_cond/{}_K/run_{}/rerun.in.init".format(name,t,j))

                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/combined.thermo.rerun".format(name,t,j)) is False:

                    # Submit job
                    substring="python {}/Automation_Scripts/lammps_submit.py -f rerun.in.init -d {}_cond^{}_K^run_{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} -path_to_exe {} -o {} --overwrite --silent"
                    substring = substring.format(c["taffi_path"],name,t,j,1,c["fom_md_cond_ppn"],c["fom_md_cond_size"],c["fom_md_cond_q"],c["fom_md_cond_sched"],c["fom_md_cond_wt"],c["acct"],\
                                                 c["lammps_exe"],"rerun_{}_{}_{}".format(name_sub,t,j))
                    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                    output = str(output,'utf-8')
                    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

    monitor_jobs(jobids,'lin1209')

    # Parse the diffusion coefficients and dielectric constants
    print("\tparsing diffusion and dielectric constant from trajectories...")
    bundle_string = ""
    bundle_count = 0
    b_count = 0
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','_')
        name_sub=name_sub.replace(')','_')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                # Check completeness, remove incomplete folders 
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/msd/avg-0.txt".format(name,t,j)) is False:
                    shutil.rmtree("FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j))

                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/eps/0_ac.txt".format(name,t,j)) is False:
                    shutil.rmtree("FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j))

                # Parse diffusivities
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)) is False:

                    # Assemble the submission string for the batch submission file
                    substring="python {}/Parsers/msd_parse_new.py -map run_{}.map -traj combined.lammpstrj -mols 0 -o msd -f_every 10 &\n"
                    substring=substring.format(c["taffi_path"],j)
                    bundle_string+="\ncd {}/FOM/{}_cond/{}_K/run_{}\n".format(os.getcwd(),name_safe,t,j)+substring

                    # Assemble the submission string for the batch submission file
                    b_count += 1
                    if b_count == int(c["fom_geoopt_ppn"]):
                        bundle_string+="\nwait\n"
                        #shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
                        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"], "standby",c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd())) ## FOR SLURM
                        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
                        output = output.split()[-1]
                        jobids += [output]
                        bundle_string = ""
                        bundle_count += 1
                        b_count = 0

                # Parse dielectric constants
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)) is False:

                    # Assemble the submission string for the batch submission file
                    #substring="python {}/Parsers/dielectric_parse_new.py -map run_{}.map -traj combined.lammpstrj -every 1 -T {} -mols 0 -folder eps &\n"
                    substring="python {}/Parsers/dielectric_parse.new.py -map run_{}.map -traj combined.lammpstrj -every 10 -T {} -mols 0 -folder eps &\n" ## 10 ps
                    substring=substring.format(c["taffi_path"],j,t)
                    bundle_string+="\ncd {}/FOM/{}_cond/{}_K/run_{}\n".format(os.getcwd(),name_safe,t,j)+substring

                    # Assemble the submission string for the batch submission file
                    b_count += 1
                    if b_count == int(c["fom_geoopt_ppn"]):
                        bundle_string+="\nwait\n"
                        #shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
                        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"], "standby",c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd())) 
                        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
                        output = output.split()[-1]
                        jobids += [output]
                        bundle_string = ""
                        bundle_count += 1
                        b_count = 0
    
    # submit remainder of jobs
    if b_count != 0:
        bundle_string+="\nwait\n"
        #shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"], "standby",c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd())) ## FOR SLURM
        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
        output = output.split()[-1]
        jobids += [output]
        bundle_string = ""
        bundle_count += 1
        b_count = 0

    # Wait until the parses complete
    monitor_jobs(jobids,'lin1209')

    # Check parse for completion
    i_msd = []
    i_eps = []
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                # Check completeness, remove incomplete folders 
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/msd/avg-0.txt".format(name,t,j)) is False:
                    i_msd += ["FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)]

                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/eps/0_ac.txt".format(name,t,j)) is False:
                    i_eps += ["FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)]

    # Exit if there are incompletes
    if len(i_msd) > 0:
        print("\nThe following msd parses failed:\n\n\t{}".format("\n\t".join([ str(_) for _ in i_msd ])))
    if len(i_eps) > 0:
        print("\nThe following eps parses failed:\n\n\t{}".format("\n\t".join([ str(_) for _ in i_eps ])))        
    if len(i_msd) > 0 or len(i_eps) > 0:
        quit()

    # Clean up submission files
    for f in glob.glob("*.submit*"):
        os.remove(f)
    for f in glob.glob("*.err"):
        os.remove(f)
    for f in glob.glob("*.out"):
        os.remove(f)

    # Parse thermodynamic figures of merit
    print("\tparsing thermodynamic averages...")
    kb=0.001987204118                           # Boltzman constant in units of kcal/mol K
    kcal_mol_to_J = 4.184 * 1000.0 / 6.0221409E23 # employed for converting the isothermal compressibility to useful units.  
    kcal_to_kJ = 4.184
    pressure = 1.0                    # for thermal expansion coeff (atm)
    atm_to_pa = 101325.0              # atmospheres to pascals
    PV_const = pressure * 1E-30 * atm_to_pa / kcal_mol_to_J # A^-3 to m^-3, conversion to pa, conversion to kcal/mol
    Data = {}
    Data_avg = {}
    if os.path.isdir("FOM/Results") is False:
        os.mkdir("FOM/Results")

    # Create header of thermo file
    if os.path.isfile("FOM/Results/thermo.txt") is False:
        with open("FOM/Results/thermo.txt",'w') as f:
            f.write("{:<40s} {:<10s} {:<16s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s}\n"\
                    .format("compound","T(k)","d(g/cm^3)","std(g/cm^3)","H_val(kJ/mol)","std(kJ/mol)","eps","std","iso(GP^-1)","std(GP^-1)","t_exp(10^-3/K)","std(10^-3/K)"))
            write_flag = True
    else:
        write_flag = False

    # For each molecule that was simulated, parse the density, dielectric, diffusivity, thermal exp coeff, and bulk modulus
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        Data[name] = {}
        Data_avg[name] = {}

        # Grab the number of atoms in this molecule
        with open(i,'r') as f:
            for lc,lines in enumerate(f):
                N_atoms = int(lines.split()[0])
                break

        # Calculate the number of molecules in the simulation
        N_mol = int(float(5000)/float(N_atoms))
        N_atoms = N_mol * N_atoms

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            Data[name][t] = {}
            for j in range(1,c["fom_n_traj"]+1):

                # Read data into dictionary
                Data[name][t][j] = {}
                Data[name][t][j]["d"],Data[name][t][j]["U"],Data[name][t][j]["V"],Data[name][t][j]["H"] = get_cols("FOM/{}_cond/{}_K/run_{}/combined.thermo.rerun".format(name,t,j),[2,4,3,8],2)
                Data[name][t][j]["U_gas"] = get_cols("FOM/{}_gas/{}_K/run_{}/thermo.equil.avg".format(name,t,j),[4],2)
                Data[name][t][j]["t_msd"],Data[name][t][j]["msd"] = get_cols("FOM/{}_cond/{}_K/run_{}/msd/avg-0.txt".format(name,t,j),[0,1],1)
                Data[name][t][j]["t_eps"],Data[name][t][j]["eps_ravg"] = get_cols("FOM/{}_cond/{}_K/run_{}/eps/0_eps.txt".format(name,t,j),[0,1],1) 

                # Convert to arrays
                Data[name][t][j]["d"] = np.array(Data[name][t][j]["d"])
                Data[name][t][j]["U"] = np.array(Data[name][t][j]["U"])
                Data[name][t][j]["V"] = np.array(Data[name][t][j]["V"])
                Data[name][t][j]["H"] = np.array(Data[name][t][j]["H"])
                Data[name][t][j]["U_gas"] = np.array(Data[name][t][j]["U_gas"])
                Data[name][t][j]["t_msd"] = np.array(Data[name][t][j]["t_msd"])
                Data[name][t][j]["msd"] = np.array(Data[name][t][j]["msd"])
                Data[name][t][j]["t_eps"] = np.array(Data[name][t][j]["t_eps"])
                Data[name][t][j]["eps_ravg"] = np.array(Data[name][t][j]["eps_ravg"])

                # Calculate isothermal compressibility
                Data[name][t][j]["iso_comp"] = np.std(Data[name][t][j]["V"])**(2.0)/(np.mean(Data[name][t][j]["V"])*kb*float(t))*1E-30*1E9/kcal_mol_to_J 

                # # Calculate thermal expansion coefficient (OLD)
                # tmp_dev_V = Data[name][t][j]["V"] - mean(Data[name][t][j]["V"])
                # tmp_dev_H = Data[name][t][j]["H"] - mean(Data[name][t][j]["H"])
                # Data[name][t][j]["therm_exp"] = mean(tmp_dev_V*tmp_dev_H)/(kb*float(t)**(2.0)*mean(Data[name][t][j]["V"]))*1E3

                # Calculate thermal expansion coefficient
                Data[name][t][j]["therm_exp"] = ( np.mean(Data[name][t][j]["V"]*Data[name][t][j]["U"]) - np.mean(Data[name][t][j]["V"])*np.mean(Data[name][t][j]["U"]) + PV_const * np.std(Data[name][t][j]["V"])**(2.0) ) \
                                         / (kb*float(t)**(2.0)*np.mean(Data[name][t][j]["V"])) * 1E3

                # Read eps from file
                with open("FOM/{}_cond/{}_K/run_{}/eps/dielectric_parse.log".format(name,t,j),'r') as f:
                    for lines in f:
                        fields = lines.split()
                        if len(fields) == 2 and fields[0] == "eps:":
                            Data[name][t][j]["eps"] = float(fields[1])

                print("{} {} {}: {} {} {}".format(name,t,j,Data[name][t][j]["iso_comp"],Data[name][t][j]["therm_exp"],Data[name][t][j]["eps"]))

                # Save thermodynamic mean and std for each run to a dedicated mean_thermo.txt file
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j)) is True:
                    os.remove("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j))
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j)) is False:
                    with open("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j),'w') as f:
                        f.write("{:<20s}  {:<20s} {:<20s}\n".format("property","mean","stdev"))
                        f.write("{:<20s} {:< 20.8f} {:< 20.8f}\n".format("density",np.mean(Data[name][t][j]["d"]),np.std(Data[name][t][j]["d"])))
                        f.write("{:<20s} {:< 20.8f} {:< 20.8f}\n".format("U",np.mean(Data[name][t][j]["U"]),np.std(Data[name][t][j]["U"])))

            # Make plot of average internal energy vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_U_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_U_composite.png".format(name,t),\
                          [ np.arange(len(Data[name][t][j]["U"]))*0.001 for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["U"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="U (kcal/mol)")                      

            # Make plot of density vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_d_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_d_composite.png".format(name,t),\
                          [ np.arange(len(Data[name][t][j]["d"]))*0.001 for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["d"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="d (g/cm^3)")                      

            # Make plot of msd vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_msd_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_msd_composite.png".format(name,t),\
                          [ Data[name][t][j]["t_msd"] for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["msd"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="MSD (A^2)")                      
            
            # Make plot of running average of eps vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_eps_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_eps_composite.png".format(name,t),\
                          [ Data[name][t][j]["t_eps"] for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["eps_ravg"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="eps")                      

        # Save averages of important quatities
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                Data[name][t][j]["d"] = np.mean(Data[name][t][j]["d"])
                Data[name][t][j]["H_vap"] = ( np.mean(Data[name][t][j]["U_gas"]) - np.mean(Data[name][t][j]["U"])/float(N_mol) + kb*float(t) ) * kcal_to_kJ


        # Calculate averages over the trajectories
        keys = ["d","H_vap","iso_comp","therm_exp","eps","t_msd","msd"]
        for t in T:
            Data_avg[name][t] = {}
            for j in range(1,c["fom_n_traj"]+1):
                for k in keys:
                    if j == 1:
                        Data_avg[name][t][k] = Data[name][t][j][k]
                        Data_avg[name][t][k+"_std"] = Data[name][t][j][k]**(2.0)
                    else:
                        Data_avg[name][t][k] = Data_avg[name][t][k] + Data[name][t][j][k]
                        Data_avg[name][t][k+"_std"] = Data_avg[name][t][k+"_std"] + Data[name][t][j][k]**(2.0)
            
            # Calculate average and std_dev
            for k in keys:
                Data_avg[name][t][k] = Data_avg[name][t][k] / float(c["fom_n_traj"])
                #Data_avg[name][t][k+"_std"] = ( Data_avg[name][t][k+"_std"] / float(c["fom_n_traj"]) - Data_avg[name][t][k]**(2.0) )**(0.5)
                Data_avg[name][t][k+"_std"] = ( np.abs(Data_avg[name][t][k+"_std"] / float(c["fom_n_traj"]) - Data_avg[name][t][k]**(2.0)) )**(0.5) ## NEED TO BE VALIDATED

            # Check convergence with respect to energy
            for j in range(1,c["fom_n_traj"]+1):
                if j == 1:
                    Data_avg[name][t]["U"] = Data[name][t][j]["U"]
                else:
                    Data_avg[name][t]["U"] = Data_avg[name][t]["U"] + Data[name][t][j]["U"]
            Data_avg[name][t]["U"] = Data_avg[name][t]["U"] / float(c["fom_n_traj"])
            if abs((np.mean(Data_avg[name][t]["U"][:1000]) - np.mean(Data_avg[name][t]["U"][-1000:]))/float(N_atoms)) > 0.01:
                print("{} {} didn't converge ({})".format(name,t,abs((np.mean(Data_avg[name][t]["U"][:1000]) - np.mean(Data_avg[name][t]["U"][-1000:]))/float(N_atoms))))

            if write_flag is True:
                with open("FOM/Results/thermo.txt",'a') as f:                
                    f.write("{:<40s} {:<10s} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f}\n"\
                            .format(name,t,Data_avg[name][t]["d"],Data_avg[name][t]["d_std"],Data_avg[name][t]["H_vap"],Data_avg[name][t]["H_vap_std"],\
                                    Data_avg[name][t]["eps"],Data_avg[name][t]["eps_std"],Data_avg[name][t]["iso_comp"],Data_avg[name][t]["iso_comp_std"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["therm_exp_std"]))

        # Calculate the derivatives
        for t in T:
            Data_avg[name][t]["msd_slope_lin"],Data_avg[name][t]["msd_std_lin"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"])
            Data_avg[name][t]["msd_slope_log"],Data_avg[name][t]["msd_std_log"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"],logx_opt=True,logy_opt=True)
            Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"] = calc_mobility(np.array(Data_avg[name][t]["msd_slope_log"]),np.array(Data_avg[name][t]["msd_slope_lin"]),start=100) 
#            print "{} {}: {} {} {} {} {}".format(name,t,Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["iso_comp"])
    quit()

# Write lammps rerun.in.init file for rerunning the combined trajectory without averaging
def write_rerun(input_file,output_file):
    
    with open(input_file,'r') as f:
        with open(output_file,'w') as o:
            for lines in f:
                fields = lines.split()
                # reset run index, necessary for backwards compatibility, since lammps doesn't like to mix rerun with read_restart
                if len(fields) > 3 and fields[0] == "variable" and fields[1] == "run" and fields[2] == "index":
                    o.write("variable    run             index   0\n")
                elif len(fields) >= 2 and fields[0] == "read_restart":
                    o.write("read_data ${data_name}\n")
                elif len(fields) > 3 and fields[0] == "fix" and fields[1] == "averages" and fields[2] == "all" and fields[3] == "ave/time":
                    o.write("fix raw all ave/time 1 1 ${thermo_freq} v_my_temp v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file combined.thermo.rerun format %20.10g\n\n"+\
                            "# Rerun trajectory\n"+\
                            "rerun combined.lammpstrj dump x y z\n")
                    break
                else:
                    o.write(lines)
    return

# Description: Calculate the mobility based on the first continuous window of N datapoints with a loglog slope greater than a_thresh 
# 
# inputs: y_log (array): slope from the loglog derivative data
#         y_lin (array): slope from the linlin derivative data
#         a_thresh (float): slope used for evaluating the diffusive regime
#         N_thresh (int): number of continuous frames with y_log slope greater than a_thresh used for calculating the derivative.
#         start (int): index to start the search at
def calc_mobility(y_log,y_lin,a_thresh=0.95,N_thresh=10,start=0):

    count = 0
    inds = []
    for i in range(len(y_log)):
        if y_log[i] > a_thresh and i > start:
            count += 1
            inds += [i]
            if count >= N_thresh:
                return np.mean(y_lin[inds])/6.0*1.0E-16*1.0E9,np.std(y_lin[inds]/6.0*1.0E-16*1.0E9),inds
        else:
            count = 0
            inds = []

# Description: Calculates the derivative of a dataset using linear least-squares regression on a running window of the data.
#
# inputs: x (array): x data for calculating the derivative
#         y (array): y data for calculating the derivative
#         abs_flag (boolean): option for taking the absolute value of the data before calculating the derivative
#         logx_opt (boolean): option for taking the log of the x data before calculating the derivative
#         logy_opt (boolean): option for taking the log of the y data before calculating the derivative
#         x_scale (float): scale factor applied to the x data before calculating the derivative
#         y_scale (float): scale factor applied to the y data before calculating the derivative
#         mode (string): forward or central windowing used for calculating the derivative
#         window (int): number of datapoints to use in the running window for calculating the derivative. In central mode, the 
#                       regression is applied to +/- this number of datapoints about the point that the derivative is being evaluated.
#                       In forward mode, + this number of datapoints about the point that the derivative is being evaluated.
def calc_derivative(x,y,abs_flag=False,logx_opt=False,logy_opt=False,x_scale=1.0,y_scale=1.0,mode="forward",window=10):

    valid_modes = ["forward","central"]
    if mode not in valid_modes:
        print("ERROR in calc_derivative: mode must be set to one of the following {}".format(valid_modes))
        quit()

    # Apply user-defined scaling options
    x_0 = deepcopy(x)
    y_0 = deepcopy(y)
    if abs_flag:
        x_0 = np.abs(x_0)
        y_0 = np.abs(y_0)
    if logx_opt:
        for ind, x_tmp in enumerate(x_0): # ADDED
            if x_tmp == 0:                # ADDED
                x_0[ind] = 1E-16          # ADDED
        x_0 = np.log10(x_0)
    if logy_opt:
        for ind, y_tmp in enumerate(y_0): # ADDED
            if y_tmp == 0:                # ADDED
                y_0[ind] = 1E-16          # ADDED
        y_0 = np.log10(y_0)
 
    # Apply scaling
    y_0 = y_0*y_scale
    x_0 = x_0*x_scale
    N = len(x_0)
    
    slopes = []
    stds = []

    # Calculate the derivative
    if mode == "central":
        for i in range(len(x_0)):

            # starting indices can only do forward windows
            if i < window:
                inds = list(range(0,i)) + list(range(i,i+window+1))

            # ending indices can only do backwards windows
            elif i+window >= N:
                inds = list(range(i-window,i)) + list(range(i,N))

            # the rest can sample the full window
            else:
                inds = list(range(i-window,i+window+1))

            slope, intercept, r_value, p_value, std_err = linregress(x_0[inds],y_0[inds])
            slopes += [slope]
            stds += [std_err]

    # Calculate the derivative
    if mode == "forward":
        for i in range(len(x_0)):

            # ending indices can only do fractional windows
            if i+window >= N:
                inds = list(range(i,N))

            # the rest can sample the full window (N-2 because derivative must involve at least 2 points
            elif i < N-2:
                inds = list(range(i,i+window+1))

            # calculate derivative
            if i < N-2:
                slope, intercept, r_value, p_value, std_err = linregress(x_0[inds],y_0[inds])

                slopes += [slope]
                stds += [std_err]

    return slopes,stds

# # Print FOM summary the densities and print summary of values
# kb=0.001987204118
# T=298
# kcal2kJ=4.184
# printf "\n%40s   %20s   %40s %20s" "Molecule" "density (g/cm^3)" "diff (1E5 * cm^2/s)" "H_vap (kcal/mol)" > summary.txt
# for i in *cond/; do
#     name=$( echo ${i} | awk -F '_cond' '{ print $1 }' )
#     density=$(${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 60000 3)
#     diff=$( awk '{ if ( $1 == "The" && $2 == "Self-Diffusion" && $3 == "Coefficient" ) { printf("%12.6f",$9*1.0E5) }}' ${name}_cond/diff/msd_parse.log)
#     U_vap=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_gas/thermo.equil.avg 1000 0 5 )
#     U_bulk=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 0 5 )
#     N_mol=$( awk '{ if ( NR == 1 ) { print $1; exit } }' ${name}.xyz )
#     N_mol=$( echo "1000/${N_mol}" | bc )
#     U_bulk=$( echo "${U_bulk} / ${N_mol}" | bc -l )
#     H_vap=$( echo "(${U_vap} - ${U_bulk} + ${kb} * ${T} ) * ${kcal2kJ}" | bc -l )
#     printf "\n%40s %12.4f %40s %23.4f" ${name} ${density} ${diff} ${H_vap} >> summary.txt
# done
# echo "" >> summary.txt

    quit()

# Wrapper function for parsing the run value from a lammps file. Used to determine job completion.
def get_run_N(input_file):
    with open(input_file,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) >= 4 and fields[0] == "variable" and fields[1] == "run" and fields[2] == "index":
                return int(fields[3])
    return None

# Wrapper function for parsing the list of simulation temperatures from the temperature file
# input: t_file: text file holding list of molecule names and temperatures
#        name:  molecule name to return the temperatures for
def get_temps(t_file,name):
    if t_file is None:
        T = ["298.00"]
    else:
        T = None
        with open(t_file,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0:
                    if fields[0] == name:
                        T = [ str(_) for _ in fields[1:] ]
                        break
        if T is None:
            T = ["298.00"]
    return T

# Function for writing shell based input file
#def shell_submit(cmd,output,p,t,q,ppn=0,a=None,run_dir=None,batch="pbs"):
def shell_submit(cmd, output, p, t, a=None, ppn=0, batch="slurm", run_dir=None):
    N_nodes = int(np.ceil(float(p)/float(ppn)))
    if "min" in str(t):
        t = t.strip("min")
        min_flag = True
    else:
        min_flag = False
    if run_dir is None:
        run_dir = os.getcwd()
    with open(output,'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#\n")

        if batch == "pbs":
            f.write("#PBS -N {} \n".format(output))
            if a is not None:
                f.write("#PBS -A {}\n".format(a))
            if ppn == 0:
                f.write("#PBS -l nodes={}\n".format(N_nodes))
            else:
                f.write("#PBS -l nodes={}:ppn={}\n".format(N_nodes,ppn))
            if min_flag is True:
                f.write("#PBS -l walltime=00:{}:00\n".format(t))
            else:
                f.write("#PBS -l walltime={}:00:00\n".format(t))
            f.write("#PBS -q {}\n".format(q))
            f.write("#PBS -S /bin/sh\n")

        if batch == "slurm":
            f.write("#SBATCH --job-name={}\n".format(output))
            f.write("#SBATCH --output={}out \n".format(output[:-6]))
            f.write("#SBATCH --error={}err \n".format(output[:-6]))
            if a is not None:
                f.write("#SBATCH -A {}\n".format(a))
            if ppn == 0:
                f.write("#SBATCH --nodes={}\n".format(N_nodes))
            else:
                f.write("#SBATCH --nodes={}\n".format(N_nodes))
                #f.write("#SBATCH --tasks-per-node={}\n".format(ppn))
                f.write("#SBATCH --ntasks-per-node={}\n".format(ppn))
            #if min_flag is True:
            #    f.write("#SBATCH --time 00:{}:()\n".format(t))
            #else:
            f.write("#SBATCH --time 4:00:00\n")
            f.write("\n# cd into the submission directory\n")
            #f.write("cd -- $SLURM_SUBMIT_DIR\n")
            f.write("cd -- {}\n".format(run_dir))
            #f.write("echo Working directory is $SLURM_SUBMIT_DIR\n")
            f.write("echo Working directory is {}\n".format(run_dir))
            f.write("echo Running on host $SLURM_SUBMIT_HOST\n")

        #f.write("\n# cd into the submission directory\n")
        #f.write("cd -- {}\n".format(run_dir))
        #f.write("echo Working directory is {}\n".format(run_dir))
        #f.write("echo Running on host `hostname`\n")
        #f.write("echo Time is `date`\n\n")
        f.write('module load gcc/9.3.0 \n')
        f.write('module load openmpi/3.1.4 \n')
        f.write('module load ffmpeg/4.2.2 \n')
        f.write('module load  openblas/0.3.8 \n')
        f.write('module load  gsl/2.4 \n\n')

        f.write("# Copy the local path\n")
        f.write("PATH={}\n".format(os.getenv("PATH")))
        f.write("sleep 2\n\n# Run Commands\n{}\nsleep 2\n".format(cmd))

    return
                

# Inserts escape characters before parentheses so that the script can                                                                                                                                                                         
# handle filenames/folders that involve these                                                                                                                                                                                                 
def esc_str(a):
    a = a.replace('(','\\(')
    a = a.replace(')','\\)')
    return a

def num_frames(traj):
    N_frames = 0
    with open(traj,'r') as f:
        for lines in f:
            if "ITEM: TIMESTEP" in lines:
                N_frames += 1
    return N_frames

def run_dilatometry(c):

    # Make dilatometry directory
    if os.path.isdir("dilatometry") is False:
        os.mkdir("dilatometry")

    # Loop over pressures and create inputs if they don't exist
    print("\tgenerating dilatometry inputs...")
    for i in c["pressures"]:

        # Make the pressure folder if it doesn't exist
        if os.path.isdir("dilatometry/{}".format(i)) is False:
            os.mkdir("dilatometry/{}".format(i))

        # Make the runs if they don't exist
        for j in range(1,int(c["n_dil_runs"])+1):

            if os.path.isdir("dilatometry/{}/run_{}".format(i,j)) is False:
                substring="python {}/FF_functions/gen_md_for_sampling.py {} -N {} -gens {} -FF {} -pair_styles {} "+\
                          "-onefourscale_coul {} -onefourscale_lj {} -t {} -t_A {} -t_ext {} -T {} -T_A {} -T_D {} "+\
                          "-N_dil {} --tail --impropers --velocities -P {} -o {} --force_read -f {}"
                substring = substring.format(c["taffi_dir"],c["xyz"],c["num_mol"],c["gens"],c["ff"],','.join(c["pair_styles"].split()),c["14coul"],c["14lj"],c["t_dil_equil"],c["t_dil_anneal"],c["t_dil_ext"],\
                                             c["temp_dil_start"],c["temp_dil_start"],c["temp_dil_end"],c["dil_intervals"],i,"dilatometry/{}/run_{}".format(i,j),c["print_freq"])
                subprocess.check_output(subproc_string(substring), encoding='utf8')

    # Run initializations
    jobids = []
    print("\trunning dilatometry initializations...")
    for i in c["pressures"]:
        for j in range(1,int(c["n_dil_runs"])+1):
            if os.path.isfile("dilatometry/{}/run_{}/extend.end.restart".format(i,j)) is False:
                substring="python {}/Automation_Scripts/lammps_submit.py -f {}.in.init -d {},{},{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} --overwrite -path_to_exe {} -o {} -shell {} --silent"
                substring = substring.format(c["taffi_dir"],"run_{}".format(j),"dilatometry",i,"run_{}".format(j),c["procs"],c["ppn"],c["size"],c["queue"],c["sched"],c["wt"],\
                                             c["acct"],c["lammps_exe"],"{}_run_{}".format(i,j),','.join(c["shell_cmd"].split()))
                output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                output = str(output,'utf-8')
                jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

    monitor_jobs(jobids,'lin1209')

    # Run extensions
    for t in range(int(c["dil_intervals"])):
        
        jobids = []
        print("\trunning dilatometry interval {}...".format(t))
        for i in c["pressures"]:
            for j in range(1,int(c["n_dil_runs"])+1):

                # If the output trajectory for this interval is not present then attempt to run the job
                if os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)) is False:

                    # If the prerequisite restart file is missing, print diagnostic to user and exit
                    if os.path.isfile("dilatometry/{}/run_{}/extend.end.restart".format(i,j)) is False:
                        print("ERROR in dilatometry phase: the initialization run ({}) did not complete.".format("dilatometry/{}/run_{}".format(i,j)))
                        quit()

                    # If the trajectory from the previous interval is missing, print diagnostic to user and exit.
                    elif t > 0 and os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t-1)) is False:
                        print("ERROR in dilatometry phase: the extension run ({}) did not complete.".format("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)))
                        quit()

                    # Else, everything is in place and the job is run.
                    else:
                        substring="python {}/Automation_Scripts/lammps_submit.py -f {}.in.init -d {},{},{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} --overwrite -path_to_exe {} -o {} -shell {} --silent"
                        substring = substring.format(c["taffi_dir"],"extend","dilatometry",i,"run_{}".format(j),c["procs"],c["ppn"],c["size"],c["queue"],c["sched"],c["wt"],\
                                                     c["acct"],c["lammps_exe"],"{}_run_{}".format(i,j),','.join(c["shell_cmd"].split()))
                        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
                        output = str(output,'utf-8')
                        jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

        monitor_jobs(jobids,'lin1209')
        
        # Check that all of the runs completed
        for i in c["pressures"]:
            for j in range(1,int(c["n_dil_runs"])):

                # If the trajectory from the previous interval is missing, print diagnostic to user and exit.
                if os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)) is False:
                    print("ERROR in dilatometry phase: the extension run ({}) did not complete.".format("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)))
                    quit()

    # Combine trajectories
    print("\tcombining dilatometry trajectories...")
    current_dir = os.getcwd()
    for i in c["pressures"]:
        for j in range(1,int(c["n_dil_runs"])+1):

            # Remove incomplete combinations if present (e.g., if the use aborted the script execution in this step)
            if os.path.isfile("dilatometry/{}/run_{}/combined.lammpstrj".format(i,j)):
                os.remove("dilatometry/{}/run_{}/combined.lammpstrj".format(i,j))
            if os.path.isfile("dilatometry/{}/run_{}/combined.thermo".format(i,j)):
                os.remove("dilatometry/{}/run_{}/combined.thermo".format(i,j))

            # Combine the trajectory and thermo files if they haven't been combined
            if os.path.isfile("dilatometry/{}/run_{}/{}_run_{}_combined.lammpstrj".format(i,j,i,j)) is False:

                # Run combine script
                os.chdir("dilatometry/{}/run_{}".format(i,j))
                substring="{}/Automation_Scripts/combine_ext.sh *.sys.lammpstrj 0 {}".format(c["taffi_dir"],int(c["dil_intervals"])-1)
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                substring="{}/Automation_Scripts/combine_ext.sh *.thermo.avg 0 {} --thermo".format(c["taffi_dir"],int(c["dil_intervals"])-1)
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                shutil.move("combined.lammpstrj","{}_run_{}_combined.lammpstrj".format(i,j))
                shutil.move("combined.thermo","{}_run_{}_combined.thermo".format(i,j))
                os.chdir(current_dir)

    # Make run folders if they don't exist 
    if os.path.isdir("Results") is False:
        os.mkdir("Results")
    if os.path.isdir("Results/dilatometry") is False:
        os.mkdir("Results/dilatometry")    

    # Process dilatometry data and calculate approximate Tg values
    data={}
    Tg={}
    print("\tprocessing dilatometry data...")
    for i in c["pressures"]:
        data[i] = {}
        for j in range(1,int(c["n_dil_runs"])+1):            
            data[i][j] = { "T":[],"d":[],"U":[],"H":[] }
            if os.path.isfile("Results/dilatometry/{}_run_{}_smooth.thermo".format(i,j)) is False:

                # Smooth the data
                substring="python {}/smooth_data.py -f {} -o Results/dilatometry/{}_run_{}_smooth.thermo -start 0 -end {} -head 2 -x_col 1 -y_col 2,4,8 -x_min {} -x_max {} -x_step 10 -w 100 "+\
                          "-labels T(K),density(g/cm^3),PE(kcal/mol),H(kcal/mol) -y_scale 1.0,{},{}"
                substring = substring.format(c["colops_dir"],"dilatometry/{}/run_{}/{}_run_{}_combined.thermo".format(i,j,i,j),i,j,int(float(c["t_dil_ext"])*float(c["dil_intervals"])/float(c["print_freq"])),\
                                             c["temp_dil_end"],c["temp_dil_start"],1.0/float(c["num_mol"]),1.0/float(c["num_mol"]))
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')

            # Collect the smoothed data
            data[i][j]["T"],data[i][j]["d"],data[i][j]["U"],data[i][j]["H"] = get_cols("Results/dilatometry/{}_run_{}_smooth.thermo".format(i,j),[0,1,2,3],1)

            # Make the Tvd TvU and TvH plots (if necessary)
            T_raw,d_raw,U_raw,H_raw = get_cols("dilatometry/{}/run_{}/{}_run_{}_combined.thermo".format(i,j,i,j),[1,2,4,8],2)
            U_raw = np.array(U_raw)/float(c["num_mol"])
            H_raw = np.array(H_raw)/float(c["num_mol"])

            if os.path.isfile("Results/dilatometry/{}_run_{}_Tvd_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_Tvd_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[d_raw,data[i][j]["d"]],x_label="T (K)",y_label="d (g/cm^3)")            
            if os.path.isfile("Results/dilatometry/{}_run_{}_TvU_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_TvU_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[U_raw,data[i][j]["U"]],x_label="T (K)",y_label="U (kcal/mol)")
            if os.path.isfile("Results/dilatometry/{}_run_{}_TvH_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_TvH_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[H_raw,data[i][j]["H"]],x_label="T (K)",y_label="H (kcal/mol)")

        # Make plots with all runs in the same figure
        if os.path.isfile("Results/dilatometry/{}_Tvd_runs.pdf".format(i)) is False:        
            make_plot("Results/dilatometry/{}_Tvd_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["d"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="d (g/cm^3)")                    
        if os.path.isfile("Results/dilatometry/{}_TvU_runs.pdf".format(i)) is False:
            make_plot("Results/dilatometry/{}_TvU_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["U"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="U (kcal/mol)")
        if os.path.isfile("Results/dilatometry/{}_TvH_runs.pdf".format(i)) is False:
            make_plot("Results/dilatometry/{}_TvH_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["H"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="H (kcal/mol)")

        # Calculate averages and stdevs across the runs            
        for k in ["d","U","H"]:
            for j in range(1,int(c["n_dil_runs"])+1):

                # Initialize arrays if this is the first being read
                if j == 1:
                    x_mean = np.array(data[i][j]["T"])
                    y_mean = np.array(data[i][j][k])
                    y_std  = np.array(data[i][j][k])**(2.0)

                # Else, add to arrays
                else:
                    x_mean += np.array(data[i][j]["T"])
                    y_mean += np.array(data[i][j][k])
                    y_std  += np.array(data[i][j][k])**(2.0)

            # Calculate quantities
            x_mean = x_mean / float(int(c["n_dil_runs"]))
            y_mean = y_mean / float(int(c["n_dil_runs"]))
            y_std  = ( y_std  / float(int(c["n_dil_runs"])) - y_mean**(2.0) )**(0.5)

            # Make plot and calculate derivatives
            if k == "d":
                x1,y1,x2,y2,x_disc = calc_discontinuity(x_mean,y_mean)
                Tg[i] = x_disc
                if os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                    make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean,x1,x2],[y_mean,y1,y2],y_errs=[y_std,None,None],x_label="T (K)",y_label="d (g/cm^3)",color_inds=[0,3,3],\
                              linewidths=[5,2,2],v_line=[x_disc],labels=["mean","lin_left","lin_right"])
            if k == "H" and os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean],[y_mean],y_errs=[y_std],x_label="T (K)",y_label="H (kcal/mol)",labels=["mean"])
            if k == "U" and os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean],[y_mean],y_errs=[y_std],x_label="T (K)",y_label="U (kcal/mol)",labels=["mean"])                    
    return
                
# Calculates the Tg based on an algorithm that splits the data in two halves and fits the 
# low temperature region to a linear function and high temperature region to a quadratic function.
# Tg is taken as the split point that minimizes the error of the fit.
def calc_discontinuity(x,y):

    # Collected sorted values as arrays
    x,y = list(zip(*sorted([ (i,y[count_i]) for count_i,i in enumerate(x) ])))
    x = np.array(x)
    y = np.array(y)
    N = len(x)

    # Calculate linear fits on each half of the data
    best_err = None
    for i in range(2,N-2):

        # Left fit (linear only)
        A = np.array([ [ 1.0,_] for _ in x[:i] ])
        intercept1,slope1 = lstsq(A,y[:i],rcond=None)[0]
        err_low = np.mean((y[:i] - (slope1*x[:i] + intercept1))**(2.0))

        # Right fit (linear and quandratic terms)
        A = np.array([ [ 1.0,_,_**(2.0)] for _ in x[i:] ])
        intercept2,slope2,quad2 = lstsq(A,y[i:],rcond=None)[0]
        err_high = np.mean((y[i:] - (slope2*x[i:] + intercept2 + quad2*x[i:]**2.0 ))**(2.0))

        # Calculate total error with weighting based on the number of points in each set
        tot_err = err_low * float(i)/float(N) + err_high * float(N-i)/float(N)
        if best_err is None or tot_err < best_err:
            best_err = tot_err
            crossing_point = (intercept2-intercept1)/(slope1-slope2)
            crossing_ind = i
            b_s1 = slope1
            b_i1 = intercept1
            b_s2 = slope2
            b_i2 = intercept2

            b_q1 = 0
            b_q2 = quad2

    # Calculate the interpolated sets
    y1 = np.array([ (b_i1 + b_s1*x[j] + b_q1*x[j]**(2.0)) for j in range(crossing_ind) ])
    y2 = np.array([ (b_i2 + b_s2*x[j] + b_q2*x[j]**(2.0)) for j in range(crossing_ind,N) ])
    x1 = np.array([ x[j] for j in range(crossing_ind) ])
    x2 = np.array([ x[j] for j in range(crossing_ind,N) ])
    
    return x1,y1,x2,y2,x[crossing_ind]

# Fits linear data with a weight function
def disc_err(params,x,y,w):
    return np.mean(((params[0]*x+params[1]) - y)**2.0 * w)

def disc_quad(params,x,y,w):
    return np.mean((params[0] + params[1]*x + params[2]*x**(2) - y)**2.0 * x )

# Reads in column data from file
def get_cols(name,cols,head):
    data = { i:[] for i in cols }
    with open(name,'r') as f:
        for lc,lines in enumerate(f):
            if lc >= head:
                fields = lines.split()
                for i in cols:
                    data[i] += [float(fields[i])]
    return [ data[i] for i in cols ]

# Formats a subprocess string for execution by subprocess.check_output
def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

"""
# Function that sleeps the script until jobids are no longer in a running or pending state in the queue
def monitor_jobs(jobids):
    
    current_jobs = check_queue()
    while True in [ i in current_jobs for i in jobids ]:
        time.sleep(60)
        current_jobs = check_queue()  
    return
"""

# Returns the pending and running jobids for the user as a list
def check_queue():

    # The first time this function is executed, find the user name and scheduler being used. 
    if not hasattr(check_queue, "user"):

        # Get user name
        check_queue.user = subprocess.check_output("echo ${USER}", shell=True, encoding='utf8').strip("\r\n")

        # Get batch system being used

        squeue_tmp = subprocess.Popen(['which', 'squeue'], stdout=subprocess.PIPE,stderr=subprocess.STDOUT, encoding='utf8').communicate()[0].strip("\r\n")
        qstat_tmp = subprocess.Popen(['which', 'qstat'], stdout=subprocess.PIPE,stderr=subprocess.STDOUT, encoding='utf8').communicate()[0].strip("\r\n")
        check_queue.sched =  None
        if "no squeue in" not in squeue_tmp:
            check_queue.sched = "slurm"
        elif "no qstat in" not in qstat_tmp:
            check_queue.sched = "pbs"
        else:
            print("ERROR in check_queue: neither slurm or pbs schedulers are being used.")
            quit()

    # Get running and pending jobs using the slurm scheduler
    if check_queue.sched == "slurm":

        # redirect a squeue call into output
        output = subprocess.check_output("squeue -l", shell=True, encoding='utf8')

        # Initialize job information dictionary
        jobs = []
        id_ind = None
        for count_i,i in enumerate(output.split('\n')):            
            fields = i.split()
            if len(fields) == 0: continue                
            if id_ind is None and "JOBID" in fields:
                id_ind = fields.index("JOBID")
                if "STATE" not in fields:
                    print("ERROR in check_queue: Could not identify STATE column in squeue -l output.")
                    quit()
                else:
                    state_ind = fields.index("STATE")
                if "USER" not in fields:
                    print("ERROR in check_queue: Could not identify USER column in squeue -l output.")
                    quit()
                else:
                    user_ind = fields.index("USER")
                continue

            # If this job belongs to the user and it is pending or running, then add it to the list of active jobs
            if id_ind is not None and fields[user_ind] == check_queue.user and fields[state_ind] in ["PENDING","RUNNING"]:
                jobs += [fields[id_ind]]

    # Get running and pending jobs using the pbs scheduler
    elif check_queue.sched == "pbs":

        # redirect a qstat call into output
        output = subprocess.check_output("qstat -f", shell=True, encoding='utf8')

        # Initialize job information dictionary
        jobs = []
        job_dict = {}
        current_key = None
        for count_i,i in enumerate(output.split('\n')):
            fields = i.split()
            if len(fields) == 0: continue
            if "Job Id" in i:

                # Check if the previous job belongs to the user and needs to be added to the pending or running list. 
                if current_key is not None:
                    if job_dict[current_key]["State"] in ["R","Q"] and job_dict[current_key]["User"] == check_queue.user:
                        jobs += [current_key]
                current_key = i.split()[2]
                job_dict[current_key] = { "State":"NA" , "Name":"NA", "Walltime":"NA", "Queue":"NA", "User":"NA"}
                continue
            if "Job_Name" == fields[0]:
                job_dict[current_key]["Name"] = fields[2]
            if "job_state" == fields[0]:
                job_dict[current_key]["State"] = fields[2]
            if "queue" == fields[0]:
                job_dict[current_key]["Queue"] = fields[2]
            if "Resource_List.walltime" == fields[0]:
                job_dict[current_key]["Walltime"] = fields[2]        
            if "Job_Owner" == fields[0]:
                job_dict[current_key]["User"] = fields[2].split("@")[0]

        # Check if the last job belongs to the user and needs to be added to the pending or running list. 
        if current_key is not None:
            if job_dict[current_key]["State"] in ["R","Q"] and job_dict[current_key]["User"] == check_queue.user:
                jobs += [current_key]

    return jobs

# Plotting script
def make_plot(name,x_sets,y_sets,y_errs=[],dims=(4.5,4),x_label="x_label",y_label="y_label",labels=None,log_x=False,log_y=False,x_min=None,x_max=None,y_min=None,y_max=None,\
              tick_scale=1.0,x_tick_inc=None,y_tick_inc=None,diagonal=False,no_legend=False,linewidths=[3.0],color_inds=[],v_line=[]):
    
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.05,0.05,0.05),(0.9,0.6,0.05),(0.9,0.5,0.7),(0.35,0.7,0.9),(0.95,0.9,0.25)]   # Supposedly more CB-friendly
    if len(x_sets) != len(y_sets):
        print("ERROR in make_plot: len(x_sets) must equal len(y_sets)")
        quit()
    if False in [ len(i) == len(y_sets[count_i]) for count_i,i in enumerate(x_sets) ]:
        print("ERROR in make_plot: dimension mismatch in x and y sets.")
        quit()
    if labels is None:
        labels=[ i for i in range(len(x_sets)) ]
    elif len(labels) != len(x_sets):
        print("ERROR in make_plot: len(labels) is not equal to len(x_sets) and len(y_sets).")
        quit()
    elif len(y_errs) != 0 and len(y_errs) != len(y_sets):
        print("ERROR in make_plot: len(y_errs) is not equal to len(x_sets) and len(y_sets).")
        quit()
    while len(color_inds) < len(x_sets):
        color_inds += [len(color_inds)]
    while len(linewidths) < len(x_sets):
        linewidths += [linewidths[-1]]

    # Generate a plot with all individual datasets
    fig = plt.figure(figsize=dims)
    ax = plt.subplot(111)
    plot_handles = []

    # Plot error region
    for count_i,i in enumerate(y_errs):
        if i is not None:
            ax.fill_between(x_sets[count_i],\
                            (np.array(y_sets[count_i])-np.array(i)),\
                            (np.array(y_sets[count_i])+np.array(i)),\
                            facecolor=color_list[color_inds[count_i]],edgecolor='none',alpha=0.25)

    for count_i,i in enumerate(x_sets):
        plot_handle, = ax.plot(np.array(i),np.array(y_sets[count_i]),\
                               linestyle="-",linewidth=linewidths[count_i],\
                               color=color_list[color_inds[count_i]],\
                               label=labels[count_i])
        
        plot_handle.set_dash_capstyle('round')
        plot_handle.set_dash_joinstyle('round')
    #ax.set_ylabel(y_label.decode('utf-8'),fontsize=32,labelpad=10,fontweight='bold')
    #ax.set_xlabel(x_label.decode('utf-8'),fontsize=32,labelpad=10,fontweight='bold')
    ax.set_ylabel(y_label,fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(x_label,fontsize=32,labelpad=10,fontweight='bold')

    # Set Logarithmic Bounds
    if log_y:
        ax.set_yscale('log')
    if log_x:
        ax.set_xscale('log')
                               
    # Set limits
    c_y_min,c_y_max = ax.get_ylim()
    c_x_min,c_x_max = ax.get_xlim()
    if x_min is not None: c_x_min = x_min
    if x_max is not None: c_x_max = x_max
    if y_min is not None: c_y_min = y_min
    if y_max is not None: c_y_max = y_max
    ax.set_xlim([c_x_min,c_x_max])
    ax.set_ylim([c_y_min,c_y_max])

    # Set tick incremement
    if x_tick_inc is not None:
        plt.xticks(np.arange(c_x_min, c_x_max+x_tick_inc/1E8, x_tick_inc))
    if y_tick_inc is not None:
        plt.yticks(np.arange(c_y_min, c_y_max+y_tick_inc/1E8, y_tick_inc))
                               
    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=24*tick_scale,pad=10,direction='out',width=3,length=6)
    ax.tick_params(axis='both', which='minor',labelsize=24*tick_scale,pad=10,direction='out',width=2,length=4)
    [j.set_linewidth(3) for j in ax.spines.values()]
                               
    # Add diagonal
    if diagonal is True:
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c=".3")

    # Add vertical lines
    for i in v_line:
        ax.plot([i,i],ax.get_ylim(), ls="--", c=".3")

    # Put a legend to the right of the current axis and use modified save call to account for the legend artist
    if no_legend is False:

        # Place legend outside to the right by default
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))
        plt.savefig(name, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')            

    # If no legend is requested then place to use ordinary tight_layout and save calls
    else:
        plt.tight_layout()
        plt.savefig(name, bbox_inches=0,dpi=300)
    plt.close(fig)
    return

# Function for keeping tabs on the validity of the user supplied inputs
def parse_configuration(args):

    # Convert inputs to the proper data type
    if os.path.isfile(args.config) is False:
        print("ERROR in python_driver: the configuration file {} does not exist.".format(args.config))
        quit()

    # Process configuration file for keywords
    keywords = [ "FOM_GEOOPT_PROCS", "FOM_GEOOPT_WT", "FOM_GEOOPT_Q", "FOM_GEOOPT_SCHED", "FOM_GEOOPT_PPN", "FOM_GEOOPT_SIZE",\
                 "FOM_MD_GAS_PROCS", "FOM_MD_GAS_WT", "FOM_MD_GAS_Q", "FOM_MD_GAS_NPP", "FOM_MD_GAS_SCHED", "FOM_MD_GAS_PPN",\
                 "FOM_MD_GAS_SIZE", "FOM_MD_COND_PROCS", "FOM_MD_COND_WT", "FOM_MD_COND_Q", "FOM_MD_COND_NPP", "FOM_MD_COND_SCHED",\
                 "FOM_MD_COND_PPN", "FOM_MD_COND_SIZE", "FOM_TEMPS", "LAMMPS_EXE", "ORCA_EXE", "TAFFI_PATH", "fom_ff", 'batch', "FOM_N_TRAJ", "acct", \
                 "fom_t_ext", "fom_t_tot", "fom_t_a","fom_restarts", "fom_t_equil","gens" ]
    keywords = [ _.lower() for _ in keywords ]

    list_delimiters = [ "," ]  # values containing any delimiters in this list will be split into lists based on the delimiter
    space_delimiters = [ "&" ] # values containing any delimiters in this list will be turned into strings with spaces replacing delimiters
    configs = { i:None for i in keywords }    
    with open(args.config,'r') as f:
        for lines in f:
            fields = lines.split()
            
            # Delete comments
            if "#" in fields:
                del fields[fields.index("#"):]

            # Parse keywords
            l_fields = [ _.lower() for _ in fields ] 
            for i in keywords:
                if i in l_fields:

                    # Parse keyword value pair
                    ind = l_fields.index(i) + 1
                    if len(fields) >= ind + 1:
                        configs[i] = fields[ind]

                        # Handle delimiter parsing of lists
                        for j in space_delimiters:
                            if j in configs[i]:
                                configs[i] = " ".join([ _ for _ in configs[i].split(j) ])
                        for j in list_delimiters:
                            if j in configs[i]:
                                configs[i] = configs[i].split(j)
                                break
                                
                    # Break if keyword is encountered in a non-comment token without an argument
                    else:
                        print("ERROR in python_driver: enountered a keyword ({}) without an argument.".format(i))
                        quit()

    # Set defaults if None
    if configs["taffi_path"] is None:
        configs["taffi_path"] = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    if configs["fom_n_traj"] is None:
        configs["fom_n_traj"] = 5
    else:
        configs["fom_n_traj"] = int(configs["fom_n_traj"])
    if configs["fom_restarts"] is None:
        configs["fom_restarts"] = 0
    else:
        configs["fom_restarts"] = int(configs["fom_restarts"])
        if configs["fom_restarts"] < 0:
            print("ERROR: fom_restarts value must be >= 0.")
            quit()

    # Collect the xyz files
    configs["xyz"] = [ f for f in os.listdir('.') if os.path.isfile(f) and ".xyz" in f ]

    # Check that batch is an acceptable system
    if configs["batch"] not in ["pbs","slurm"]:
        print("ERROR in FOM_driver: only pbs and slurm are acceptable arguments for the batch variable.")
        quit()
    elif configs["batch"] == "pbs":
        configs["sub_cmd"] = "qsub"
    elif configs["batch"] == "slurm":
        configs["sub_cmd"] = "sbatch"

    return configs

if __name__ == "__main__":
    main(sys.argv[1:])
