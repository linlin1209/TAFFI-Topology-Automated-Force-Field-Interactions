
import sys,os,argparse
import shutil
from numpy import *


def main(argv):

    parser = argparse.ArgumentParser(description='Restarts  in a .xyz file and generates smallest atomtype-consistent fragments for parametrizing each missing force-field mode.')

    #required (positional) arguments                                                                                                  
    parser.add_argument('folder', help = 'Folder for recursive search of flexible dihedral files in need of resubmission.')
   
    # optional arguments
    parser.add_argument('--post', dest='post_process', default=False, action="store_const", const=True,\
                        help = 'Runs in post-processing mode.')

    parser.add_argument('--verbose', dest='verbose', default=False, action="store_const", const=True,\
                        help = 'When present, the program will print the names of each file that is being restarted.')

    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    

    # Find the dihedral folders
    dihedral_folders = [ dp+'/'+dirs for dp,dn,fn in os.walk(args.folder) for dirs in dn if dirs in ["scan-forward","scan-reverse"] ] 
    
    # Check for job completion
    incomplete_list = []
    for i in dihedral_folders:
        output_name = next( j for j in os.listdir(i) if j.split('.')[-1] == "out" )
        base_name   = '.'.join(output_name.split('.')[:-1])
        completed_flag = False
        with open("{}/{}".format(i,output_name),'r') as f:
            for lines in f:
                if "****ORCA TERMINATED NORMALLY****" in lines:
                    completed_flag = True
        if completed_flag is False:
            incomplete_list += [(i,base_name)]

    # Print the files being restarted if the flag is present.
    if args.verbose is True:
        for i in incomplete_list:
            print(("{}/{}.in".format(i[0],i[1])))
        
    # Run pre-processing/job preparation
    if args.post_process is False:        

        # Loop over the incomplete files/folders and generate the resubmission jobs
        for i in incomplete_list:

            # Check for pre-requisites
            if os.path.isfile("{}/{}.in".format(i[0],i[1])) is False:
                print("ERROR: Could not restart {}/{}.in because the input file is missing.".format(i[0],i[1]))
            if os.path.isfile("{}/{}.out".format(i[0],i[1])) is False:
                print("ERROR: Could not restart {}/{}.in because the output file is missing.".format(i[0],i[1]))
                
            # Perform pre-processing
            pre_process("{}/{}.in".format(i[0],i[1]))

    # Run post-processing/merging of new results
    else:

        # Loop over the incomplete files/folders and process the results
        for i in incomplete_list:
            
            # Perform pre-processing
            post_process("{}/{}.in".format(i[0],i[1]))


    return

# Description: This function performs a job restart for flexible dihedral scans
def pre_process(filename,out_name="RESTART"):

    # Process the name and path
    name = filename.split('/')[-1]
    path_to_file = '/'.join(filename.split('/')[:-1])
    if len(path_to_file) != 0:
        path_to_file = path_to_file+'/'

    # Check the existence of the path
    if os.path.isfile(path_to_file+name) is False:
        print("ERROR in pre_process: {} does not exist. Exiting...".format(filename))
        quit()

    # Create the restart folder (out_name) if required
    if os.path.isdir(path_to_file+out_name) is False:        
        os.mkdir(path_to_file+out_name)
    else:
        if os.path.isfile(path_to_file+out_name+'/tmp.out') is True and os.path.isfile(path_to_file+out_name+"/RESTART.in") is True:
            return
        elif os.path.isfile(path_to_file+out_name+'/tmp.out') is True:
            os.remove(path_to_file+out_name+'/tmp.out')
        elif os.path.isfile(path_to_file+out_name+"/RESTART.in") is True:
            os.remove(path_to_file+out_name+'/RESTART.in')

    # Parse the run details from the input file
    flag = 0
    with open(path_to_file+name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) == 9 and fields[0] == "D":
                flag = 1                                            # toggle indicates constraint data was found                                                
                D_ind = (fields[1],fields[2],fields[3],fields[4])   # indices of the constrained dihedral
                start = float(fields[6].strip(','))                 # starting angle for the scan
                end = float(fields[7].strip(','))                   # ending angle for the scan
                N_steps = int(fields[8].strip('.'))                 # number of steps in the scan
                inc = (end - start)/(float(N_steps)-1.0)            # step size, in degrees, for the scan
                D_line = lc                                         # line count for the constraint information
                break
            if len(fields) > 1 and fields[0] == "%base":
                base_name = fields[1].strip("'").strip("\"")        # base_name for the output files from this job

    # Parse how far the job progressed from the output file, and the break line for the last complete run.
    cycles,break_line = number_complete(path_to_file+'.'.join(name.split('.')[:-1])+'.out',break_line_opt=True)

    # If all the cycles are complete then there is no need for a resubmission.
    if cycles == N_steps:
        return

    # Grab new geometry
    E,G = xyz_parse(path_to_file+base_name+".{}.xyz".format(str(cycles).zfill(3)))

    # Parse information from the original job
    orca_dict = orca_in_parse(path_to_file+name)

    # Write the new input file
    with open(path_to_file+out_name+'/{}.in'.format(out_name),'w') as f:        
        f.write("#Run DFT dihedral scan\n")
        f.write("! "+orca_dict["0"]["header_commands"])
        f.write("\n%base \"{}\"\n\n".format(base_name))
        f.write("* xyz {:d} {:d}\n".format(int(orca_dict["0"]["charge"]),int(orca_dict["0"]["multiplicity"])))
        for count_i,i in enumerate(E):
            f.write("  {:<20s} {:<20.8f} {:<20.8f} {:<20.8f}\n".format(i,G[count_i,0],G[count_i,1],G[count_i,2]))
        f.write("*\n")                

    # Add dihedral constraint block to the input file
    with open(path_to_file+'.'.join(name.split('.')[:-1])+'.in','r') as f:
        with open(path_to_file+out_name+'/{}.in'.format(out_name),'w') as o:        
            for lc,lines in enumerate(f):
                fields = lines.split()
                if lc == D_line:
                    o.write("    D {} {} {} {} = {}, {}, {}\n".format(D_ind[0],D_ind[1],D_ind[2],D_ind[3],inc*float(cycles)+start,end,N_steps-cycles))
                elif len(fields) == 4 and fields[0] == "*" and fields[1] == "xyz":
                    o.write("* xyz {:d} {:d}\n".format(int(orca_dict["0"]["charge"]),int(orca_dict["0"]["multiplicity"])))           
                    for count_i,i in enumerate(E):
                        o.write("  {:<20s} {:<20.8f} {:<20.8f} {:<20.8f}\n".format(i,G[count_i,0],G[count_i,1],G[count_i,2]))
                    o.write("*\n")                
                    break        
                else:
                    o.write(lines)

    # Write the truncated output file
    with open(path_to_file+'.'.join(name.split('.')[:-1])+'.out','r') as f:
        with open(path_to_file+'{}/tmp.out'.format(out_name),'w') as o:
            for lc,lines in enumerate(f):
                o.write(lines)
                if lc == break_line:
                    break
    
    return

# Description: Used to post process a restarted dihedral scan. This function is used
#              after the script has already been used in the pre-processing mode. 
def post_process(filename,out_name="RESTART"):

    # Process the name and path
    name = filename.split('/')[-1]
    path_to_file = '/'.join(filename.split('/')[:-1])
    if len(path_to_file) != 0:
        path_to_file = path_to_file+'/'

    # Check the existence of the path
    if os.path.isfile(path_to_file+name) is False:
        print("ERROR: {} does not exist. Exiting...".format(filename))
        quit()

    # Create the out_name folder if required
    if os.path.isdir(path_to_file+out_name) is False:        
        print("ERROR in post_process: {}{} does not exist. You need to run this script in pre-process mode first.".format(path_to_file,out_name))

    # Parse the run details from the input file
    flag = 0
    break_line = 0
    with open(path_to_file+name,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) == 9 and fields[0] == "D":
                flag = 1                                            # toggle indicates constraint data was found                                                
                D_ind = (fields[1],fields[2],fields[3],fields[4])   # indices of the constrained dihedral
                start = float(fields[6].strip(','))                 # starting angle for the scan
                end = float(fields[7].strip(','))                   # ending angle for the scan
                N_steps = int(fields[8].strip('.'))                 # number of steps in the scan
                inc = (end - start)/(float(N_steps)-1.0)            # step size, in degrees, for the scan
                D_line = lc                                         # line count for the constraint information
                break
            if len(fields) > 1 and fields[0] == "%base":
                base_name = fields[1].strip("'").strip("\"")        # base_name for the output files from this job

    # Determine the number of cycles that completed in the original calculation
    cycles         = number_complete(path_to_file+out_name+'/tmp.out')

    # Determine the number of cycles that completed in the second calculation
    cycles_restart = number_complete(path_to_file+out_name+'/RESTART.out')

    # Rename the files
    for i in range(1,cycles_restart+1):
       if os.path.isfile("{}{}.{}.xyz".format(path_to_file,base_name,str(i).zfill(3))):
           shutil.move("{}{}/{}.{}.xyz".format(path_to_file,out_name,base_name,str(i).zfill(3)),"{}{}.{}.xyz".format(path_to_file,base_name,str(i+cycles).zfill(3)))

       if os.path.isfile("{}{}.{}.gbw".format(path_to_file,base_name,str(i).zfill(3))):
           shutil.move("{}{}/{}.{}.gbw".format(path_to_file,out_name,base_name,str(i).zfill(3)),"{}{}.{}.gbw".format(path_to_file,base_name,str(i+cycles).zfill(3)))


    # Update the output file
    counter = 0
    with open("{}{}/RESTART.out".format(path_to_file,out_name),'r') as f:
        with open("{}{}/tmp.out".format(path_to_file,out_name),'a') as o:
            for lc,lines in enumerate(f):
                fields = lines.split()
                if len(fields) == 7 and fields[0] == "*" and fields[1] == "RELAXED" and fields[2] == "SURFACE" and fields[3] == "SCAN" and fields[4] == "STEP":

                    if counter == 0:
                        o.write("         {}\n".format("*"*61))
                
                    o.write("         {}\n".format("*{:^59s}*".format("RELAXED SURFACE SCAN STEP {:>3s}".format(str(counter+cycles+1)))))
                    counter += 1
                    continue
                if counter > 0:
                    o.write(lines)

    shutil.move("{}{}/tmp.out".format(path_to_file,out_name),path_to_file+'.'.join(name.split('.')[:-1])+'.out')
    shutil.rmtree(path_to_file+out_name)
    return

# parses the number of cycles completed before termination
def number_complete(filename,break_line_opt=False):

    # Parse the number of complete steps in the scan, based on the output file.
    cycles = 0
    with open(filename,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) == 5 and "*** OPTIMIZATION RUN DONE ***" in lines:
              cycles += 1

    # If the break_line is not needed then just return the cycles value
    if break_line_opt is False:
        return cycles

    # Parse last line that should be kept
    cycle_count = 0
    break_line  = 0
    with open(filename,'r') as f:
        for lc,lines in enumerate(f):
            fields = lines.split()

            if len(fields) == 5 and "*** OPTIMIZATION RUN DONE ***" in lines:
              cycle_count += 1

            if cycle_count == cycles:
                if len(fields) == 7 and fields[0] == "*" and fields[1] == "RELAXED" and fields[2] == "SURFACE" and fields[3] == "SCAN": 
                    break_line = lc - 2

    # If the break line wasn't assigned then if is because the job never even started the cycles+1 iteration
    # In this case the whole file is retained.
    if break_line == 0:
        break_line = lc

    return cycles,break_line

# Description: Simple wrapper function for grabbing the coordinates and
#              elements from an xyz file
#
# Inputs      input: string holding the filename of the xyz
# Returns     Elements: list of element types (list of strings)
#             Geometry: Nx3 array holding the cartesian coordinates of the
#                       geometry (atoms are indexed to the elements in Elements)
#
def xyz_parse(input):

    # Iterate over the remainder of contents and read the
    # geometry and elements into variable. Note that the
    # first two lines are considered a header
    with open(input,'r') as f:
        for lc,lines in enumerate(f):
            fields=lines.split()

            # Parse header
            if lc == 0:
                if len(fields) < 1:
                    print("ERROR in xyz_parse: {} is missing atom number information".format(input))
                    quit()
                else:
                    N_atoms = int(fields[0])
                    Elements = ["X"]*N_atoms
                    Geometry = zeros([N_atoms,3])
                    count = 0

            # Parse body
            if lc > 1:

                # Skip empty lines
                if len(fields) == 0:
                    continue            

                # Write geometry containing lines to variable
                if len(fields) > 3:

                    # Consistency check
                    if count == N_atoms:
                        print("ERROR in xyz_parse: {} has more coordinates than indicated by the header.".format(input))
                        quit()

                    # Parse commands
                    else:
                        Elements[count]=fields[0]
                        Geometry[count,:]=array([float(fields[1]),float(fields[2]),float(fields[3])])
                        count = count + 1

    # Consistency check
    if count != len(Elements):
        print("ERROR in xyz_parse: {} has less coordinates than indicated by the header.".format(input))
        
    return Elements,Geometry

# Description: Parses keywords and geometry block from an orca input file
#
# Inputs        input: string holding the filename of the orca input file
# Returns       orca_dict: dictionary holding the run information for each job in the input file
#                          the first key in the dictionary corresponds to the job number (i.e.,
#                          orca_dict["0"] references the information in the first job. The job info
#                          can be accessed with content specific keys ("header_commands", "geo", 
#                          "elements", "constraints", "N_proc", "job_name", "charge", "multiplicity",
#                          "content", "geom_block" )
def orca_in_parse(input):

    # Iterate over the contents and return a dictionary of input components indexed to each job in the input file
    job_num = 0
    orca_dict = {str(job_num):{"header_commands": "","content": "","elements": None, "geo": None, "constraints": None, "geo_opts_block": None, "job_name": None}}
    geo_opts_flag = 0
    geo_block_flag = 0
    con_flag  = 0
    
    # Open the file and begin the parse
    with open(input,'r') as f:
        for lc,lines in enumerate(f):

            # Grab fields 
            fields = lines.split()            
            
            # Update the "content" block, which contains everything
            orca_dict[str(job_num)]["content"] += lines

            # If a new job is encountered reset all flags and update the job_num counter            
            if len(fields) > 0 and fields[0] == "$new_job":
                job_num += 1
                con_flag = 0
                geo_opts_flag = 0
                geo_block_flag = 0
                orca_dict[str(job_num)] = {"header_commands": "","content": "","elements": None, "geo": None, "constraints": None, "geo_opts_block": None, "N_proc": orca_dict[str(job_num-1)]["N_proc"]}

            # Component based parse commands
            if len(fields) > 0 and fields[0] == "!":
                orca_dict[str(job_num)]["header_commands"] += " ".join(fields[1:]) + " "
                if "PAL" in lines:
                    orca_dict[str(job_num)]["N_proc"] = int([ i.split("PAL")[1] for i in fields if "PAL" in i ][0])
                elif job_num != 0:
                    orca_dict[str(job_num)]["N_proc"] = orca_dict[str(job_num-1)]["N_proc"]
                else:
                    orca_dict[str(job_num)]["N_proc"] = 1                    
            if len(fields) > 0 and fields[0] == "%base":
                orca_dict[str(job_num)]["job_name"] = fields[1]
                
            # Check for turning on flags
            if len(fields) > 0 and fields[0] == "%geom":
                geo_opts_flag = 1
                orca_dict[str(job_num)]["geo_opts_block"] = ""                
                continue
            if len(fields) > 0 and fields[0] == "Constraints":
                if geo_opts_flag == 1:
                    orca_dict[str(job_num)]["geo_opts_block"] += lines                
                con_flag = 1
                orca_dict[str(job_num)]["constraints"] = ""
                continue
            if len(fields) >= 2 and fields[0] == "*" and fields[1] == "xyz":
                geo_block_flag = 1
                orca_dict[str(job_num)]["charge"] = float(fields[2])
                orca_dict[str(job_num)]["multiplicity"] = int(fields[3])
                orca_dict[str(job_num)]["geo"] = []
                orca_dict[str(job_num)]["elements"] = []
                continue
            if len(fields) >= 2 and fields[0] == "*" and fields[1] == "xyzfile":
                orca_dict[str(job_num)]["charge"] = float(fields[2])
                orca_dict[str(job_num)]["multiplicity"] = int(fields[3])                
                orca_dict[str(job_num)]["geo"] = None
                orca_dict[str(job_num)]["elements"] = None
                continue

            # Checks for turning off flags
            if con_flag == 1 and len(fields) > 0 and fields[0] == "end":
                con_flag = 0
                continue
            if geo_opts_flag == 1 and len(fields) > 0 and fields[0] == "end":
                geo_opts_flag = 0            
                continue
            if geo_block_flag == 1 and len(fields) > 0 and fields[0] == "*":
                geo_block_flag = 0            
                continue
            
            # Flag based parse commands
            if geo_opts_flag == 1:
                orca_dict[str(job_num)]["geo_opts_block"] += lines
            if con_flag == 1:
                orca_dict[str(job_num)]["constraints"] += lines
            if geo_block_flag == 1:
                orca_dict[str(job_num)]["geo"] += [ [ float(i) for i in fields[1:] ] ]
                orca_dict[str(job_num)]["elements"] += [ str(fields[0]) ]
            
    return orca_dict

# Main sentinel
if __name__ == "__main__":
    main(sys.argv[1:])
