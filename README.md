# TAFFI-Topology-Automated-Force-Field-Interactions
 A Framework for Developing Transferable Force Fields
### This is the code for https://doi.org/10.1021/acs.jcim.1c00491 and has only been uploaded for company's viewing
## Software dependency:
1. ORCA
2. LAMMPS
3. SLURM submission system (driver for other job scheduler system such as PBS can be requested)
## How to run:
1. put the xyz of all the molecules that you want to parametrize in your executing folder
2. prepare your configuration file (there is an example config.txt in the Automation_scripts folder)
3. execute "driver.py -c your_config_file"
4. Wait for all jobs to complete, then you should find your desire parameters in Param_for_batch.db
