#!/bin/env python                                                                                                                                                             

import sys,argparse,os,datetime,fnmatch,re
import random,math
import numpy as np


def main(argv):

     parser = argparse.ArgumentParser(description='Reads in lammps data file then transform it into gro file for viewing')


     #required (positional) arguments                                                                                                  
     #parser.add_argument('data_file', help = 'lammps data file')
                                            

     #optional arguments    
     parser.add_argument('-o', dest='outputname', default='',
                        help = 'output file name')

     parser.add_argument('--boxlo', dest='halfbox', default=0, action='store_const', const=1,
                        help = 'move molecule from lower box boundary to 0  (default: off)')

     # Make parse inputs
     args=parser.parse_args(argv)

     xi,v = get_v_xi('back.lammpstrj',700)
     f = open('v_xi_time.txt','w')
     f.write('{:<20s} {:<20s} {:<20s}\n'.format('time','xi','v'))
     for t in xi: 
         f.write('{:< 20.5f} {:< 20.5f} {:< 20.5f}\n'.format(t,xi[t],v[t]))
     f.close()
     quit()
     #"""
         
         
             
         
     for i in range(0,num_jobs):
         folder = str(i)
         if os.path.isdir(folder) is False:
            os.makedirs(folder)
         os.chdir(folder)
         v = sample_velocity(126.90447,298)
         write_data('../data_from_lammpstrj.data',-v,700)
         os.chdir('..')
         print(v)
     write_submit(128,24,num_jobs)
     quit()

     if childCount < childTrajectories:

         
         # Seed the random number generator
         self.initializeRandomNumberGenerator()

         # Generate initial position using transition state geometry
         # XXX: need to be filled, this is currently done manually
         
         # Equilibrate parent trajectory while constraining to dividing surface
         logging.info('Equilibrating parent trajectory for {0:g} ps...'.format(equilibrationSteps * self.dt * 2.418884326505e-5))
         # XXX: need to be filled, this is currently done manually
         
         logging.info('Finished equilibrating parent trajectory.')
         logging.info('')
     
         # Continue evolving parent trajectory, interrupting to sample sets of
         # child trajectories in order to update the recrossing factor
         parentIter = 0
         while childCount < childTrajectories:
             
             logging.info('Sampling {0} child trajectories at {1:g} ps...'.format(childrenPerSampling, parentIter * childSamplingSteps * self.dt * 2.418884326505e-5))
 
             # Sample a number of child trajectories using the current parent
             # configuration
             results = []
             for child in range(childrenPerSampling / 2):
                 p_child,v_child = sample_momentum(mass,T)
                 
                 childCount += 2
 
             for child in range(childrenPerSampling / 2):
                 # This line will block until the child trajectory finishes
                 if pool:
                     num, denom = results[child].get()
                 else:
                     num, denom = results[child]
                 # Update the numerator and denominator of the recrossing factor expression
                 kappa_num += num
                 kappa_denom += denom
         
             logging.info('Finished sampling {0} child trajectories at {1:g} ps.'.format(childrenPerSampling, parentIter * childSamplingSteps * self.dt * 2.418884326505e-5))
             
             self.saveRecrossingFactor(recrossingFilename, kappa_num, kappa_denom, childCount,
                 childTrajectories, equilibrationSteps, childSamplingSteps, childEvolutionSteps, childrenPerSampling)
             
             logging.info('Current value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
             logging.info('')
                             
             # Further evolve parent trajectory while constraining to dividing
             # surface and sampling from Andersen thermostat
             logging.info('Evolving parent trajectory to {0:g} ps...'.format((parentIter+1) * childSamplingSteps * self.dt * 2.418884326505e-5))
             result = system.equilibrate(0, p, q, childSamplingSteps, self.xi_current, self.potential, 0.0, True, saveParentTrajectory)
             while result != 0:
                 q = numpy.asfortranarray(q0.copy())
                 p = self.sampleMomentum()            
                 result = system.equilibrate(0, p, q, equilibrationSteps, self.xi_current, self.potential, 0.0, True, saveParentTrajectory)
             
             parentIter += 1
         
         logging.info('Finished sampling of {0:d} child trajectories.'.format(childCount))
         logging.info('')

def runRecrossingTrajetory(equil_data,v,target_index):
    """
    Run an individual pair of recrossing factor child trajectories, returning
    the contributions to the numerator and denominator of the recrossing factor
    from this trajectory pair. We use pairs of trajectories so that we always
    sample in the positive and negative directions of the initial sampled
    momenta.
    """
    # Trajectory for the negative of the sampled momenta
    v1 = v*(-1)
    kappa_num1, kappa_denom1 = recrossing_trajectory(equil_data,v1,target_index)

    # Trajectory for the positive of the sampled momenta
    v2 = v
    kappa_num2, kappa_denom2 = recrossing_trajectory(equil_data,v2,target_index)

    return kappa_num1 + kappa_num2, kappa_denom1 + kappa_denom2

#def recrossing_trajetory(equil_data,v1,target_index):

def get_v_xi(data_name,target_index):
   v = {}
   xi = {}
   print('reading {}...'.format(data_name),end='\r')
   with open(data_name,'r') as f:
      flag = 0
      for lc,lines in enumerate(f):
         fields = lines.split()
         if flag == 0 and len(fields) < 2: continue
         if fields[0] == 'ITEM:' and fields[1] == 'TIMESTEP':
            flag = 1
            continue

         if fields[0] == 'ITEM:' and fields[1] == 'ATOMS' and fields[2] == 'id':
            flag  = 2
            continue

         if flag == 1:
            time = float(fields[0])
            flag = 0
            continue 
            
         if flag == 2:
            index = int(fields[0])
            if index == target_index:
               xi[time] = float(fields[4])
               v[time] = float(fields[7])
            continue

   return xi,v 
   

# replace the ion velocity to the newly generated one
def write_data(equil_data,v,target_index):   
   if os.path.isfile(equil_data) is False:
      print('equil data after constrained equilibration not found')
      quit()
   g = open('init.data','w')
   with open(equil_data,'r') as f:
      flag = 0
      for lc,lines in enumerate(f):
         fields = lines.split()
         if flag == 0: 
            g.write(lines)

         if len(fields)>0 and fields[0] == 'Velocities':
            flag = 1
            tag = 0
            continue
            
         if flag == 1 :
            if len(fields) == 0: 
               g.write('\n')
               tag += 1
            else:
               index = int(fields[0])
               if index == target_index: 
                  fields[3] = '{:20.13e}'.format(v)
                  #for count_i,i in enumerate(range(1,4)):
                  #   fields[i] = '{:20.13e}'.format(v[count_i])
               g.write('{}\n'.format(' '.join(fields)))
               continue
            if tag == 2: 
               flag = 0
            continue
   
   return  

def write_submit(total_cpu,cpu,num_jobs):
   ori_dir = os.getcwd()
   jobs_per_node = int(total_cpu/cpu)
   N_nodes = int(num_jobs/jobs_per_node)
   job = 0
   for i  in range(N_nodes):
      with open('submit.{}.sh'.format(i),'w') as f:
         f.write("#!/bin/bash\n")
         f.write("#\n")
         f.write("#SBATCH --job-name submit_{}\n".format(i))
         f.write("#SBATCH -o submit.{}.out\n".format(i))
         f.write("#SBATCH -e submit.{}.err\n".format(i))
         f.write("#SBATCH -A standby\n")
         f.write("#SBATCH -N 1\n")
         f.write("#SBATCH -n {}\n".format(int(cpu*jobs_per_node)))
         f.write("#SBATCH -t 04:00:00\n")
         for j in range(jobs_per_node):
            f.write('\ncd {}/{}\n'.format(ori_dir,job))
            f.write('mpirun -np {} /depot/bsavoie/apps/lammps/exe/lmp_mpi_180501 -in npt.in.init > npt.out &\n\n'.format(cpu))
            job += 1
         f.write('wait\n')
   if job < num_jobs:
      with open('submit.{}.sh'.format(N_nodes),'w') as f:
         f.write("#!/bin/bash\n")
         f.write("#\n")
         f.write("#SBATCH --job-name submit_{}\n".format(i))
         f.write("#SBATCH -o submit.{}.out\n".format(N_nodes))
         f.write("#SBATCH -e submit.{}.err\n".format(N_nodes))
         f.write("#SBATCH -A standby\n")
         f.write("#SBATCH -N 1\n")
         f.write("#SBATCH -n {}\n".format(int(cpu*jobs_per_node)))
         f.write("#SBATCH -t 04:00:00\n")
         while job < num_jobs:
            f.write('\ncd {}/{}\n'.format(ori_dir,job))
            f.write('mpirun -np {} /depot/bsavoie/apps/lammps/exe/lmp_mpi_180501 -in npt.in.init > npt.out &\n\n'.format(cpu))
            job += 1
         f.write('wait\n')
   return
      

def sample_velocity(mass,T):
   
   mass = mass/1000 # kg/mol
   p = [0,0,0]
   kB = 1.3806504e-23*6.02e23 # J/(molK)
   beta = 1 / (kB*T)
   v = math.sqrt(-2*kB*T/mass*np.log(random.uniform(0,1)))
   v = v*1e-5 # m/s -> A/fs
   return v



def sample_momentum(mass,T):
   
   p = [0,0,0]
   kB = 1.3806504e-23*6.02e23/4.184
   beta = 1 / (kB*T)
   dp = math.sqrt(mass/beta)
   for i in range(0,3):
      p[i] = randomn()*dp # momentum
   v = [ i/mass for i in p] # velocity

   return p,v
      
def randomn():

   # Marsaglia polar method
   x = 1
   y = 1
   s = 1
   while (s >= 1):
        x = random.uniform(0,1)
        y = random.uniform(0,1)
        s = x**2 + y**2
   ins = math.sqrt(-2.0 * math.log(s) / s)
   g1 = x * ins
   g2 = y * ins

   return g1



class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gro_to_lmp.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

    def flush(self):
        pass


if __name__ == "__main__":
   main(sys.argv[1:])
