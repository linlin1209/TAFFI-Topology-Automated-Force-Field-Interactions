#!/bin/env python                                                                                                                                                             
import sys,os,subprocess,argparse,time
from subprocess import PIPE

class monitor_jobs() :
      
   def qstat(self):

        #parser = argparse.ArgumentParser(description='This script interprets a qstat call and prints with better formatting/control.')

           #optional arguments
        #parser.add_argument('-w', dest='name_width', default=100,
        #                   help = 'The width of the jobname field (default: 100)')
        #parser.add_argument('-u', dest='user', default="",
        #                   help = 'user name (default: none)')

        # Assign stdin arguments to variables
        #args=parser.parse_args(arglist)
        #name_width = int(args.name_width)
        #user=str(args.user)
        
        name_width = int(100)
        user = str(self.user)


        # redirect a qstat call into output
        p = subprocess.Popen("qstat -f".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        output=str(output)
        # Initialize job information dictionary
        dict = {}
        jobs = []
        for count_i,i in enumerate(output.split(r'\n')):
                fields = i.split()
                if len(fields) == 0: continue
                if "Job Id" in i:
                  current_key = i.split()[2]
                  dict[current_key] = { "State":"NA" , "Name":"NA", "Walltime":"NA", "Queue":"NA", "User":"NA"}
                  continue
                if "Job_Name" == fields[0]:
                  dict[current_key]["Name"] = fields[2]
                if "job_state" == fields[0]:
                  dict[current_key]["State"] = fields[2]
                if "queue" == fields[0]:
                  dict[current_key]["Queue"] = fields[2]
                if "Resource_List.walltime" == fields[0]:
                  dict[current_key]["Walltime"] = fields[2]        
                if "Job_Owner" == fields[0]:
                  dict[current_key]["User"] = fields[2].split("@")[0]
        self.dict = dict
        #template_string="{:<40s} {:<"+str(name_width)+"s} {:<20s} {:20s} {:<20s}"
        #print (template_string.format("ID","Name","User","Walltime","Status","Queue"))
        for i in dict.keys():
                if dict[i]["State"] in ["R","Q"] and (user == "" or dict[i]["User"] == user):
                   jobs += [i]
                   #print (template_string.format(i,dict[i]["Name"],dict[i]["User"],dict[i]["Walltime"],dict[i]["State"],dict[i]["Queue"]))

        return jobs

   def main(self):
      current_jobs = self.qstat()
      print("waiting for jobs {} to complete....".format(self.jobids))
      while True in [ i in current_jobs for i in self.jobids ]:
        time.sleep(60)
        current_jobs = self.qstat()  
      return

   def __init__(self,jobids,user):
      self.dict = {}
      self.jobids = jobids
      self.user = user
      self.main()


if __name__ == "__main__":
   #monitor_jobs(sys.argv[1:],'kim2096') 
   monitor_jobs(['5610892.brown-adm.rcac.purdue.edu', '5610893.brown-adm.rcac.purdue.edu'],'kim2096') #for test
