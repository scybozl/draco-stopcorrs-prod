import os
import sys
import time

def listdirs(d):

	top = []
	for f in os.listdir(d):
	  if f.find("pwg-")!=-1 and f.find("-NLO.top")!=-1 and f.find(".outlier") == -1:
	    top += [f]
	return top


path = sys.argv[1]
cwd = os.getcwd()
l = len(os.listdir(path))

for i,dirs in enumerate(os.listdir(path)):

  curdir = cwd+'/'+path+dirs
#  os.chdir(curdir)
  hand = dirs.split('-')[0]
  PRmss = float(dirs.split('-Mst.')[1].split('-M')[0])
  PRcim = float(dirs.split('-M1.')[1].split('-mt')[0])

  arg = "Results.stop_bino_mssm_LO_CT14_tst"+hand.lower()+"-"+str(PRmss)+"-"+str(PRcim)

  print i+1, "/", l, " ---> Cross-sections for ", hand, ", mstop = ", PRmss, ", mchi = ", PRcim
#  print "srun "+cwd+"/calcxsec " + arg +" -f RunStops_MSSM_CT14_tst"+hand.lower()+"-merged.dat PRparamcard:=t"+hand+"."+str(int(PRmss))+"-M1."+str(int(PRcim))+"-mt.173.1.dat PRmss:="+str(PRmss)+" PRcim:="+str(PRcim)+" PRshower:=CSS PRtgerr:=999"

  runname = hand+'-'+str(PRmss)+'-'+str(PRcim)
  batchname = 'qsub/submit-'+runname+'.sh'
  batchfile = open(batchname, 'w')
  batchfile.write('#!/bin/bash -l\n')
  batchfile.write('#SBATCH --error="'+cwd+'/qsub/out/slurm-'+runname+'.err"\n')
  batchfile.write('#SBATCH --output="'+cwd+'/qsub/out/slurm-'+runname+'.out"\n')
  batchfile.write('#SBATCH -J stop-xsecs\n')
  batchfile.write('#SBATCH --partition=express\n')
  batchfile.write('#SBATCH --time=00:30:00\n')
  batchfile.write('cd '+curdir+'\n')
  batchfile.write('source '+cwd+'/setup.sh\n')
  batchfile.write("srun "+cwd+"/calcxsec " + arg +" -f RunStops_MSSM_CT14_tst"+hand.lower()+"-merged.dat PRparamcard:=t"+hand+"."+str(int(PRmss))+"-M1."+str(int(PRcim))+"-mt.173.1.dat PRmss:="+str(PRmss)+" PRcim:="+str(PRcim)+" PRshower:=CSS PRtgerr:=999\n")
  batchfile.close()

  os.system('sbatch '+batchname)

  while True:
    njobs = int(os.popen("squeue -u lscyboz -n 'stop-xsecs' | wc -l").read())
    print 'njobs = ', njobs-1, ' ...'
    if (njobs < 20): break
    time.sleep(20)

#  print "srun "+cwd+"/calcxsec " + arg +" -f RunStops_MSSM_CT14_tst"+hand.lower()+"-merged.dat PRparamcard:=t"+hand+"."+str(int(PRmss))+"-M1."+str(int(PRcim))+"-mt.173.1.dat PRmss:="+str(PRmss)+" PRcim:="+str(PRcim)+" PRshower:=CSS PRtgerr:=999"
#  os.system("srun "+cwd+"/calcxsec " + arg +" -f RunStops_MSSM_CT14_tst"+hand.lower()+"-merged.dat PRparamcard:=t"+hand+"."+str(int(PRmss))+"-M1."+str(int(PRcim))+"-mt.173.1.dat PRmss:="+str(PRmss)+" PRcim:="+str(PRcim)+" PRshower:=CSS PRtgerr:=999")
