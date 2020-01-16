from subprocess import Popen, PIPE #, check_output
from sys import argv, exit
from glob import glob
from datetime import date
import time
import os

mode = argv[1] # grid or event or resubmit
resubmitFolder = argv[2]
ps = 'CSS'

submit = 'slurm' # choose ['serial','slurm','condor']
nstart = 1 # first job number
n = 32 # number of jobs

if mode == "grid":
    err = 0.01
elif mode == "event":
    err = 0.1
elif mode == "rivet":
    err = "err"

analyses = ['-a MC13TeV_SPINCORR_emu_TOP']
today = date.today().strftime("%b-%d-%Y")
cwd = os.getcwd()

if mode == "grid":

  if os.path.exists('qsub/'+today): exit()
  os.system('mkdir qsub/'+today+'\n')
  os.system('mkdir qsub/out/'+today+'\n')

  inputlist = open("ParamCards/tL.half.list", 'r')
  for pcard in inputlist:
    pcard = pcard.split('\n')[0]
#  for pcard in glob("ParamCards/tL*.dat"):
    for hand in ['L', 'R']:

#        paramcard = pcard.split('/')[1]
        paramcard = pcard.replace('tL', 't'+hand)
        mstop = paramcard.split('.')[1].split('-')[0]
        mchi  = paramcard.split('-')[1].split('M1.')[1]

        print " ------- mstop = ", mstop, " - mchi = ", mchi, " -- ", hand, " ------- "

        runname = hand+'-Mst.'+mstop+'-M1.'+mchi+'-mt.173.1'
        dirname = cwd+'/'+today+'/'+runname
        filename = 'RunStops_MSSM_CT14_tst'+hand.lower()+'-merged.dat'

        os.system('mkdir -p '+dirname)
        os.system('cp '+'ParamCards/'+paramcard+' '+dirname)
        os.system('cp base/'+filename+' '+dirname)

        sherpaopts = [ "PRtgerr:="+str(err),
                       "PRmss:="+str(float(mstop)),
                       "PRcim:="+str(float(mchi)),
                       "PRshower:="+ps,
                       "PRparamcard:=t"+hand+"."+mstop+"-M1."+mchi+"-mt.173.1.dat",
                       "-f ", filename ]

        batchname = 'qsub/'+today+'/submit-'+runname+'.sh'
        batchfile = open(batchname, 'w')
        batchfile.write('#!/bin/bash -l\n')
        batchfile.write('#SBATCH --error="'+cwd+'/qsub/out/'+today+'/slurm-'+runname+'.err"\n')
        batchfile.write('#SBATCH --output="'+cwd+'/qsub/out/'+today+'/slurm-'+runname+'.out"\n')
        batchfile.write('#SBATCH -J stop-prod\n')
        batchfile.write('#SBATCH --partition=short\n')
        batchfile.write('#SBATCH --nodes=1\n')
        batchfile.write('#SBATCH --ntasks-per-node=32\n')
        batchfile.write('#SBATCH --time=04:00:00\n')
        batchfile.write('cd '+dirname+'\n')
        batchfile.write('mkdir OutHepMC-stop_bino_mssm_LO_CT14_tst'+hand.lower()+'\n')
        batchfile.write('module load gcc/8\n')
        batchfile.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
        batchfile.write('export PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/bin/:$PATH\n')
        batchfile.write(' '.join(['srun', 'Sherpa'] + sherpaopts)+' -e 10\n')
        batchfile.close()

#        os.system('sbatch '+batchname)

        while True:
          njobs = int(os.popen("squeue -u lscyboz -n 'stop-prod' | wc -l").read())
          print 'njobs = ', njobs-1, ' ...'
          if (njobs < 20): break
          time.sleep(20)

  inputlist.close()

if mode == "resubmit":

  if os.path.exists('qsub/resub-'+resubmitFolder): exit()
  os.system('mkdir qsub/resub-'+resubmitFolder+'\n')
  os.system('mkdir qsub/out/resub-'+resubmitFolder+'\n')

  inputlist = open("ParamCards/tL.half.list", 'r')
  for pcard in inputlist:
    pcard = pcard.split('\n')[0]
#  for pcard in glob("ParamCards/tL*.dat"):
    for hand in ['L', 'R']:

#        paramcard = pcard.split('/')[1]
        paramcard = pcard.replace('tL', 't'+hand)
        mstop = paramcard.split('.')[1].split('-')[0]
        mchi  = paramcard.split('-')[1].split('M1.')[1]

        print " ------- mstop = ", mstop, " - mchi = ", mchi, " -- ", hand, " ------- "

        runname = hand+'-Mst.'+mstop+'-M1.'+mchi+'-mt.173.1'
        dirname = cwd+'/resub-'+resubmitFolder+'/'+runname
        filename = 'RunStops_MSSM_CT14_tst'+hand.lower()+'-merged.dat'

	if not os.path.exists('qsub/out/'+resubmitFolder+'/slurm-'+runname+'.err'):
	  print "Output files not found!"
        if os.stat('qsub/out/'+resubmitFolder+'/slurm-'+runname+'.err').st_size == 0: continue

        os.system('mkdir -p '+dirname)
        os.system('cp '+'ParamCards/'+paramcard+' '+dirname)
        os.system('cp base/'+filename+' '+dirname)

        sherpaopts = [ "PRtgerr:="+str(err),
                       "PRmss:="+str(float(mstop)),
                       "PRcim:="+str(float(mchi)),
                       "PRshower:="+ps,
                       "PRparamcard:=t"+hand+"."+mstop+"-M1."+mchi+"-mt.173.1.dat",
                       "-f ", filename ]

        batchname = 'qsub/resub-'+resubmitFolder+'/submit-'+runname+'.sh'
        batchfile = open(batchname, 'w')
        batchfile.write('#!/bin/bash -l\n')
        batchfile.write('#SBATCH --error="'+cwd+'/qsub/out/resub-'+resubmitFolder+'/slurm-'+runname+'.err"\n')
        batchfile.write('#SBATCH --output="'+cwd+'/qsub/out/resub-'+resubmitFolder+'/slurm-'+runname+'.out"\n')
        batchfile.write('#SBATCH -J stop-prod\n')
        batchfile.write('#SBATCH --partition=short\n')
        batchfile.write('#SBATCH --nodes=1\n')
        batchfile.write('#SBATCH --ntasks-per-node=32\n')
        batchfile.write('#SBATCH --time=04:00:00\n')
        batchfile.write('cd '+dirname+'\n')
        batchfile.write('mkdir OutHepMC-stop_bino_mssm_LO_CT14_tst'+hand.lower()+'\n')
        batchfile.write('module load gcc/8\n')
        batchfile.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
        batchfile.write('export PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/bin/:$PATH\n')
        batchfile.write(' '.join(['srun', 'Sherpa'] + sherpaopts)+' -e 10\n')
        batchfile.close()

        os.system('sbatch '+batchname)

        while True:
          njobs = int(os.popen("squeue -u lscyboz -n 'stop-prod' | wc -l").read())
          print 'njobs = ', njobs-1, ' ...'
          if (njobs < 20): break
          time.sleep(20)

  inputlist.close()

if mode == "grid-unitywidth":

  if os.path.exists('qsub/resub-unitywidth-'+resubmitFolder): exit()
  os.system('mkdir qsub/resub-unitywidth-'+resubmitFolder+'\n')
  os.system('mkdir qsub/out/resub-unitywidth-'+resubmitFolder+'\n')
  os.system('mkdir ' + cwd + '/resub-unitywidth-' + resubmitFolder )

  inputlist = open("ParamCards/tL.half.list", 'r')
  for pcard in inputlist:
    pcard = pcard.split('\n')[0]
#  for pcard in glob("ParamCards/tL*.dat"):
    for hand in ['L', 'R']:

#        paramcard = pcard.split('/')[1]
        paramcard = pcard.replace('tL', 't'+hand)
        mstop = paramcard.split('.')[1].split('-')[0]
        mchi  = paramcard.split('-')[1].split('M1.')[1]

        print " ------- mstop = ", mstop, " - mchi = ", mchi, " -- ", hand, " ------- "

        runname = hand+'-Mst.'+mstop+'-M1.'+mchi+'-mt.173.1'
        dirname = cwd+'/resub-unitywidth-'+resubmitFolder+'/'+runname
        filename = 'RunStops_MSSM_CT14_tst'+hand.lower()+'-merged.dat'

        badfile = open(cwd+'/'+resubmitFolder+'/'+runname+"/status.xs", 'r')
        if (badfile.readline() == "GOOD"): 
          badfile.close()
          continue;
        badfile.close()

        os.system("rsync -azv --exclude Out* --exclude Status* --exclude Results* " + cwd +'/' + resubmitFolder + '/' + runname + " " + cwd+'/resub-unitywidth-'+resubmitFolder)
        time.sleep(5)
        os.system("rm -rf " + dirname + "/Results*.db*")
        os.system("sed -i 's/^DECAY 1000006.*$/DECAY 1000006 1.000000e+00/g' " + dirname + '/' + paramcard)

        sherpaopts = [ "PRtgerr:=0.01",
                       "PRmss:="+str(float(mstop)),
                       "PRcim:="+str(float(mchi)),
                       "PRshower:="+ps,
                       "PRparamcard:=t"+hand+"."+mstop+"-M1."+mchi+"-mt.173.1.dat",
                       "-f ", filename ]

        batchname = 'qsub/resub-unitywidth-'+resubmitFolder+'/submit-'+runname+'.sh'
        batchfile = open(batchname, 'w')
        batchfile.write('#!/bin/bash -l\n')
        batchfile.write('#SBATCH --error="'+cwd+'/qsub/out/resub-unitywidth-'+resubmitFolder+'/slurm-'+runname+'.err"\n')
        batchfile.write('#SBATCH --output="'+cwd+'/qsub/out/resub-unitywidth-'+resubmitFolder+'/slurm-'+runname+'.out"\n')
        batchfile.write('#SBATCH -J stop-prod\n')
        batchfile.write('#SBATCH --partition=short\n')
        batchfile.write('#SBATCH --nodes=1\n')
        batchfile.write('#SBATCH --ntasks-per-node=32\n')
        batchfile.write('#SBATCH --time=04:00:00\n')
        batchfile.write('cd '+dirname+'\n')
        batchfile.write('mkdir OutHepMC-stop_bino_mssm_LO_CT14_tst'+hand.lower()+'\n')
        batchfile.write('module load gcc/8\n')
        batchfile.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
        batchfile.write('export PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/bin/:$PATH\n')
        batchfile.write(' '.join(['srun', 'Sherpa'] + sherpaopts)+' -e 10\n')
        batchfile.close()

        os.system('sbatch '+batchname)

        while True:
          njobs = int(os.popen("squeue -u lscyboz -n 'stop-prod' | wc -l").read())
          print 'njobs = ', njobs-1, ' ...'
          if (njobs < 20): break
          time.sleep(20)

  inputlist.close()

if mode == "event":

  os.system('mkdir qsub/event-'+resubmitFolder+'\n')
  os.system('mkdir qsub/out/event-'+resubmitFolder+'\n')

  inputlist = open("ParamCards/tL.half.list", 'r')
  for pcard in inputlist:
    pcard = pcard.split('\n')[0]
#  for pcard in glob("ParamCards/tL*.dat"):
    for hand in ['L', 'R']:

#        paramcard = pcard.split('/')[1]
        paramcard = pcard.replace('tL', 't'+hand)
        mstop = paramcard.split('.')[1].split('-')[0]
        mchi  = paramcard.split('-')[1].split('M1.')[1]

        print " ------- mstop = ", mstop, " - mchi = ", mchi, " -- ", hand, " ------- "

        runname = hand+'-Mst.'+mstop+'-M1.'+mchi+'-mt.173.1'
        dirname = cwd+'/'+resubmitFolder+'/'+runname
        filename = 'RunStops_MSSM_CT14_tst'+hand.lower()+'-merged.dat'

        statusf = open(dirname+"/status.xs",'r')
        if statusf.readline() == "FAILED":
          print "No final integrated xs"
          continue
        statusf.close()

        sherpaopts = [ "PRtgerr:=999",
                       "PRmss:="+str(float(mstop)),
                       "PRcim:="+str(float(mchi)),
                       "PRshower:="+ps,
                       "PRparamcard:=t"+hand+"."+mstop+"-M1."+mchi+"-mt.173.1.dat",
                       "-f ", filename ]

        batchname = 'qsub/event-'+resubmitFolder+'/submit-'+runname+'.sh'
        batchfile = open(batchname, 'w')
        batchfile.write('#!/bin/bash -l\n')
        batchfile.write('#SBATCH --error="'+cwd+'/qsub/out/event-'+resubmitFolder+'/slurm-'+runname+'.err"\n')
        batchfile.write('#SBATCH --output="'+cwd+'/qsub/out/event-'+resubmitFolder+'/slurm-'+runname+'.out"\n')
        batchfile.write('#SBATCH -J prod-st\n')
        batchfile.write('#SBATCH --partition=short\n')
        batchfile.write('#SBATCH --nodes=1\n')
        batchfile.write('#SBATCH --ntasks-per-node=32\n')
        batchfile.write('#SBATCH --time=04:00:00\n')
        batchfile.write('cd '+dirname+'\n')
        batchfile.write('mkdir OutHepMC-stop_bino_mssm_LO_CT14_tst'+hand.lower()+'\n')
        batchfile.write('module load gcc/8\n')
	batchfile.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
        batchfile.write('export PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/bin/:$PATH\n')
        batchfile.write(' '.join(['srun', 'Sherpa', ' PRran:=1$SLURM_JOB_ID'] + sherpaopts)+' -e 50000\n')
        batchfile.close()

        os.system('sbatch '+batchname)

        while True:
          njobs = int(os.popen("squeue -u lscyboz -n 'prod-st' | wc -l").read())
          print 'njobs = ', njobs-1, ' ...'
          if (njobs < 20): break
          time.sleep(20)

  inputlist.close()

elif mode == "test":
    if submit == 'serial':
        with open( "sherpa"+fileidentifier+".log", "a+" ) as outfile:
            process = Popen( ["Sherpa"] + sherpaopts, stdout=outfile )
    elif submit == 'slurm':
     for i in range(n):
        process = Popen( [ "sbatch"], stdin = PIPE )
        process.stdin.write('#!/bin/bash -l\n')
#        process.stdin.write('#SBATCH --array='+str(nstart)+'-'+str(nstart+n-1)+'\n')
        process.stdin.write('#SBATCH --error="slurm-%A_%a.err"\n')
        process.stdin.write('#SBATCH --output="slurm-%A_%a.out"\n')
        process.stdin.write('#SBATCH --job-name="evt_'+fileName+'_'+str(i)+'"\n')
        process.stdin.write('#SBATCH --partition=express\n')
        process.stdin.write('#SBATCH --ntasks=1\n')
        process.stdin.write('#SBATCH --time=00:10:00\n')
	process.stdin.write('module load gcc/6.3\n')
#	process.stdin.write('source /ptmp/lscyboz/Rivet-2.5.4/rivetenv.sh\n')
	process.stdin.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
	process.stdin.write('export PATH=/ptmp/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/bin/:$PATH\n')
#	process.stdin.write('srun Sherpa -f Run_wwbb.dat -e 1000000 -A 06.04.17.\n')
        process.stdin.write(' '.join(['srun', 'Sherpa',' PRran:=1$SLURM_JOB_ID'] + sherpaopts)+'\n')
        process.stdin.write('exit 0\n')
        process.stdin.close()
    elif submit == 'condor':
        process = Popen( [ "condor_submit"], stdin = PIPE )
        process.stdin.write('#!/usr/bin/env condor_submit\n')
        process.stdin.write('Universe = vanilla\n')
        process.stdin.write('Executable = '+str(check_output(['which','Sherpa']))+'\n')
        process.stdin.write('Notification = Never\n')
        process.stdin.write('Input = /dev/null\n')
        process.stdin.write('Output = condor'+fileidentifier+'.$(Cluster).$(Process)-stdout.txt\n')
        process.stdin.write('Error = condor'+fileidentifier+'.$(Cluster).$(Process)-stderr.txt\n')
        process.stdin.write('Log = condor'+fileidentifier+'.$(Cluster).$(Process)-stdlog.txt\n')
        process.stdin.write('Getenv = true\n')
        process.stdin.write('Requirements = Arch == "X86_64" && Pool == "Theory" && Distro == "openSUSE" && (Machine != "pcl343.mppmu.mpg.de")\n')
        process.stdin.write('Initialdir = '+getcwd()+'\n')
        for ran in range(nstart, nstart + n ):
            process.stdin.write(' '.join(['Arguments  = "-o 0 random:=1'+str(ran)] + sherpaopts + ['"'])+'\n')
            process.stdin.write('Queue\n\n')
        process.stdin.close()

elif mode == "rivet":
  for f in glob("*.hepmc2g"):

    process = Popen( [ "sbatch"], stdin = PIPE )
    process.stdin.write('#!/bin/bash -l\n')
#        process.stdin.write('#SBATCH --array='+str(nstart)+'-'+str(nstart+n-1)+'\n')
    process.stdin.write('#SBATCH --error="slurm-%A_%a.err"\n')
    process.stdin.write('#SBATCH --output="slurm-%A_%a.out"\n')
    process.stdin.write('#SBATCH --job-name="stopcorr"\n')
    process.stdin.write('#SBATCH --partition=express\n')
    process.stdin.write('#SBATCH --time=00:30:00\n')
    process.stdin.write('module load gcc/8\n')
    process.stdin.write('source  /ptmp/lscyboz/Rivet-2.5.4/rivetenv.sh\n')
    process.stdin.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
    process.stdin.write(' '.join(['srun', 'rivet', f, '-o ' + f.split(".hepmc2g")[0] + '.yoda' ] + analyses)+'\n')
    process.stdin.write('exit 0\n')
    process.stdin.close()


