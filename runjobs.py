from subprocess import Popen, PIPE #, check_output
from sys import argv, exit
from glob import glob
from datetime import date
import os

mode = argv[1] # grid or event
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

  for pcard in glob("ParamCards/tL*.dat"):
    for hand in ['L', 'R']:

        paramcard = pcard.split('/')[1]
        mstop = paramcard.split('.')[1].split('-')[0]
        mchi  = paramcard.split('-')[1].split('M1.')[1]

        print " ------- mstop = ", mstop, " - mchi = ", mchi, " -- ", hand, " ------- "

        runname = hand+'-Mst.'+mstop+'-M1.'+mchi+'-mt.173.1'
        dirname = cwd+'/'+today+'/'+runname
        filename = 'RunStops_MSSM_CT14_tst'+hand.lower()+'-merged.dat'

        os.system('mkdir -p '+dirname)
        os.system('cp '+pcard+' '+dirname)
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
        batchfile.write('mkdir OutHepMC-stop_bino_mssm_LO_CT14_tstl\n')
        batchfile.write('module load gcc/8\n')
        batchfile.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
        batchfile.write('export LD_LIBRARY_PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
        batchfile.write('export PATH=/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations/bin/:$PATH\n')
        batchfile.write(' '.join(['srun', 'Sherpa'] + sherpaopts)+' -e 10\n')
        batchfile.close()

        os.system('sbatch '+batchfile)

        while True:
          njobs = os.system("squeue -u lscyboz -n 'stop-prod' | wc -l")
          print 'njobs = ', njobs-1, ' ...'
          if (njobs < 20): break


elif mode == "event":
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
#        process.stdin.write('#SBATCH --nodes=1\n')
#        process.stdin.write('#SBATCH --ntasks-per-node=32\n')
        process.stdin.write('#SBATCH --time=00:30:00\n')
	process.stdin.write('module load gcc/8\n')
#	process.stdin.write('source /ptmp/lscyboz/Rivet-2.5.4/rivetenv.sh\n')
	process.stdin.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/LHAPDF-6.2.3/lib:$LD_LIBRARY_PATH\n')
	process.stdin.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
	process.stdin.write('export PATH=/ptmp/lscyboz/SHERPA-MC-2.2.8_stopCorrelations/bin/:$PATH\n')
        process.stdin.write('export LD_LIBRARY_PATH=/ptmp/lscyboz/SHERPA-MC-2.2.8_stopCorrelations/lib/SHERPA-MC/:$LD_LIBRARY_PATH\n')
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
	process.stdin.write('export PATH=/ptmp/lscyboz/SHERPA-MC-2.2.8_stopCorrelations/bin/:$PATH\n')
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
    ids = []
    idsr = []
    for f in glob("*tstl*.hepmc"):
#        ids += [ f.split('.hepmc')[0].split('/')[-1].split('_')[0] ]
        ids += [f]

#    ids = list(set(ids))
#    print ids
#    for g in ids:
#        ls = glob("Output_CT14_emu_SMEARATLAS/"+g+"_*.hepmc2g")
#        ls1, ls2 = split_list(ls)

    process = Popen( [ "sbatch"], stdin = PIPE )
    process.stdin.write('#!/bin/bash -l\n')
#        process.stdin.write('#SBATCH --array='+str(nstart)+'-'+str(nstart+n-1)+'\n')
    process.stdin.write('#SBATCH --error="slurm-%A_%a.err"\n')
    process.stdin.write('#SBATCH --output="slurm-%A_%a.out"\n')
    process.stdin.write('#SBATCH --job-name="stopcorr"\n')
    process.stdin.write('#SBATCH --partition=short\n')
    process.stdin.write('#SBATCH --time=04:00:00\n')
    process.stdin.write('#SBATCH --nodes=1\n')
    process.stdin.write('#SBATCH --ntasks-per-node=32\n')
    process.stdin.write('module load gcc/8\n')
    process.stdin.write('source  /ptmp/lscyboz/Rivet-2.5.4/rivetenv.sh\n')
    process.stdin.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
    process.stdin.write(' '.join(['srun', 'rivet', '-o Output_yodas-new/stop_bino-L-noshower.yoda' ] + analyses + ids)+'\n')
    process.stdin.write('exit 0\n')
    process.stdin.close()


    for i in glob("*tstr*.hepmc"):
        idsr += [i]

    process = Popen( [ "sbatch"], stdin = PIPE )
    process.stdin.write('#!/bin/bash -l\n')
#        process.stdin.write('#SBATCH --array='+str(nstart)+'-'+str(nstart+n-1)+'\n')
    process.stdin.write('#SBATCH --error="slurm-%A_%a.err"\n')
    process.stdin.write('#SBATCH --output="slurm-%A_%a.out"\n')
    process.stdin.write('#SBATCH --job-name="stopcorr"\n')
    process.stdin.write('#SBATCH --partition=short\n')
    process.stdin.write('#SBATCH --time=04:00:00\n')
    process.stdin.write('#SBATCH --nodes=1\n')
    process.stdin.write('#SBATCH --ntasks-per-node=32\n')
    process.stdin.write('module load gcc/8\n')
    process.stdin.write('source  /ptmp/lscyboz/Rivet-2.5.4/rivetenv.sh\n')
    process.stdin.write('export RIVET_ANALYSIS_PATH=/ptmp/lscyboz/RivetCustomJan_new:$RIVET_ANALYSIS_PATH\n')
    process.stdin.write(' '.join(['srun', 'rivet', '-o Output_yodas-new/stop_bino-R-noshower.yoda' ] + analyses + idsr)+'\n')
    process.stdin.write('exit 0\n')
    process.stdin.close()

