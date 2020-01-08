#!/bin/bash -l
#SBATCH --error="/draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/qsub/out/slurm-R-220.0-80.0.err"
#SBATCH --output="/draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/qsub/out/slurm-R-220.0-80.0.out"
#SBATCH -J stop-xsecs
#SBATCH --partition=express
#SBATCH --time=00:30:00
cd /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/../PROD-Dec-17-2019/R-Mst.220-M1.80-mt.173.1
source /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/setup.sh
srun /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/calcxsec Results.stop_bino_mssm_LO_CT14_tstr-220.0-80.0 -f RunStops_MSSM_CT14_tstr-merged.dat PRparamcard:=tR.220-M1.80-mt.173.1.dat PRmss:=220.0 PRcim:=80.0 PRshower:=CSS PRtgerr:=999
