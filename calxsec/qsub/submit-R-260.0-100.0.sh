#!/bin/bash -l
#SBATCH --error="/draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/qsub/out/slurm-R-260.0-100.0.err"
#SBATCH --output="/draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/qsub/out/slurm-R-260.0-100.0.out"
#SBATCH -J stop-xsecs
#SBATCH --partition=express
#SBATCH --time=00:30:00
cd /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/../PROD-Dec-17-2019/R-Mst.260-M1.100-mt.173.1
source /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/setup.sh
srun /draco/u/lscyboz/SHERPA-MC-2.2.8_stopCorrelations-new/PROD/calxsec/calcxsec Results.stop_bino_mssm_LO_CT14_tstr-260.0-100.0 -f RunStops_MSSM_CT14_tstr-merged.dat PRparamcard:=tR.260-M1.100-mt.173.1.dat PRmss:=260.0 PRcim:=100.0 PRshower:=CSS PRtgerr:=999
