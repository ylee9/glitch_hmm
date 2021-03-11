#!/bin/bash
#SBATCH --job-name=matlab_glitch
#SBATCH --out=matlab_out.dat
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=5000
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=skylake

module load gcc/7.3.0
module load python/3.6.4
module load openmpi/3.0.1
module load matlab/2017b
source /fred/oz022/ldunn_venvs/glitches/bin/activate

#export DATA_DIR=~/NCIdy4_scratch/ldunn_glitch/sofia_data_20180726_1e-8Hz_1e-18/
export DATA_DIR=/fred/oz022/ldunn_glitch/ldunn_final_det_alarm_rocs_paper_params_lower_noise/dfp
export MATLAB_CODE_DIR=/fred/oz022/ldunn_glitch/matlab_code_logs

cd $DATA_DIR
mkdir results_data
cp $MATLAB_CODE_DIR/process_single_uuid.m results_data
export OMP_NUM_THREADS=1
mpirun python3 $MATLAB_CODE_DIR/run_matlab.py $DATA_DIR $MATLAB_CODE_DIR
