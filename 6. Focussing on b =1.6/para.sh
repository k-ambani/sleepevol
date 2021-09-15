#!/bin/sh
#SBATCH --job-name=sample_array
#Name of the partition/queue
#SBATCH --partition=gpu
## how much time is requested [ upper bound ]
#SBATCH --time=48:00:00
#specify the number of processors needed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#specify memory needed (MB)
#SBATCH --mem=20000
#SBATCH --array=38


##    notify me about this job
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=krutika.ambani@monash.edu
##    specify output file
###SBATCH --parajob.out

module load  matlab/r2019a
ulimit -s 60000


###export PARAMS =`head -n $SLURM_ARRAY_TASK_ID list | tail -1`
echo "SLURM_SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
### on this next line is where you actually run your program
./run_Evol_code.sh /usr/local/matlab/r2019a $SLURM_ARRAY_TASK_ID

### change "/usr/local/matlab/r2019a" to the location of the matlab that you
### are using... simply replace 'r2019a' with the version of the
### matlab module use are using
