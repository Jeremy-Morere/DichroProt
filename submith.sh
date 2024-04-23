#!/bin/bash
#SBATCH -N 1
#SBATCH -p std
#SBATCH -J CD-amber
#SBARCH -n 120
#SBATCH -t 16-00:00:00

module purge
module load amber/16
module load gaussian/09/pgi

for rep in R-*/ ; do
   cd $rep
   echo "Go into $rep"
   sander -O -i mdgau-min -o test.out -c protein.rst -p 1pit_box.prmtop -r test.rst -x test.rst &
   cd ..
done

wait 

echo "done"
