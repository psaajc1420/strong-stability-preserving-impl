#/bin/sh

# script to run all .job files
for i in {1..8}
do 
   sbatch t$i.job		
done

read -p "Press enter to continue"

for i in {9..16}
do 
   sbatch t$i.job
done


