#!/usr/bin/env python

# import necessary packages

# write .job files for each core
for i in range(1,17):
	with open('t{}.job'.format(i),'w') as f:
		f.write('#!/bin/bash\n')
		f.write('#SBATCH -A SMU-Math-6370 \t{}\n'.format('# account name'))
		f.write('#SBATCH -J driver_omp_t{} \t{}\n'.format(i,'# job name'))
		f.write('#SBATCH -o t{}_out.%j \t\t{}\n'.format(i,'# output file'))
		f.write('#SBATCH -e t{}_err.%j \t\t{}\n'.format(i,'# error file'))
		f.write('#SBATCH -N 1 \t\t\t{}\n'.format('# total nodes requested'))
		f.write('#SBATCH -n 1 \t\t\t{}\n'.format('# total MPI tasks requested'))
		f.write('#SBATCH -p normal \t\t{}\n'.format('# queue name'))
		f.write('#SBATCH -t 00:10:00 \t\t{}\n\n'.format('# total time requested <hh:mm:ss>'))
		f.write('export OMP_NUM_THREADS={}'.format(i))
		f.write('\nibrun ./WENOConvergence.exe 300 300 S \n')

