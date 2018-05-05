#!/usr/bin/env python

# read in results data
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

times = []
with open('time_data.txt','r') as f:
	for line in f:
		times.append(float(line))

with open('results.txt','w') as f:
	 f.write('{:^8} {:^20} {:<20} {:^20}\n\n'.format('Threads','Execution time','Parallel Speedup','Parallel Efficiency'))
	 print('{:^8} {:^20} {:<20} {:^20}'.format('Threads','Execution time','Parallel Speedup','Parallel Efficiency'))
	 for p,t in enumerate(times,start=1):
	 	f.write('{:^8} {:^20.4f} {:<20.16f} {:<20.16f}\n'.format(p,t,times[0]/t,times[0]/t/p))
	 	print('{:^8} {:^20.4f} {:<20.16f} {:<20.16f}'.format(p,t,times[0]/t,times[0]/t/p))



plt.loglog(np.array([ _ for _ in range(1,17)]),np.array(times))
plt.xlabel('Processors')
plt.ylabel('Time')
plt.title('Strong Scaling Plot')
savefig('Stongscaling.pdf')
plt.show()
