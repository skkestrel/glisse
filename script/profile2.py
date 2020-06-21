import shutil
import random
import subprocess
import sys
import math
import os

CONFIG = """
Initial-Time 0
Time-Step 1e-4
Final-Time {2}
Time-Block-Size {0}
Split-Track-File 0
Track-Interval 0
Log-Interval 128
Status-Interval 1
Resync-Interval {4}
CPU-Thread-Count 1
Output-Folder temp-data
Limit-Particle-Count {1}
Dump-Interval 25000
Input-File temp-state.in
Write-Binary-Output 0
Read-Binary-Input 0
"""

nparts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192, 9216, 10240, 11264, 12288, 13312, 14436, 15360, 16384]
nparts += range(16384, 16384 * 2, 2048)
nparts += range(16384 * 2, 16384 * 4, 2048)
nparts += range(16384 * 4, 133120, 2048)
nparts += [133120]

binary = "bin/glisse"

def randmul():
	return 1 + random.uniform(-1, 1) * 0.00000001

def genstate(npl, npa):
	with open('temp-state.in', 'w') as stateout:
		stateout.write("{0}\n".format(npl + 1))
		stateout.write("1\n0 0 0\n0 0 0\n0\n")
		for i in range(npl):
			r = i + 1
			v = math.sqrt(1. / (i + 1))
			stateout.write("0.0000001\n{0} {1} 0\n{2} {3} 0\n{4}\n"
					.format(r, 0, 0, v * randmul(), i + 1))
		stateout.write("{0}\n".format(npa))
		for i in range(npa):
			stateout.write("{0} 0 0\n0 {1} 0\n{2} 0 0\n"
					.format(npl + 1, randmul() * math.sqrt(1. / (npl + 1)), i + 1))

def run():
	shutil.rmtree('temp-data', True)
	subprocess.call([binary, "temp-config.in"])
	with open('temp-data/time.out', 'r') as outin:
		for line in outin:
			line = line.strip().split()
			if len(line) == 0: continue
			if line[0] == "time":
				maxtime = float(line[1]) * 60.
	return maxtime

genstate(4, 133120)
# number timesteps
nstep = 16384 * 8
for T in nparts * 3:
	with open('temp-config.in', 'w') as cfgout:
		cfgout.write(CONFIG.format(120, T, nstep * 1e-4, 1, 4))
	maxtime = run()



	mp_count = 13
	max_tpb = 1024

	threads_per_mp = (T + mp_count - 1) // mp_count
	allowed_blocks = (threads_per_mp + max_tpb - 1) // max_tpb
	threads_per_block = (threads_per_mp + allowed_blocks - 1) // allowed_blocks
	block_size = (threads_per_block + 31) // 32 * 32
	grid_size = (T + block_size - 1) // block_size

	with open('prof.csv', 'a') as profout:
		profout.write('{0},{1},{2},{3},{4},{5}\n'.format(T, int(math.ceil(nstep / 120) * 120), maxtime, block_size, grid_size, block_size * grid_size))
