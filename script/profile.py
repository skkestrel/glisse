import shutil
import subprocess
import sys
import math
import os

CONFIG = """
Initial-Time 0
Time-Step 1e-4
Final-Time {2}
Time-Block-Size {0}
Resolve-Encounters 0
Split-Track-File 0
Track-Interval 0
Log-Interval 32
Status-Interval 1
Enable-GPU 1
Output-Folder temp-data
Limit-Particle-Count {1}
Dump-Interval 25000
Input-File temp-state.in
Write-Binary-Output 0
Read-Binary-Input 0
"""

mode = int(sys.argv[1])

tbsizes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192, 9216, 10240, 11264, 12288, 13312, 14436, 15360, 16384]
nparts = [1, 10, 100, 500, 1000, 3000, 5000, 10000, 13312, 13313, 15000, 17500, 20000, 24000, 26626, 26627, 39936, 53248, 66560, 93184, 106496, 119808, 133120]
npls = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]

def genstate(npl, npa):
	with open('temp-state.in', 'w') as stateout:
		stateout.write("{0}\n".format(npl + 1))
		stateout.write("1\n0 0 0\n0 0 0\n0\n")
		for i in range(npl):
			r = i + 1
			v = math.sqrt(1. / (i + 1))
			stateout.write("0.00001\n{0} {1} 0\n{2} {3} 0\n{4}\n"
					.format(r, 0, 0, v, i + 1))
		stateout.write("{0}\n".format(npa))
		for i in range(npa):
			stateout.write("{0} 0 0\n{1} 0 0\n{2} 0 0\n"
					.format(npl + 1, math.sqrt(1. / (npl + 1)), i + 1))

def run():
	shutil.rmtree('temp-data', True)
	subprocess.call(["bin/sr", "temp-config.in"])
	with open('temp-data/time.out', 'r') as outin:
		for line in outin:
			line = line.strip().split()
			if len(line) == 0: continue
			if line[0] == "time":
				maxtime = float(line[1]) * 60.
	return maxtime

if mode == 0:
	genstate(4, 133120)
	nstep = 16384
	for T in reversed(tbsizes * 3):
		with open('temp-config.in', 'w') as cfgout:
			cfgout.write(CONFIG.format(T, 999999, nstep * 1e-4))
		maxtime = run()
		with open('timeblock-prof.csv', 'a') as profout:
			profout.write('{0},{1},{2},{3},{4}\n'.format(4, 133120, T, int(math.ceil(nstep / T) * T), maxtime))
elif mode == 1:
	genstate(4, 133120)
	nstep = 16384 * 8
	for T in nparts * 3:
		with open('temp-config.in', 'w') as cfgout:
			cfgout.write(CONFIG.format(16384, T, nstep * 1e-4))
		maxtime = run()
		with open('particle-prof.csv', 'a') as profout:
			profout.write('{0},{1},{2},{3},{4}\n'.format(4, T, 16384, int(math.ceil(nstep / 16384) * 16384), maxtime))
elif mode == 2:
	nstep = 16384 * 2
	for T in npls * 3:
		genstate(T, 133120)
		with open('temp-config.in', 'w') as cfgout:
			cfgout.write(CONFIG.format(16384, 133120, nstep * 1e-4))
		maxtime = run()
		with open('planet-prof.csv', 'a') as profout:
			profout.write('{0},{1},{2},{3},{4}\n'.format(T, 133120, 16384, int(math.ceil(nstep / 16384) * 16384), maxtime))