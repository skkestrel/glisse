#!/usr/bin/env python3
# import numpy as np
# import matplotlib.pyplot as plt
import sys
import time
import os

CONFIG = """
Input-File state.in
Initial-Time 0
Time-Step 100
Final-Time 5000000
Time-Block-Size 1000
Output-Folder output
Log-Interval 10
Energy-Interval 10
Track-Interval 10
Resync-Interval 1
Dump-Interval 10
Keep-All-Dumps 1
Write-Barycentric-Track 0
Split-Track-File 0
Write-Split-Output 0
Write-Binary-Output 0
Read-Split-Input 0
Read-Binary-Input 0
Max-Kepler-Iterations 10
Big-G 1
Hill-Radius-Factor 1
Particle-Inner-Boundary 0.005
Particle-Outer-Boundary 1500
"""

os.chdir("benchmark")
os.system('rm timelog.out')

num = 100
numlist = []
while(num <= 200000):
    numlist.append(num)
    if(num < 1000):
        num += 100
    elif(num < 50000):
        num += 1000
    elif(num <= 200000):
        num += 10000
print(numlist)

for lps in numlist:
    ENDLINE = "Limit-Particle-Count " + str(lps)
    CONFIG_NEW = CONFIG + ENDLINE
    # print(CONFIG_NEW)

    os.system('rm -r output')
    os.system('rm config.in')
    fin = open("config.in", "w")
    fin.write(CONFIG_NEW)
    fin.close()

    now = time.time()
    os.system('../bin/glisse config.in')
    difference = time.time() - now
    print(difference)

    fout = open("timelog.out", "a")
    fout.write(str(lps) + " " + str(difference) + "\n")
    fout.close()