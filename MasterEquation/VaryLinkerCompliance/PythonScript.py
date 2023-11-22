"""
 * Copyright (C) 2023 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
 """

import subprocess
import time
import os

filepath = "input.txt"
args = []
with open(filepath) as fp:
    line = fp.readline()
    
    while(line):
        split_string = line.split()
        
        args.append([])
        
        for s in split_string:
            args[len(args)-1].append(s)
        
        line = fp.readline()

max_num_parallel_jobs = 8

count = 0

running_jobs = set()

#https://stackoverflow.com/questions/4992400/running-several-system-commands-in-parallel
curr_directory = os.getcwd()

for a in args:
    #running_jobs.add(subprocess.Popen(["D:\\dmr518\\Martin\\3d-Simulations\\MembraneKMC-3d-Sphere-6-EndoToPoint.exe", a[0], a[1], a[2], a[3]]))]
    running_jobs.add(subprocess.Popen([curr_directory + "//MasterEqMethod", a[0]]))
    
    count += 1
    
    while len(running_jobs) >= max_num_parallel_jobs or (count == len(args) and len(running_jobs) > 0):
        time.sleep(10)
        # p.poll() returns None if process is running, so not None will be a job that has finished
        running_jobs.difference_update([p for p in running_jobs if p.poll() is not None])