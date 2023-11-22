"""
 * Copyright (C) 2022 Lehigh University.
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

# looks at input.txt and makes n_jobs temporary copies to be run with main with different output names
filepath = "input.txt"
n_jobs = 20000
max_num_parallel_jobs = 10

input_file_str = ""
org_output_name = ""
with open(filepath) as fp:
    line = fp.readline()
    
    while(line):
        split_string = line.split()
        
        param_name = split_string[0]
        
        if param_name == 'outputName':
            org_output_name = split_string[1]
        else:
            input_file_str += line
        
        
        line = fp.readline()

tmp_output_file_names = []
for i in range(0, n_jobs):
    tmp_output_name = filepath + '_tmp_' + str(i)
    tmp_output_file_names.append(tmp_output_name)
    
    tmp_output = open(tmp_output_name, 'w')
    
    tmp_output.write(input_file_str)
    
    tmp_output.write('{0} {1}\n'.format('outputName', org_output_name + '_' + str(i)))
    
    tmp_output.close()

count = 0

running_jobs = set()

#https://stackoverflow.com/questions/4992400/running-several-system-commands-in-parallel
curr_directory = os.getcwd()

with open(os.devnull, 'w') as devnull:
    for i in range(0, n_jobs):
        print("Running job {0}".format(tmp_output_file_names[i]))
        running_jobs.add(subprocess.Popen([curr_directory + "//main", tmp_output_file_names[i]], stdout=devnull))
        
        count += 1
        
        while len(running_jobs) >= max_num_parallel_jobs or (count == n_jobs and len(running_jobs) > 0):
            time.sleep(10)
            running_jobs.difference_update([p for p in running_jobs if p.poll() is not None])
        
for i in range(0, n_jobs):
    os.remove(tmp_output_file_names[i])
