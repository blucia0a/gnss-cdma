#!/bin/python
#
# This scripts runs a vtune amplifier program x amount of times and finds the
# mean, std, max, and min of the CPU Times found
#
# Created by Jake Johnson on March 10, 2018
# Last Modified March 10, 2018
#
# This program does not use advanced-hotspots.  When using advanced-hotspots, the
# CPU time increases up to the elapsed time for some reason.  But advanced-hotspots
# does give the CPU cycle count
#
# Make sure to run this in a directory where multiple vtune results can be stored
#
# In Linux environments, to enable profiling, please set
# echo 0 > /proc/sys/kernel/yama/ptrace_scope
#
#

## Relevant dependencies
import subprocess
import numpy
import re

## Script Configuration
#--- Provide list of available executables to profile
exe_list = ["avx_si32_avx2", "nom_si32_avx2", "avx_fl32_avx2", "nom_fl32_avx2",
"avx_si32_avx512", "nom_si32_avx512", "avx_fl32_avx512", "nom_fl32_avx512",
"reg_standalone"]
#--- Specify the maximum number of times a given executable will run
num_exec = 50
#--- Specify output results path
work_path = "-app-working-dir /home/gnssi9/dami7269/gnss-intrinsics/install/"
#--- Specify executable path
exe_path = "/home/gnssi9/dami7269/gnss-intrinsics/install/"
#--- Specify Vtune Amplifier command line executable
vtune_amp_path = "/opt/intel/vtune_amplifier_2018.1.0.535340/bin64/amplxe-cl"
#--- Profiler configurations
prof_conf = """ -collect hotspots -quiet -knob sampling-interval=1 """
#--- Store profiling times
cpu_times = []

## Script Execution
f= open("guru99.txt","w+")
for j in range(0,len(exe_list)):
    print("Running: ")
    print(exe_list[j])
    print("\n")
    # Reset CPU Times to store only latest results
    cpu_times = []
    for i in range(0,num_exec):
        command = vtune_amp_path + prof_conf + work_path + " -- " + exe_path + exe_list[j]
        output = subprocess.check_output(command, shell=True)

        # Go through output add desired values to array
        for line in output.splitlines():
            if re.search('CPU Time:', line):
                cpu_times.append(float(line.split()[2]))
                print('Analysis ' + str(i) + ' complete: CPU Time: ' + line.split()[2] + '\n')

    mean_cpu_times = numpy.mean(cpu_times)
    std_cpu_times = numpy.std(cpu_times)
    maximum_cpu_times = numpy.amax(cpu_times)
    minimum_cpu_times = numpy.amin(cpu_times)

    print("cpu_times: ")
    print(cpu_times)
    print("mean_cpu_times: " + str(mean_cpu_times))
    print("std_cpu_times: " + str(std_cpu_times))
    print("maximum_cpu_times: " + str(maximum_cpu_times))
    print("minimum_cpu_times: " + str(minimum_cpu_times))

    f.write(exe_list[j])
    f.write(str(cpu_times))
    f.write("\n\n")

f.close()
