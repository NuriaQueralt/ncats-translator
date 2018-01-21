####
# @author: Núria Queralt Rosinach
# @date: 01-02-2018
# @version: v1
# @usage: workflow manager

import subprocess

query = "q1"

# run part 1: retrieve seed pairwise counts summary table for arguments combination
cmd =  "python3 {}_driver.py".format(query)
print(cmd)
subprocess.call(cmd, shell = True)

# run part 2: retrieve extended seed pairwise count summary table and the detailed paths for the selected arguments combination
cmd =  "python3 {}_1_driver.py".format(query)
print(cmd)
subprocess.call(cmd, shell = True)

cmd =  "python3 {}_2_driver.py".format(query)
print(cmd)
subprocess.call(cmd, shell = True)