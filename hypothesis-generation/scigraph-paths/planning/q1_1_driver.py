####
# @author: Núria Queralt Rosinach
# @date: 12-12-2017
# @version: v1
# @usage: read selected args combo file from data dir and run q1_1_cypher_to_hypotheses.py to get paths

import pandas as pd
import subprocess

# query
query = "q1"

# read file query rows
df = pd.read_csv('./data/{}_1.tsv'.format(query), header = 0, sep = '\t', keep_default_na = False)
for index, row in df.iterrows():
    # control execution flow
    if ( row.max_pathway_degree == 0 and row.max_disease_degree == 0 ):
        output_filename = "{}_1_in".format(query) + str(index) + "_pwdl50_phdl20_paths"
        cmd = "python3 {}_1_cypher_to_hypotheses.py -o ".format(query) + output_filename
    elif ( row.max_pathway_degree == 0):
        output_filename = "{}_1_in".format(query) + str(index) + "_pwdl50_phdl" + str(row.max_disease_degree) + "_paths"
        cmd = "python3 {}_1_cypher_to_hypotheses.py -phdl ".format(query) + str(row.max_disease_degree) + " -o " + output_filename
    elif ( row.max_disease_degree == 0 ):
        output_filename = "{}_1_in".format(query) + str(index) + "_pwdl" + str(row.max_pathway_degree) + "_phdl20_paths"
        cmd = "python3 {}_1_cypher_to_hypotheses.py -pwdl ".format(query) + str(row.max_pathway_degree) + " -o " + output_filename
    else:
        output_filename = "{}_1_in".format(query) + str(index) + "_pwdl" + str(row.max_pathway_degree) + "_phdl" + str(row.max_disease_degree) + "_paths"
        cmd = "python3 {}_1_cypher_to_hypotheses.py -pwdl ".format(query) + str(row.max_pathway_degree) + " -phdl " + str(row.max_disease_degree) + " -o " + output_filename

    # execute cypher hypothesis generation script
    print(cmd)
    subprocess.call(cmd, shell = True)