####
# @author: Núria Queralt Rosinach
# @date: 06-12-2017
# @version: v1 [q1_paths]


from neo4j.v1 import GraphDatabase, basic_auth
import argparse
import sys,os
import json
import yaml


# Core function to parse neo4j results
def parsePath( path ):
    out = {}
    out['Nodes'] = []
    for node in path['path'].nodes:
        n = {}
        n['idx'] = node.id
        n['label'] = list(node.labels)[0]
        n['id'] = node.properties['id']
        n['preflabel'] = node.properties['preflabel']
        n['description'] = node.properties['description']
        out['Nodes'].append(n)
    out['Edges'] = []
    for edge in path['path'].relationships:
        e = {}
        e['idx'] = edge.id
        e['start_node'] = edge.start
        e['end_node'] = edge.end
        e['type'] = edge.type
        e['preflabel'] = edge.properties['property_label']
        e['references'] = edge.properties['reference_uri']
        out['Edges'].append(e)
    return out

# Parse and validate command line arguments
parser = argparse.ArgumentParser(description='process user given parameters')
parser.add_argument("-f", "--format", required = False, dest = "format", help = "yaml/json", default="json")
parser.add_argument("-pwdl", "--pathwayDegreeLimit", required = False, dest = "pwDegree", help = "maximum node degree for pathway type nodes (default = 50)", default = "50")
parser.add_argument("-phdl", "--phenotypeDegreeLimit", required = False, dest = "phDegree", help = "maximum node degree for phenotype type nodes (default = 20)", default = "20")
parser.add_argument("-o", "--output", required = True, dest = "output", help = "tab output format")

args = parser.parse_args()

if(args.format not in ["json","yaml"]):
   sys.stderr.write("Invalid format. Exiting.\n")
   exit()

# initialize neo4j
driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "xena"))

# seed nodes to pairwise
seed = list( [
    #'OMIM:615273', # NGLY1 deficiency
    #'CHEBI:506227', # GlcNAc
    'NCBIGene:55768', # NGLY1 human gene
    'NCBIGene:358', # AQP1 human gene
    'NCBIGene:11826', # AQP1 mouse gene
    'NCBIGene:4779', # NRF1 human gene* Ginger: known as NFE2L1. http://biogps.org/#goto=genereport&id=4779
    'NCBIGene:64772', # ENGASE human gene
    'NCBIGene:360', # AQP3 human gene
    'NCBIGene:282679' # AQP11 human gene
 ] )

# query topology
query_topology = """
(source:GENE)-[:`RO:HOM0000020`]-(:GENE)--(ds:DISO)--(:GENE)-[:`RO:HOM0000020`]-(g1:GENE)--(pw:PHYS)--(target:GENE)
"""

# ask the driver object for a new session
with driver.session() as session:
    # create outdir
    if not os.path.isdir('./out'): os.makedirs('./out')
    sys.path.insert(0,'.')

    # run query
    outputAll = list()
    for gene1 in seed:
        for gene2 in seed:
            if gene2 == gene1:
                continue
            source=gene1
            target=gene2
            query = """
            MATCH path=""" + query_topology + """
            MATCH (source { id: '""" + source + """'}), (target { id: '""" + target + """'})
            // no loops or only one pass per node
            WHERE ALL(x IN nodes(path) WHERE single(y IN nodes(path) WHERE y = x))
            WITH g1, ds, pw, path,
                 // count animal models
                 size( (source)-[:`RO:HOM0000020`]-() ) AS source_ortho,
                 size( (g1)-[:`RO:HOM0000020`]-() ) AS other_ortho,
                 // count node degree
                 max(size( (pw)-[]-() )) AS pwDegree,
                 max(size( (ds)-[]-() )) AS dsDegree,
                 // mark general nodes to filter path that contain them out
                 [n IN nodes(path) WHERE n.preflabel IN ['cytoplasm','cytosol','nucleus','metabolism','membrane','protein binding','visible','viable','phenotype']] AS nodes_marked,
                 // mark promiscuous edges to filter path that contain them out,
                 [r IN relationships(path) WHERE r.property_label IN ['interacts with','in paralogy relationship with','in orthology relationship with','colocalizes with']] AS edges_marked
            // condition to filter paths that do contain marked nodes and edges out
            WHERE size(nodes_marked) = 0 AND size(edges_marked) = 0 AND pwDegree <= """ + args.pwDegree + """ AND dsDegree <= """ + args.phDegree + """
            RETURN path
            """
            result = session.run(query)
            pair = {}
            pair['source'] = source
            pair['target'] = target
            # parse query results
            output = list()
            counter = 0
            for record in result:
                path_dct = parsePath(record)
                output.append(path_dct)
                counter += 1
                if (counter % 100000 == 0):
                    sys.stderr.write("Processed " + str(counter) + "\n")
            pair['paths'] = output
            outputAll.append(pair)

            # print output
            if(args.format == "yaml"):
                print(yaml.dump(outputAll, default_flow_style=False))
            elif(args.format == "json"):
                with open('./out/{}.json'.format(args.output), 'w') as f:
                    json.dump(outputAll, f)
                #print(json.dumps(outputAll))
            else:
                sys.stderr.write("Error.\n")


sys.stderr.write("{} QUERIES completed.\n".format(len(output)))


