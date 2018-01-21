####
# @author: Núria Queralt Rosinach
# @date: 06-12-2017
# @version: v4 [q1_args_used]


from neo4j.v1 import GraphDatabase, basic_auth
import argparse
import sys,os


# Parse and validate command line arguments
parser = argparse.ArgumentParser(description='process user given parameters')
parser.add_argument("-pwdl", "--pathwayDegreeLimit", required = False, dest = "pwDegree", help = "maximum node degree for pathway type nodes (default = 50)", default = "50")
parser.add_argument("-phdl", "--phenotypeDegreeLimit", required = False, dest = "phDegree", help = "maximum node degree for phenotype type nodes (default = 20)", default = "20")
parser.add_argument("-o", "--output", required = True, dest = "output", help = "tab output format")

args = parser.parse_args()

if len(args.output) == 0:
   sys.stderr.write("Not output filename provided. Exiting.\n")
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
    # run query
    output = list()
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
            RETURN count(distinct path) AS paths, count(distinct pw) as pathways, count(distinct ds) as diseases, count( distinct source_ortho ) AS source_ortho, count( distinct other_ortho ) AS other_ortho
            ORDER BY source_ortho, other_ortho DESC
            """
            result = session.run(query)
            pair = {}
            pair['source'] = source
            pair['target'] = target
            for record in result:
                pair['paths'] = record['paths']
                pair['pathways'] = record['pathways']
                pair['diseases'] = record['diseases']
                pair['source_ortho'] = record['source_ortho']
                pair['other_ortho'] = record['other_ortho']
            output.append(pair)

    # create outdir and output file
    if not os.path.isdir('./out'): os.makedirs('./out')
    sys.path.insert(0,'.')
    with open('./out/{}.tsv'.format(args.output), 'w') as f:
        f.write('source\ttarget\tpaths\tpathways\tdiseases\tsource_orthologs\tother_orthologs\n')
        for pairwise in output:
            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pairwise['source'], pairwise['target'], pairwise['paths'], pairwise['pathways'], pairwise['diseases'], pairwise['source_ortho'], pairwise['other_ortho']))

sys.stderr.write("{} QUERIES completed.\n".format(len(output)))


