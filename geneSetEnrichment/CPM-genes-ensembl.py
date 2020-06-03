import sys

'''
Script takes as input (1) ensembl id --> gene id table 
(2) DE genes output table from edgeR analysis

USAGE: "python hnscc-DE-genes-ensembl.py /path/to/ensembl_id_table /path/to/de_genes_table"
For each ensembl id listed in the DE genes table, script prints
gene id alias and gene name/function & sources along with FC and
p-values as written in the DE genes table
'''

ensembl = {}
#open ensembl id table 
with open(sys.argv[1], 'r') as f:
    #skip header line
    f.readline()
    # for each line
    for line in f:
        # split line into fields delimited by tabls
        line = line.rstrip().split('\t')
        # capture field indices with info of interest (ensembl id, geneid,
        # gene desc)
        ensembl_id = line[0]
        name = line[1]
        desc = ''
        # desc only available if line contained 4+ fields
        if len(line) >= 4:
            desc = line[3]
        # assign ensembl id & name/desc as k,v pair in ensembl dictionary
        ensembl[ensembl_id] = [name, desc]

files = list()
group = list()
# open Targets.txt
with open(sys.argv[4], 'r') as t:
    #skip header
    t.readline()
    # for each line in table
    for line in t:
        line = line.rstrip().split()
        files.append(line[0]) # get all file names
        group.append(line[1]) # get corresponding groups

counter = 0
# open CPM genes table
with open(sys.argv[2], 'r') as f1:
    with open(sys.argv[3], 'w') as output:
        output.write('\t'.join(["", "", ""]+group[0:])+"\n")
        output.write('\t'.join(["ensembl_id", "gene_name", "gene_description"]+files[0:])+"\n")
        # skip header
        f1.readline()
        # for each line in table
        for line in f1:
            # split line on space delimiter
            line = line.rstrip().split()
            # strip version # from ensembl id
            cpm_ensembl_id = line[0].split('.')[0]
            # if ensembl id exists in ensembl dictionary
            if cpm_ensembl_id in ensembl:
                # print line with gene id and desc and cpms
                output.write('\t'.join([cpm_ensembl_id, ensembl[cpm_ensembl_id][0], ensembl[cpm_ensembl_id][1]]+line[1:]))
                output.write('\n')
            # handle non-match case (not reached)
            else:
                print (counter, end=" ")
                print (cpm_ensembl_id, end="")
                print ("not found in Ensembl db")
            counter = counter + 1
        
    
