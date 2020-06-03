import sys

'''
Script takes as input (1) ensembl id --> gene id table 
(2) DE genes output table from edgeR analysis

USAGE: "python hnscc-DE-genes-ensembl.py /path/to/ensembl_id_table /path/to/de_genes_table"
For each ensembl id listed in the DE genes table, script prints
gene id alias and gene name/function & sources along with FC and
p-values as written in the DE genes table
'''

ia = list()
#open ia genes 
with open(sys.argv[1], 'r') as f:
    #skip header line
    f.readline()
    # for each line
    for line in f:
        # split line into fields delimited by tabls
        line = line.rstrip().split('\t')
        # capture field indices with info of interest (ensembl id, geneid,
        # gene desc)
        ia.append(line[0])

# open CPM genes table
with open(sys.argv[2], 'r') as f1:
    with open("phenotype.cls", 'w') as output:
        output.write("#numeric"+"\n")
        # skip 2 header lines
        targets = f1.readline()
        targets = targets.rstrip().split("\t")
        if "Normal" in targets:
            print("\nERROR: YOU FORGOT TO REMOVE NORMALS!!!! Please remove normals and re-run")
            exit(1)
        f1.readline()
        # for each line in table
        for line in f1:
            # split line on space delimiter
            line = line.rstrip().split("\t")
            # strip version # from ensembl id
            cpm_ensembl_id = line[0]
            # if ensembl id exists in ia list
            if cpm_ensembl_id in ia:
                # print line with gene id and desc and cpms
                output.write("#"+line[1]+"\n")
                output.write("\t".join(line[3:])+"\n")
        
    
