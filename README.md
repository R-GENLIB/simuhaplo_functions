# simuhaplo_scripts
standalone scripts to be used with the output of gen.simuhaplo 

'percent_IBD.py' can be used to find the percent shared IBD between all pairs of probands in a 'probands_haplotypes.txt' output file. To use this script:

'''
>python3 ./percent_IBD.py ./probands_haplotypes.txt
'''
the script will create an output file in the same directory called 'IBD_pairs.csv' that has details for IBD sharing for every pair of proband, for every simulation in the 'probands_haplotypes.txt' file. 

The 'traceback.py' script can be used to traceback a segment in a proband up it's path of inheritance. To use it you need to pass in the genealogy file, the 'All_nodes_haplotypes.txt' output file, the 'probands_haplotypes.txt' output file, the ID of the proband you want to trace back from, and the ID of the founder of the segment being traced back,
