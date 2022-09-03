# simuhaplo_scripts
standalone scripts to be used with the output of gen.simuhaplo 

'percent_IBD.py' can be used to find the percent shared IBD between all pairs of probands in a 'probands_haplotypes.txt' output file. To use this script:


>python3 ./percent_IBD.py ./probands_haplotypes.txt

the script will create an output file in the same directory called 'IBD_pairs.csv' that has details for IBD sharing for every pair of proband, for every simulation in the 'probands_haplotypes.txt' file. 

The 'traceback.py' script can be used to traceback a segment in a proband up it's path of inheritance. To use it you need to pass in the genealogy file, the 'All_nodes_haplotypes.txt' output file, the 'probands_haplotypes.txt' output file, the ID of the proband you want to trace back from, and the ID of the founder from which the segment originates. The script will output two files: 'traceback.csv' and 'pathlist.txt'. For any simulation where the specified proband has a segment from the specified founder, it will be traced up the path of innheritance. 'Traceback.csv' contains a row for each such occurance, with the path and the length of the segment. 'pathlist.txt' contains the ID's of all individuals in a particular path. Note a value of '0' for the path means that the segment was not inherited down a single path, but instead comes from a combined path due to a homozygous recombination at an internal ancestor. 

E.g.: tracing back any instances of proband 1 inheriting a segment from individual 100:

>python3 ./traceback.py ./genealogy_file.csv ./all_nodes_haplotypes.txt ./proband_haplotypes.txt 1 100

The process_Genlib.pl script can be used to convert the haplotype representations in 'proband_haplotypes.txt' into sequence data. To use need to pass the proband_haplotypes file, a .map file, a .hap file, and the length of the simulated region in BP:

>process_Genlib.pl ./proband_haplotypes.txt ./hapfile.hap ./mapfile.map 100000000

the hapfile should contain one line for every founder. Each line should have the founder ID then the string representing the first haplotype, then a string representing the second haplotype. Seperated by whitespace. The string representing the haplotypes can be in any formatting, so long as all the haplotype strings are the same length for every individual. Each character in the string corresponds to a specific BP position, but any character can be used to represent the sequence information (only a single character, however). And the map file should contain the BP location for every 'marker' in the hap file. One line for every position. If you wish to represent missing data for some of the haplotypes you can use any character you wish, just so long as every sequence is the same length.
