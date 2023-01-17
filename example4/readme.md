## Using reconstruct.pl

For this example we will use simulated haplotypes for a simple pedigree of 2 siblings (probands) and their grandparents (founders).  
![img](sample_gen.png)  
To use `reconstruct.pl` we need to provide the `Proband_haplotypes.txt` output from `gen.simuHaplo()`. The output file in this example is from a hypothetical 100,000,000 BP chromosomal segment, simulated with the following parameters:
```
require(GENLIB)
sample_data <- data.frame(ind=c(1,2,3,6,7,8,4,5), 
                          mother=c(3,3,4,7,0,0,0,0),
                          father=c(6,6,5,8,0,0,0,0),
                          sex =c(2,1,2,1,2,1,2,1))

sample_gen <- gen.genealogy(sample_data)
gen.simuHaplo(sample_gen, model = 1, simulNo = 10, model_params=c(1,1), cM_len=c(100,100), BP_len = 100000000)
```  
We also need to provide a file (e.g.: `founders_seq_data.hap`) that contains sequences for each founder chromosome. This file should contain a line for each founder chromosome, each line should be the ID followed by a string (encoding the genotype), in the same formatting as the example file. The sequences should all be strings of the same length, and can be composed of any characters, where each character represents the genotype at a particular BP location. The loci for each of these genotypes is contained in a secondary file (e.g.: `seq_data.map`). This file should only contain a sequence of integers, each on a new line, where each line is the BP position of the corresponding character in the strings contained in the other file.  
To use the script, pass the three required files in as arguments, and the BP length of the simulated chromosome/segment:  
 ```
 >reconstruct.pl Proband_haplotypes.txt founders_seq_data.hap seq_data.map 100000000
 ```
 The proband haplotypes will be converted to sequence data according to the founder haplotypes, and the output will be saved as `Probands.hap`