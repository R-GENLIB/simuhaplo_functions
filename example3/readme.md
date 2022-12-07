
## Using **`gen.simuHaplo_compare_IBD`**

This function is meant to analyze the output of `gen.simuHaplo`.
Specifically, given a pair of probands, it will return details on the
diploid IBD sharing between the pair’s haplotypes for every simulation.
Since `gen.simuHaplo` describes the haplotypes in terms of segments (BP
positions) along the chromosome (each segment labelled by the ID of the
founder chromosome it is inherited from), then haplotypes can be
assessed to be shared IBD if they have the same ancestral origin at the
same region of the chromosome. To identify IBD segments this function
looks for any overlaps between segments from the same unique founder
chromosome.

First: load GENLIB, and the genealogy:

``` r
library(GENLIB)
library(readr)
example_gen_table <- read_csv("../example_genealogy.csv")
```

``` r
gen_obj           <- gen.genealogy(example_gen_table)
gen_obj
```

    ##  GENLIB: Genealogical Data Software version  1 
    ## 
    ##   Number of individuals :  58942 
    ##   Number of parent-child relations :  99694 
    ##  Number of men  :  29342 
    ##   Number of women :  29600 
    ##   Number of subjects :  227 
    ## 
    ## 
    ##  Genealogical depth : 17
    ## 
    ##  Created on : Wed Dec 07 11:08:36 2022

Then run the simulation:

``` r
# we do not set all_nodes = 1, because we only need the proband haplotypes to compare
BP_len = 290000000
gen.simuHaplo(gen_obj ,pro=c(200,201), simulNo=200, model_params = c(1.98, 3.28), cM_len = c(198, 328), BP_len = BP_len)
```

    ## No map function specified to convert genetic distance to physical. Assumed constant along length of chromosome

    ## seed: 1670429316

    ## output file:  C:/Users/Mohan/OneDrive - University of Ottawa/simuhaplo_functions/example3/Proband_Haplotypes.txt

We use the same parameters as the other examples, but we only simulate
the haplotypes for 2 arbitrarily selected probands (ID: 200, 201). We do
not simulate the rest of the genealogy because it would be unnecessary
for this illustration. However, `gen.simuHaplo_compare_IBD` will still
work on output for more than 2 probands, it will just ignore all
probands other than the specified pair.  
To use `gen.simuHaplo_compare_IBD` we only need to pass in the ID’s for
our proband pair, the BP length of the chromosome, and the path to the
“Proband_Haplotypes.txt” file created from `gen.simuHaplo`. The function
will print to the R console the exact BP positions for each IBD segment
(for each haploid chromosome). This can be disabled by suppressing
messages in R.  
The function also returns a dataframe with the following columns:
**simulNo**, **len_Shared_IBD**, **num_seg**, **mean_len**

**simulNo** is the simulation number. The dataframe will only contain
rows for simulations where the pair have non-zero IBD sharing.  
**num_seg** is the total number of non-contiguous IBD segments. Divide
by 2 to get an average for the pair. This number may be odd if one of
the individuals is homozygous over a region shared IBD with the other
individual. In this case the homozygous individual will have 2 IBD
segments (one from each chromosome), while the other individual will
have one (potentially longer) segment.  
**pIBD** is the percent of the diploid chromosome shared IBD between the
pair (percent IBD sharing). It is the total length shared IBD averaged
over all haploid chromosomes.  
**mean_len** the mean length of the IBD segments (in BP)

``` r
x <- gen.simuHaplo_IBD_compare(200,201,BP_len,'Proband_Haplotypes.txt')
## simulation: 21
## pro. 1, chr. 1: 133955228->139683890 
## pro. 2, chr. 1: 133955228->139683890 
## simulation: 23
## pro. 1, chr. 1: 5005746->5233782 
## pro. 2, chr. 2: 5005746->5233782 
## simulation: 28
## pro. 1, chr. 2: 255403506->257844978 
## pro. 2, chr. 2: 255403506->257844978 
## simulation: 31
## pro. 1, chr. 1: 38720302->40230808   
## pro. 2, chr. 1: 38720302->40230808   
## simulation: 38
## pro. 1, chr. 1: 0->13061311  
## pro. 2, chr. 2: 0->13061311  
## simulation: 40
## pro. 1, chr. 1: 35573197->47137119   56259874->56956810  
## pro. 2, chr. 1: 35573197->47137119   56259874->56956810  
## simulation: 53
## pro. 1, chr. 1: 22615994->29274771   
## pro. 2, chr. 1: 22615994->29274771   
## simulation: 56
## pro. 1, chr. 2: 34679627->43403846   
## pro. 2, chr. 2: 34679627->43403846   
## simulation: 61
## pro. 1, chr. 2: 240022538->264278333 
## pro. 2, chr. 1: 240022538->264278333 
## simulation: 68
## pro. 1, chr. 2: 51879916->52633802   
## pro. 2, chr. 2: 51879916->52633802   
## simulation: 73
## pro. 1, chr. 2: 51310548->54609818   
## pro. 2, chr. 2: 51310548->54609818   
## simulation: 74
## pro. 1, chr. 1: 151175282->152924909 
## pro. 2, chr. 1: 151175282->152924909 
## simulation: 91
## pro. 1, chr. 2: 14872374->18464880   
## pro. 2, chr. 2: 14872374->18464880   
## simulation: 93
## pro. 1, chr. 2: 6237617->10809260    18146219->19121423  
## pro. 2, chr. 2: 6237617->10809260    18146219->19121423  
## simulation: 110
## pro. 1, chr. 1: 212210832->217335229 
## pro. 2, chr. 1: 212210832->217335229 
## simulation: 141
## pro. 1, chr. 1: 192766492->195416431 
## pro. 1, chr. 2: 277773793->290000000 
## pro. 2, chr. 1: 192766492->195416431 
## pro. 2, chr. 2: 277773793->290000000 
## simulation: 153
## pro. 1, chr. 2: 43215707->44760433   54782766->59627955  
## pro. 2, chr. 2: 43215707->44760433   54782766->59627955  
## simulation: 161
## pro. 1, chr. 1: 203539222->203849380 
## pro. 2, chr. 2: 203539222->203849380 
## simulation: 163
## pro. 1, chr. 2: 113118980->121223570 
## pro. 2, chr. 2: 113118980->121223570 
## simulation: 171
## pro. 1, chr. 2: 9016517->9794927 278179772->280246309    
## pro. 2, chr. 2: 9016517->9794927 278179772->280246309    
## simulation: 173
## pro. 1, chr. 1: 140920297->160619844 
## pro. 2, chr. 2: 140920297->160619844 
## simulation: 178
## pro. 1, chr. 2: 183266057->185015495 145152645->151192994    
## pro. 2, chr. 1: 183266057->185015495 
## pro. 2, chr. 2: 145152645->151192994 
## simulation: 182
## pro. 1, chr. 2: 255559019->256630856 
## pro. 2, chr. 2: 255559019->256630856 
## simulation: 189
## pro. 1, chr. 1: 185464619->185832373 
## pro. 2, chr. 2: 185464619->185832373 
## simulation: 190
## pro. 1, chr. 2: 235550564->242508387 
## pro. 2, chr. 2: 235550564->242508387 
## simulation: 191
## pro. 1, chr. 2: 287170150->289456345 
## pro. 2, chr. 2: 287170150->289456345 
head(x)
##   simulNo n_seg       pIBD mean_seg_len
## 1      21     2 0.98770028      5728662
## 2      23     2 0.03931655       228036
## 3      28     2 0.42094344      2441472
## 4      31     2 0.26043206      1510506
## 5      38     2 2.25195003     13061311
## 6      40     4 2.11394095      6130429
```
