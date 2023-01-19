### Using gen.simuHaplo_traceback

The traceback function allows the user to specify a proband and a
founder, and given the output of `gen.simuHaplo`, the traceback function
will go through all the simulations, identify the simulations in which
the proband inherits any segment from the specified founder, and then
identify the path of inheritance.

For this example we will use the same “genea140” sample genealogy

``` r
library(GENLIB)
```

    ## Loading required package: Rcpp

``` r
data(genea140)
gen_obj <- gen.genealogy(genea140)
gen_obj
```

    ##  GENLIB: Genealogical Data Software version  1 
    ## 
    ##   Number of individuals :  41523 
    ##   Number of parent-child relations :  68248 
    ##  Number of men  :  20773 
    ##   Number of women :  20750 
    ##   Number of subjects :  140 
    ## 
    ## 
    ##  Genealogical depth : 18
    ## 
    ##  Created on : Thu Jan 19 04:33:27 2023

For this example we will simulate a hypothetical chromosome 1, with a
length of 290,000,000 BP and genetic length of 198cM for males, and
328cM for females. We will use a Poisson process as the model of
meiosis, and we will produce the all_nodes output which is required for
the traceback.

``` r
gen.simuHaplo(gen_obj, model = 1, simulNo = 200, model_params = c(1.98, 3.28), 
              cM_len = c(198, 328), BP_len = 290000000, all_nodes = 1)
```

Now that we have the simulation results we can use the traceback
function to investigate the relationship between any of the simulated
probands and a select founder. For this illustration we arbitrarily pick
a proband (can see a list of probands with `gen.pro()`), and one of
their founders (can get founders with `gen.findFounders()`). For this
example we will use proband 408992 and founder 863241. First let us
visualize the relationship:

``` r
sub_genealogy <- gen.branching(gen_obj, 408992)
gen.graph(sub_genealogy, cex=0.5)
```

!<img src="gengraph.png" width=750><!-- -->

Now we can use the **gen.simuHaplo_traceback** function to identify
every simulation where proband 408992 inherits a segment from founder
863241 and identify the path of inheritance. The function takes the
following parameters:  
**gen: ** The genealogy object that was used for the simulation **proID:
** The ID of the proband **ancestorID: ** ID of the founder. Any founder
from this ancestor appearing in the proband haplotype will be traced
back. **all_nodes_path: ** The path to the ‘All_nodes_haplotypes.txt’
file created from the simulation. **proband_haplotypes_path: ** The path
to the ‘Proband_Haplotypes.txt’ file created from the simulation.

``` r
traceback <- gen.simuHaplo_traceback(gen_obj, proID =408992, ancestorID = 863241, 
                                     all_nodes_path = 'All_nodes_haplotypes.txt',
                                     proband_haplotypes_path = 'Proband_Haplotypes.txt')
## input file paths:
## All_nodes_haplotypes.txt
## Proband_Haplotypes.txt
## 
## path: 1 863840 863843 827218 827223 863220 863222 863226 863229 863232 863236 863238 863241 
## path: 2 863839 
## path: 3 863839 863842 827218 827223 863220 863222 863226 863229 863232 863236 863238 863241 
## path: 4 863840 863844 827220 827228 863220 863222 863226 863229 863232 863236 863238 863241 
## path: 5 863839 863841 827220 827228 863220 863222 863226 863229 863232 863236 863238 863241 
## path: 6 863840 863843 827219 827224 863220 863222 863226 863229 863232 863236 863238 863241 
## path: 7 863839 863842 827219 827224 863220 863222 863226 863229 863232 863236 863238 863241
```

The traceback function will print all the paths that it encounters in
the simulation results, in the order that it encounters them. Note that
if a path does not end with the ancestor ID that means the path
bifurcates at the last internal ancestor that was printed. If an
internal ancestor is homozygous for the founder segment it is possible
for it to transmit a longer segment that has been stitched together
through a recombination. In this case the traceback is terminated at
this internal ancestor.

The function returns a dataframe that specifies the path, and the length
of the segment.

``` r
head(traceback)
```

    ##   simulNo seg_length path_n
    ## 1      11    3237764      1
    ## 2      60   42080693      2
    ## 3      60     660504      3
    ## 4      60   28581328      4
    ## 5      84     837747      5
    ## 6      87   13843592      5
