
### Using `gen.simuHaplo`

This example will show how to use `gen.simuHaplo` to drop uniquely
labelled founder genotypes down a genealogy to create proband
haplotypes.

First load `GENLIB`, the genealogy file, and create a genealogy object.

``` r
library(GENLIB)
library(readr)
genealogy_file <- read_csv("../example_genealogy.csv", col_types = 'iiii')
```

The genealogy file must use integer ID’s, and contain the ‘ind’,
‘mother’, ‘father’ columns. Use the `gen.genealogy` function to create a
genealogy object

``` r
gen_obj <- gen.genealogy(genealogy_file)
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
    ##  Created on : Wed Dec 07 11:01:12 2022

Now we can simulate transmission of haplotypes from founders to probands
using `gen.simuhaplo`. The function can only simulate one chromosome at
a time, so for multiple chromosomes will require multiple function calls
(as different chromosomes will have different length, genetic length,
and recombination patterns). Note that multiple function calls will
over-write the previous output files, so consider using the `dir.create`
function or something similar in a loop so that the output is written to
different folders. For simulating sub-chromosomal regions, treat the
region as a ‘mini chromosome’ and use appropriate parameters. The
function takes the following parameters:

**`gen`** : The genealogy object created using `gen.genealogy`. Required.  
**`pro`** : A vector of proband ID’s to simulate. Default is all individuals
without children in the specified genealogy.  
**`ancestors`** : A vector of founder ID’s to include. Default 
to all founders of the specified probands. If some of the proband's
founders are not specified in this vector then they will have their
haplotypes represented with ID 0.  
**`simulNo`** : Number of simulations to run. Defaults to 1.  
**`model`** : Integer used to specify the model used to simulate meiosis.
1 = Poisson model. 2 = Zero-truncated Poisson model. 3 = Gamma renewal
process. Both the zero-truncated Poisson, and Gamma models are not
appropriate for sub-chromosomal regions. Defaults to 1.  
**`model_params`** : Vector of length 2, where the first element specifies the
parameter for males, and the second for females. Each model requires a single
parameter but it is different for each model. For details on the model parameters
see the paper (cite). Required parameter, no default.  
**`cM_len`** : A vector of length 2. Length of the chromosome to be
simulated in centimorgans. The first element should be the length for
males, and the second for females. Required parameter, no default.  
**`BP_len`** : An integer value for the length of the chromosome in base
pairs. Required parameter, no default value.    
**`physical_map_female`** : An optional dataframe which specifies the relationship between
genetic length and physical length. The dataframe must have 2 columns
named “BP” and “cM”. Each row of this dataframe is a pair of points, and the map is
constructed by linear interpolation of the specified points. The first
row should be (0,0), and the last row (BP_len, cM_len). If not specified
the map will be assumed to be linear across the chromosome.  
**physical_map_male** : Same as above but for males.  
**seed** : An optional random integer seed that can be used for replication purposes.   
**all_nodes** : Either 1 or 0. If 1 the function will produce an
additional output file for the haplotypes of all the internal
individuals of the genealogy (not just the probands). Default is 0 (no internal haplotypes output).    
**outDir** : A path string to the directory to write output files to. Default is the current
working directory.

For this example we will run 100 simulations on the above genealogy. We
will simulate a hypothetical chromosome 1 with a length of 290,000,000
BP. We use the Poisson process model of meiosis, so the parameters are
the expected number of recombinations.

``` r
#print processor specs for the benchmark
specs <- read.csv(text = system("systeminfo /FO csv", intern = TRUE))
print(specs$Processor.s.)
```

    ## [1] "1 Processor(s) Installed.,[01]: Intel64 Family 6 Model 165 Stepping 3 GenuineIntel ~2904 Mhz"

``` r
start_time <- Sys.time()
gen.simuHaplo(gen_obj, model = 1, simulNo = 100, model_params = c(1.98, 3.28), 
              cM_len = c(198, 328), BP_len = 290000000, all_nodes = 1)
end_time <- Sys.time()
print(end_time - start_time)
```

    ## Time difference of 1.114251 mins

The function will create the file “Proband_Haplotypes.txt” in the
specified directory. If **all_nodes** was set to 1 then it will also
create the file “All_nodes_haplotypes.txt”.

### Output files

First we will look at the structure of the “Proband_Haplotypes.txt”
output file:

``` r
f <- file(paste(getwd(),"Proband_Haplotypes.txt",sep="/"),"r")
print(readLines(f,n=3))
```

    ## [1] "100;227"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    ## [2] "{1;1;0}{0;3025.1;30616597;14301.1;34541036;14301.2;63893366;14298.2;99760071;33673.2;118644715;33673.1;125298348;44766.2;142428525;13592.1;145532117;33673.1;151882042;40595.2;164731238;40575.1;172690859;40575.2;189046528;40575.1;201917179;6023.1;211331071;1037.1;214368507;1040.2;214881050;8012.2;220298080;8011.1;226297983;8011.2;238498257;1040.2;239730179;40595.2;240374735;3025.2;240382486;40595.2;244143132;40575.1;255043550;40595.2;258038121;39067.2;265522069;33009.2;268147599;40595.2;290000000}{0;1025.2;3824201;2710.1;5167829;25397.1;11950741;14305.2;23247716;5978.2;37596722;35659.1;38992543;1084.2;42601986;33673.1;56731610;33673.2;94527849;57624.1;111663110;2928.1;113422872;57624.1;116281766;2928.1;116582891;1084.2;133833596;1082.2;134432004;6228.2;135269211;6228.1;136180957;6229.2;137167420;1082.2;147858715;6229.2;149072779;2928.2;151809477;2929.2;154083167;7413.2;154636496;31415.1;167237008;6075.1;185557240;2528.1;207481426;5790.1;211236267;2528.1;213299054;7321.1;215762193;5790.1;224941897;14413.1;225932513;13845.1;226706675;20172.2;228265312;20210.1;234384167;20210.2;235860700;1038.1;239510975;1038.2;244605147;1022.1;247761557;20210.2;252336128;7415.1;269025112;13843.2;276743254;20210.2;279080146;14412.1;280333452;1026.1;282444517;1027.2;285197015;2555.2;290000000}"
    ## [3] "{1;2;0}{0;6410.2;682604;29024.2;30592273;5785.2;55339862;29027.1;62320895;14331.1;76333132;14330.2;79800317;8407.1;91419917;8044.1;92645434;8612.1;119555872;63929.1;122881428;63929.2;136948228;64293.2;142140700;64293.1;152751488;64284.1;153600371;54735.1;194210518;2871.2;231118111;64293.1;290000000}{0;33587.2;25539494;17305.1;33357600;29913.1;67888407;33429.1;85209360;33428.1;93038230;33429.1;97632249;19858.1;109701027;16750.2;109722112;6171.1;125230871;46227.2;153849828;18852.2;183004681;18852.1;188959471;6954.1;195757578;15018.1;205511134;14643.1;215223415;13268.1;227425403;15019.2;228248406;6171.2;234954443;6171.1;238641904;7776.2;254410574;7776.1;254740193;7196.2;268999215;17625.2;273073292;17625.1;282772551;15018.1;290000000}"

The first line of the output is the number of simulations, and the
number of probands. In the rest of the file there will be one line for
every proband for every simulation. So there will be num_simul \*
num_proband lines in the file. Each line contains information contained
in a set of 3 curly braces. In the first curly bracket is the simulation
number, the proband ID, and a placeholder spot. The following two sets
of brackets contain the haplotype information for the diploid
chromosome. The first of these two is the paternal chromosome copy, and
the second is the maternal.

The haplotypes are encoded **{position ; origin_ID; position; origin_ID;
… origin_ID; position}** where the positions are BP positions along the
chromosome, and origin_ID refers to the ID of the founder from which the
segment between those positions originated. The origin_ID are founder ID
numbers with “.1” or “.2” appended to them, to distinguish between the
diploid copies of the founder chromosome.

The **“All_nodes_haplotypes.txt”** file is formatted in a similar way.
But the first set of brackets contains information on the recombination
positions. The next two brackets have the same format as the other file.

``` r
f <- file(paste(getwd(),"All_nodes_haplotypes.txt",sep="/"),"r")
print(readLines(f,n=3))
```

    ## [1] "100;227;49847"                                                                                                                                                    
    ## [2] "{1;326;0,1,0;4,1,71175636,88550123,117190213,163509740}{0;1017.2;290000000}{0;1016.2;71175636;1016.1;88550123;1016.2;117190213;1016.1;163509740;1016.2;290000000}"
    ## [3] "{1;327;1,0,194115368;2,0,84049807,274408942}{0;1019.1;194115368;1019.2;290000000}{0;1018.1;84049807;1018.2;274408942;1018.1;290000000}"

The first line is the number of simulations, the number of probands, and
the number of individuals. The number of individuals does not include
the founders because they are not output in the
“All_nodes_haplotypes.txt” file (because the founder haplotypes are
trivial, they are all assumed unique).

In the first set of curly braces is the recombination history. This is
used internally to trace back upwards the simulation after the fact. In
these brackets is four fields separated by semicolons. The first field
is the simulation number, the second is the individual ID, the third and
fourth contain the recombination history. These third and fourth fields
are both formatted the same way, to represent the recombinations in the
paternal, and maternal copies of the chromosome respectively.

The recombination history is a sequence of numbers separated by commas.
The **first** number is the “start” of the chromosome. So for example,
if we are looking at the recombination history of the paternal copy of
the chromosome, then if this first value is 0 that means the paternal
chromosome starts with the father’s paternal copy, if it is 1 then it
starts with the fathers maternal copy. When recombinations are
simulated, they are done by creating a new haplotype by stitching
together the two copies of that parent’s chromosome, swapping at every
recombination position. But this operation can start with either of the
parental copy, chosen randomly. The **second** value is the number of
recombinations that occured during that meiosis. The **third** and
following values are the positions of all the recombinations.
