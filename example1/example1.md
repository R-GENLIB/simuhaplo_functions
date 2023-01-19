### Using `gen.simuHaplo`

This example will show how to use `gen.simuHaplo` to drop uniquely
labelled founder genotypes down a genealogy to create proband
haplotypes.

First load `GENLIB`

``` r
library(GENLIB)
```

For this example we will use the example genealogy “genea140” that comes
with GENLIB.

``` r
data(genea140)
head(genea140)
```

    ##     ind father mother sex
    ## 1 10086      0      0   1
    ## 2 10087      0      0   2
    ## 3 10102      0      0   1
    ## 4 10103      0      0   2
    ## 5 10128      0      0   1
    ## 6 10129      0      0   2

The genealogy file must use integer ID’s, and contain the ‘ind’,
‘mother’, ‘father’ columns.  
Use the `gen.genealogy` function to create a genealogy object

``` r
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
    ##  Created on : Thu Jan 19 03:52:03 2023

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

**gen** : genealogy object  
**pro** : vector of proband ID’s to simulate. If not specified will
default to all probands in the specified genealogy.  
**ancestors** : Vector of founder ID’s to include. If not specified will
default to all founders of the specified probands. Founders in the
genealogy that are not included in this vector will have their
haplotypes represented by “0”.  
**simulNo** : Number of simulations to run.  
**model** : Integer used to specify the model used to simulate meiosis.
1 = Poisson model. 2 = Zero-truncated Poisson model. 3 = Gamma renewal
process. Both the zero-truncated Poisson, and Gamma models are not
appropriate for sub-chromosomal regions.  
**model_params** : All the models require 1 parameter. **model_params**
should be a vector of length 2, where the first element specifies the
parameter for males, and the second for females. For details on the
model parameters see the paper (cite).  
**cM_len** : A vector of length 2. Length of the chromosome to be
simulated in centimorgans. The first element should be the length for
males, and the second for females.  
**BP_len** : Integer value for the length of the chromosome in base
pairs. Only a single value as it is the same for males and females.  
**physical_map_female** : Option to specify the relationship between
genetic length and physical length. A dataframe with required columns
named “BP” and “cM”, where each row is a pair of points, and the map is
constructed by linear interpolation of the specified points. The first
row should be (0,0), and the last row (BP_len, cM_len). If not specified
the map will be assumed to be linear across the chromosome.  
**physical_map_male** : Same as above but for males.  
**seed** : can seed the random number generator used, for replication
purposes.  
**all_nodes** : Either 1 or 0. If 1 the function will produce an
additional output file for the haplotypes of all the internal
individuals of the genealogy (not just the probands).  
**outDir** : Directory to write output files to. Default is the current
working directory.

For this example we will run 100 simulations on the above genealogy. We
will simulate a hypothetical chromosome 1 with a length of 290,000,000
BP. We use the Poisson process model of meiosis, so the parameters are
the expected number of recombinations. In this case the model parameters
are the same as the Morgan length of the region we are simulating, but
for other models of meiosis this is not the case so we need to specify
the genetic length of the region separately (`cM_len` parameter).

``` r
#print processor specs for the benchmark
specs <- read.csv(text = system("systeminfo /FO csv", intern = TRUE))
print(specs$Processor.s.)
```

    ## [1] "1 Processor(s) Installed.,[01]: Intel64 Family 6 Model 165 Stepping 3 GenuineIntel ~2904 Mhz"

``` r
start_time <- Sys.time() #get speed of simulation

gen.simuHaplo(gen_obj, model = 1, simulNo = 100, model_params = c(1.98, 3.28), 
              cM_len = c(198, 328), BP_len = 290000000, all_nodes = 1)
end_time <- Sys.time()

print(end_time - start_time)
```

    ## Time difference of 45.1532 secs

The function will create the file “Proband_Haplotypes.txt” in the
specified directory. If the **all_nodes** parameter were set to 1 then
it would also create the file “All_nodes_haplotypes.txt”.

### Output files

First we will look at the structure of the “Proband_Haplotypes.txt”
output file:

``` r
f <- file(paste(getwd(),"Proband_Haplotypes.txt",sep="/"),"r")
print(readLines(f,n=3))
```

    ## [1] "100;140"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    ## [2] "{1;217891;0}{0;18272.1;7111508;18271.2;7570900;65329.1;16859174;18271.2;20011437;18338.2;37440836;18324.1;104542440;18321.2;118121360;18324.1;145187125;18436.2;156753221;18426.1;173515573;18427.2;186385556;18426.1;188037576;18442.2;198372301;18442.1;250094195;18440.2;276244650;18441.2;288746686;18443.2;289648440;18323.2;290000000}{0;827275.1;29561408;18540.1;42411197;38248.2;42749531;74789.2;48385140;32373.2;52744823;32377.1;79802352;10151.2;92762851;827275.2;121483188;62183.1;133233024;38645.1;178112488;38649.1;180696781;38650.2;187243498;38649.1;189538403;18513.2;190717230;27146.2;190828291;18513.2;209979245;18516.1;210019816;27146.2;213066858;18452.1;214709802;18311.1;223296470;29668.1;232535696;18311.1;241012090;27146.2;243125227;827275.1;244792319;827274.1;256721502;27144.1;275955956;18513.2;289664766;33432.1;290000000}"                                                                                                                                                                                                                            
    ## [3] "{1;218089;0}{0;120704.2;2894617;120703.2;18931532;42732.2;26089399;42733.2;40095943;59016.2;80885261;39118.2;84429263;39117.1;85523587;42734.1;119126835;32296.1;120416107;32296.2;121127121;43338.2;124960218;32296.2;141074250;44787.2;145972312;43391.2;175597065;32309.2;177973388;32831.1;178744601;32309.2;188531218;32309.1;195111996;215120.2;218584485;206713.2;220735587;215120.2;220915167;10228.1;245739200;436116.1;252929840;39131.1;259499050;436116.1;284660434;32263.2;286185049;32297.2;290000000}{0;13638.2;7806860;13637.1;30686058;13638.2;31155228;24230.2;34033439;16902.2;44394051;751268.1;68765796;212999.1;107604636;32155.1;115622637;18717.1;124380013;56769.2;128312976;27740.2;135288183;32142.1;144559057;32142.2;159853002;34679.2;166704844;33545.1;167237106;32142.2;169946006;16901.2;171359855;16900.2;193821520;24943.1;213936807;32155.2;214045720;32158.1;223773728;32158.2;242753638;32158.1;245667459;32157.2;248392998;63404.2;252003212;63405.2;255075907;30682.2;263581322;18769.1;271052827;58728.1;279112664;40878.1;285629830;43650.2;290000000}"

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

    ## [1] "100;140;34124"                                                                                                                                                                                                                                                  
    ## [2] "{1;33724;2,1,13673949,87542067;5,0,91187884,92503179,125141569,179333883,181242195}{0;10086.2;13673949;10086.1;87542067;10086.2;290000000}{0;10087.1;91187884;10087.2;92503179;10087.1;125141569;10087.2;179333883;10087.1;181242195;10087.2;290000000}"        
    ## [3] "{1;10064;3,1,209793633,255441943,265470175;4,0,104085682,105737689,220393619,226718482}{0;10102.2;209793633;10102.1;255441943;10102.2;265470175;10102.1;290000000}{0;10103.1;104085682;10103.2;105737689;10103.1;220393619;10103.2;226718482;10103.1;290000000}"

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
