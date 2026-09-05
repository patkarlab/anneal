# Generation of Error model for point mutations
A site and mutation specific error model as described in Walkes et. al. Haematologica. 2017;102(9):1549-1557  
[PMID: 28572161](https://doi.org/10.3324/haematol.2017.169136 "doi link") was constructed using Biological Negative Controls (BNCs) for: 
- Single strand consensus sequences using Singleton correction (SSCS using SC)
- Duplex consensus sequence using Singleton correction (DCS using SC)  
Input files and scripts used for model generation can be found in the ``ErrorModel`` and ``ErrorModelDCS`` folder respectively.

## Steps to generate the error model
1. Variant calling for Biological Negative Controls (BNCs).  
Generation of error model requires variant calling data of BNCs. Library preparation was performed using xGen Duplex Seq Adapters followed by a capture protocol. The panel consists of 182 probes. 8 BNCs samples were sequenced and their reads were collapsed to obtain SSCS(SC) and DCS(SC) bam files. Variant calling was performed using [pisces](https://github.com/Illumina/Pisces) variant caller.  
The following steps 2-4 are to be carried out for individual samples.

2. Exclusion of variants with VAF > 0.2  
Sites with VAF > 0.2 were excluded following the methods mentioned in the [article](https://doi.org/10.3324/haematol.2017.169136 "doi link"). This was acheived using the remove_variants_gtr_20.pl script. Command used was 
	```
	perl remove_variants_gtr_20.pl BNC.vcf.gz > BNC.real_removed
	```  
3. Separating multiple variants at the same position  
Positions with more than one ALT variant were split into separate lines using the following command.  
	```
	perl print_multiple_variants_at_same_location.pl BNC.real_removed > BNC_combined.real_removed
	```

4. Adding alt count values for positions without variants  
This step will generate a file of altcount of each base (A/T/G/C) and tagcount (depth for each probe) for all positions covered by probes.  
This step requires 2 input files
	- A text file containing a **non overlapping list of probes**
	- Output of Step 3  

	The following command produces an output file with an extension of .filled  
	```
	perl fill_empty_mips.pl mips_mrd_bal210125_nooverlap.txt BNC_combined.real_removed
	```  

5. Generation of beta matrix  
	Files with *.filled extension for 8 BNC samples were used in this step to generate the beta matrix
	```
	./beta_distribution.py --samples *.filled --output beta_matrix.txt
	```
	This script assumes a fixed error rate of 1/15,000 for sites with no variants detected in one or more BNCs as mentioned in the supplementary methods of the [article](https://doi.org/10.3324/haematol.2017.169136 "doi link"). This can be modified by altering the value of `default_error_rate` variable in this script  

	Two beta matrix files were obtained, one each for SSCS(SC) and DCS(SC). The error rates mentioned for the SmallDeep panel in the supplementary table 5 of this [article](https://doi.org/10.1093/nar/gkz474) were used. The value of `default_error_rate` variable for the SSCS(SC) was taken to be 1/15,000 ( default value was used, the above article mentions a values of 1/13,000) while for the DSC(SC) the value was 1/200,000.