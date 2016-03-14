# yapima

## Table of Contents
- [Description](#description)
- [Execution](#execution)
- [Input](#input)
    - [Requirements](#requirements)
	    - [A few rules for the SampleSheet](#sample_sheet)
    - [Optional inputs](#optional)
- [Analysis steps](#analysis)
    - [Preprocessing](#preprocessing)
    - [Batch correction](#combat)
	- [Copy Number Variation](#cnv)
	- [Probe Selection](#selection)
	- [Differential Methylation Analysis](#dma)
		- [A few words about the Surrogate Variable Analysis](#sva)
- [Output](#output)
	- [Basic output](#basic)
		- [Reproducible research](#reproduce)
	- [Batch correction](#batch_output)
	- [Copy Number Variation](#cnv_output)
	- [Probe Selection](#selection_output)
	- [Differential Methylation Analysis](#dma_output)
- [Computational requirements](#requirements)
- [References](#references)

<a name="description"></a>
## Description
*yapima* is a pipeline developed to simplify and standardize the preprocessing of the Illumina InfiniumHumanMethylation450k arrays in order to easily produce tables with beta-values and/or M-values that can be later used for further analysis.  
For the sake of simplicity in the usage and comprehension of the pipeline most of the options are hardcoded, and the interaction with it almost reduced to input/output definition and processes switches.

<a name="execution"></a>
## Execution
*yapima* is executed in a very similar way to other pipelines, using a bash configuration file as single parameter:

`sh /path/to/yapima/run_yapima.sh -c /path/to/config_yapima.sh`

<a name="input"></a>
## Input
<a name="requirements"></a>
### Requirements
*yapima* uses tools that rely on the standard Illumina structures to import the data, and although it limits the flexibility of the tool, allows the user to provide the data without much need of manipulation since all the required files are provided by the Core Facility. There are two requirements regarding input data.

* **IDAT\_DIR**: the path to the directory with all the idat files. Initially the core facility provides a compressed file (suffix \_idat.7z) for each array with all the idat files for that array. It is necessary to uncompress these files so that in the end all the idat files are located in the same directory  
* **SAMPLE_ANNOTATION**: path to a SampleSheet file
* **NON\_SPECIFIC\_CG and NON\_SPECIFIC\_CH**: to remove probes mapping at multiple locations the files ['48639-non-specific-probes-Illumina450k\_nonspecific-cg-probes.csv' and '48639-non-specific-probes-Illumina450k\_nonspecific-ch-probes.csv' (*/icgc/ngs_share/general/arrayAnnotations/chen/48639-non-specific-probes-Illumina450k_nonspecific-ch-probes.csv*)](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx) need to be provided (*/icgc/ngs_share/general/arrayAnnotations/chen/48639-non-specific-probes-Illumina450k_nonspecific-cg-probes.csv* and */icgc/ngs_share/general/arrayAnnotations/chen/48639-non-specific-probes-Illumina450k_nonspecific-ch-probes.csv*).

<a name="sample_sheet"></a>
####A few rules for the SampleSheet
The SampleSheet will always be provided by the core facility, one for each array, so in case a project requires multiple arrays the files need to be concatenated and only should contain the samples related to the project. The explanations here are to understand the final structure of the Sample Sheet and besides concatenating the files and adding additional information in general nothing else needs to be done.

The basic SampleSheet is a comma delimited file that always must have the following columns:

* Sample\_Name
* Sample\_Well
* Sample\_Plate
* Sample\_Group
* Pool\_ID
* Sentrix\_ID: array barcode must correspond to its actual value for proper identification of the idat files
* Sentrix\_Position: array position (6x2 disposition of the 12 samples) must correspond to the actual position of the sample in the array in order to properly identify the idat files; it must follow the following pattern: R0[1-8]C0[12]

It's mandatory that the columns 'Sample\_Name', 'Sentrix\_ID' and 'Sentrix\_Position' contain proper values. For the other columns, if no values are provided they must be empty (i.e. no 'NA' values), so the final structure of the file must be something like this:

```
Sample_Name,Sample_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position  
Normal1,,,,,1234567890,R01C01  
Normal2,,,,,1234567890,R02C01  
Normal3,,,,,1234567890,R03C01  
Normal4,,,,,1234567890,R04C01  
Normal5,,,,,1234567890,R05C01  
Normal6,,,,,1234567890,R06C01  
Tumor1,,,,,1234567890,R01C02  
Tumor2,,,,,1234567890,R02C02  
Tumor3,,,,,1234567890,R03C02  
Tumor4,,,,,1234567890,R04C02  
Tumor5,,,,,1234567890,R05C02  
Tumor6,,,,,1234567890,R06C02  
```

Optionally, other columns of interest for the differential methylation analysis or the correction of batch effects can be added. In such case, if the value of a sample is unknown 'NA' must be provided.

```
Sample_Name,Sample_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position,Group,Hospital  
Normal1,,,,,1234567890,R01C01,Normal,AAA  
Normal2,,,,,1234567890,R02C01,Normal,AAA  
Normal3,,,,,1234567890,R03C01,Normal,BBB  
Normal4,,,,,1234567890,R04C01,Normal,NA  
Normal5,,,,,1234567890,R05C01,Normal,BBB  
Normal6,,,,,1234567890,R06C01,Normal,BBB  
Tumor1,,,,,1234567890,R01C02,Tumor,BBB  
Tumor2,,,,,1234567890,R02C02,Tumor,BBB  
Tumor3,,,,,1234567890,R03C02,Tumor,AAA  
Tumor4,,,,,1234567890,R04C02,Tumor,NA  
Tumor5,,,,,1234567890,R05C02,Tumor,AAA  
Tumor6,,,,,1234567890,R06C02,Tumor,BBB  
```

### Optional inputs
* **BLACKLIST**: stores the path to a file with a set of probes that should be discarded from the very beginning of the analysis for whatever reason. The probes should be listed per line and identified by the Illumina identifier that generally starts with 'cg'
* **POLYMORFIC**: stores the path to ['48640-polymorphic-CpGs-Illumina450k\_polymorphic-SBE.csv'](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx), necessary to derive the allelic frequency of the SNPs found at the single base extension (SBE) position (/icgc/ngs_share/general/arrayAnnotations/chen/48640-polymorphic-CpGs-Illumina450k_polymorphic-SBE.csv)

<a name="analysis"></a>
## Analysis steps
Since the pipeline is implemented in R and the config file is just a bash front-end for the pipeline, the values passed to the variables that control the processess switching need valid R logical values or expressions that can be evaluated in R context, not in bash. However, to make it easier for users of our usual pipelines, it's also possible to switch on/off the processes with 0/1.

<a name="preprocessing"></a>
### Preprocessing
The data is loaded into the R environment using the *minfi* [[1]](#minfi_ref) package from bioconductor. Initially, ambiguous probes as described by Chen et al [[2]](#chen_ref) are removed and afterwards *Noob* background correction and *SWAN* normalization [[3]](#swan_ref) are performe if required. Finally, if selected, probes with SNPs at the SBE position that could be present in at least one sample regarding to its allele frequency in european populations as estimated by Chen et al [[2]](#chen_ref) will be removed from the analysis.

<a name="combat"></a>
### Batch correction
This step should be run if some batch effect is detected in the samples. Some typical examples of batch issues are samples hybridized in different array batches or effects related to the position of the samples in the array due to consistent habits during array preparation.

The batch correction is done using the implementation of *ComBat* [[4]](#combat_ref) from the *sva* [[5]](#sva_ref) package and it cannot correct all the effects at once, but needs to be applied one after the other. Therefore, it would be convenient to first correct for the strongest effect and leave the mildest at the end. Because the variables to be corrected are taken in the order they are given to **BATCH\_VARS** it is convenient to put first the variable with the strongest effect and at the end the one with the mildest effect.

It is important to notice that in case some batch variable is strongly related to a variable of interest, the results on the analyses of such variable cannot be trusted and the experimental design should be reviewed.

<a name="selection"></a>
### Probe selection
In some cases it is of interest to try to define some classifier lacking of biological information. To achieve that the pipelines provide the most stable clustering among different probe sets based on bootstraping clustering using the *pvclust* [[6]](#pvclust_ref) package.

This is currently the most challenging (and maybe controversial) part of the pipeline and quite probably needs to be discussed.

This is by far the longest step and may take several hours, processors and a good amount of memory to accomplish it, so unless it is of real interest it may be more convenient to switch it off.

<a name="cnv"></a>
### Copy Number Variation
The CNV analysis is done using the *conumee* [[7]](#conumee_ref) package from bioconductor. This is one of the longest parts of the pipelines. In general, for 20 samples the basic execution of the pipeline (normalization, batch correction and differential methylation analysis) may take around 10 minutes, but the CNV analysis may need other 45 minutes, so when executing the pipeline for exploratory purposes it may be convenient to switch off this part.

This analysis may provide a good idea about chromosomal aberrations although it shouldn't replace an analysis at genomic level, and it proofs to be of value to detect sample swaps when compared to CNV analysis using sequencing data.

<a name="dma"></a>
### Differential Methylation Analysis
There are two different implementations of the differential methylation analysis, and the choice of the approach depends on the value of the variable **SURROGATE\_CORRECTION**.

* When **SURROGATE\_CORRECTION*+ is **switched off** the analysis is done using an implementation of the [*limma*](https://bioconductor.org/packages/release/bioc/html/limma.html) [[8]](#limma_ref) package from bioconductor
* When **SURROGATE\_CORRECTION** is **switched on** the analysis is done using the [*missMethyl*](http://bioconductor.org/packages/release/bioc/html/missMethyl.html) [[9]](#missMethyl_ref) package, that also relies on *limma* but that additionally can estimate surrogate variables and correct the results based on it

Both approaches work in the same way. They are limited to an unpaired analysis and report tables with P-values and adjusted P-values (Benjamini & Hochberg [[10]](#bh_ref)), and they consider as variables of interest all the extra columns of the SampleSheet except the provided in the **BATCH\_VARS** variable.

<a name="sva"></a>
#### A few words about the Surrogate Variable Analysis
Very often, and specially in cancer, there may be sources of variation between samples not related to the specific biological questions and that may mask the results and that cannot be identified, either technical (extractions by pathologists, ...) or biological (cellular heterogeneity, tissue of origen, gender, ...).

There are currently a good number of packages that estimate possible sources of variation and provide ways to correct them. Although at the moment the pipeline does not provide a way to report corrected values, it can include the estimated surrogate variables to correct the differential methylation analysis

<a name="output"></a>
## Output
<a name="basic"></a>
### Basic output
When the pipeline is run without switching on any of the steps the following files are produced:

* **450k\_annotation.txt**: tab delimited file with annotation for each probe in InfiniumHumanMethylation450k arrays
* **filtered\_raw\_betas.gz**: compressed file with a matrix of beta values before any transformation
* **genotyping\_betas.txt**: matrix of beta values of 65 genotyping probes present in the array
* **filtered\_normalized\_meth.RData**: an object of class 'MethylSet' produced after transforming data with *minfi* functions
* **filtered\_normalized\_betas.gz**: compressed file with a matrix of transformed beta-values, ready for analysis
* **filtered\_normalized\_M.gz**: compressed file with a matrix of transformed M-values, suited for statistical analysis
* **extended\_sample\_sheet.csv**: SampleSheet with an additional column for the predicted gender of the samples

<a name="reproduce"></a>
#### Reproducible research
The pipeline also produces some files to explain how the analysis was done and to help to reproduce it:

* **48639-non-specific-probes-Illumina450k\_nonspecific-cg-probes.csv**: comma separated file with probes mapping at multiple locations, as described by Chen et al (Chen, 2013)
* **48639-non-specific-probes-Illumina450k\_nonspecific-ch-probes.csv**: comma separated file with probes mapping at multiple locations, as described by Chen et al (Chen, 2013)
* **48640-polymorphic-CpGs-Illumina450k_polymorphic-SBE.csv**: comma separated file with probes mapping a SNP at the single base extension (SBE), as described by Chen et al (Chen, 2013); only copied if SNPs present in european populations are removed
* **methods.txt**: text file describing the analysis as a paragraph that could be added to the 'Methods' section of a paper; it also contains the commit id for the version of the pipeline that was executed
* **citations.txt**: text file with all the citations that must be mentioned in case the analysis is included in a paper
* **session.pdf**: pdf file giving all the necessary information on R and bioconductor packages version in order to precisely reproduce the environment for the analysis
* **script.R**: R script that reproduces a transparent analysis to repeat the results produced by the pipeline, with all the conditions evaluated by the pipeline already removed

#### QC
By default, the following QC files are also produced:

* **qc/beta\_distribution.pdf**: density plots of beta values for each sample, befor and after transformation
* **qc/beta\_distribution_heatmap.pdf**: density heatmaps of beta values, before and after transformation; these plots help better to identify the samples with abnormal beta values distribution
* **qc/raw\_batch\_PCA\_\***: PCA plots of samples before transformation identifying the samples by groups of batch variables
* **processed\_batch\_PCA\_\***: PCA plots of samples after transformation identifying the samples by groups of batch variables
* **qc/samples\_correlation.pdf**: heatmaps showing sample correlations using beta values after transformation
* **qc/450k\_genotypes.txt**: inferred genotypes of the samples using the beta values of 65 genotyping probes

Additionally, if there are columns in the SampleSheet other than the basic ones or identified as batch effects the following plots are produced:

* **qc/raw\_PCA\_\***: PCA plots of samples before transformation identifying the samples by variables of interest
* **qc/processed\_PCA\_\***: PCA plots of samples after transformation identifying the samples by variables of interest

<a name="batch_output"></a>
### Batch correction
When the batch correction is applied the following files are produced:

* **batch\_corrected.RData**: contains two matrices, one with corrected M-values and another with corrected beta-values
* **batch\_normalized\_M.gz**: compressed file with a matrix with corrected M-values
* **batch\_normalized\_betas.gz** compressed file with a matrix with corrected beta-values, calculated from corrected M-values

<a name="cnv_output"></a>
### Copy Number Variation
This is mostly a QC step that produces two files per sample:

* **qc/CNV\_report/$pid\_CNV\_report.pdf**: pdf file with one plot showing all the aberrations in the genome and one plot per chromosome, showing the aberrations present in each chromosome in more detail
* **qc/CNV\_report/$pid\_CNV\_report.txt**: tab delimited file with information on the regions of the genome with chromosomal aberrations

<a name="selection_output"></a>
### Probe Selection
* **pvclust.RData**: stores an R list with all the resulting clusters
* **betas\_top\_?\_variable\_probes.txt**: matrix with the top variable probes used to produce the cluster with the highest score

#### QC
* **qc/pvclust\_clusters.pdf**: pdf file with all the clusters produced
* **qc/sample\_correlation\_top\_probes.pdf**: pdf with a heatmap showing the correlation between the samples using the beta values of the top variable probes that produced the cluster with the highest score
* **qc/top\_?\_variable\_probes\_PCA\_\*.png**: PCA plots using only the top variable probes that produced the cluster with the highest score, shown by variables of interest

<a name="dma_output"></a>
### Differential Methylation Analysis
Each column in the SampleSheet that does not correspond to a basic column or a batch variables is used to define a differential methylation analysis, and for each of them the following files are produced:

* **differentialMethylation\_PCA\_?.png**: PCA plot using the probes with an adjusted P-value smaller than 0.05
* **?\_diffMeth.gz**: compressed tab delimited file with some statistical parameters (*limma*; *missMethyl*) like the log fold change ('logFC'; ), average methylation ('AveExpr'; ), the P-value ('P.value'; 'p.bayes') and the adjusted P-value('adj.P.Value'; 'p.ebayes.BH'), together with array annotation

<a name="requirements"></a>
## Computational requirements
Depending on the type of analysis the requirements for the analysis may vary a lot. Here is just a short list of examples considering 20 samples each:

* **Raw**: walltime=00:05:00,nodes=1:ppn=1,mem=4gb
* **Background correction**: walltime=00:07:00,nodes=1:ppn=1,mem=4gb
* **Background correction + normalization**: walltime=00:08:00,nodes=1:ppn=1,mem=4gb
* **Full Correction + CNV analysis**: walltime=00:45:00,nodes=1:ppn=1,mem=8gb
* **Full Correction + batch adjustment**: walltime=00:10:00,nodes=1:ppn=1,mem=4gb
* **Full Correction + probe selection**: walltime=03:00:00,nodes=1:ppn=6,mem=30g (user time=01:45:00, cpu time=10:00:00, mem=27gb)
* **Full Correction + probe selection (65 samples)**: walltime=24:00:00,nodes=1:ppn=6,mem=75g (user time=20:15:00, cpu time=120:00:00, mem=65gb)
* **Full Correction + differential methylation analysis**: walltime=00:10:00,nodes=1:ppn=1,mem=5gb (for 1 variable)
* **Full Correction + probe selection + differential methylation**: walltime=03:00:00,nodes=1:ppn=6,mem=30g (user time=01:55:00, cpu time=10:30:00)

The most resources consuming processes are the CNV analysis and the probe selection. In the first case it just means a considerable increase on the total time needed for the analysis but without much penalization on the computational requirements. Regarding the probe selection, the computational requirements scale up depending on the number of processors used for the analysis, since for each processor a copy of the contents in the memory is produced, multiplying the amount of memory needed for the analysis. Additionally, this is a extremely slow step, so taken all together, unless it's of real interest this option souldn't be switched on.

<a name="references"></a>
## References
<a name="minfi_ref"></a>[[1] Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD and Irizarry RA (2014). *Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays.* \_Bioinformatics\_, \*30\*(10), pp. 1363-1369.](http://doi.org/10.1093/bioinformatics/btu049)  
<a name="chen_ref"></a>[2] Chen Y, Lemire M, Choufani S, Butcher DT, Grafodatskayak D, Zanke BW, Gallinger S, Hudson TJ and Weksberg R (2013). *Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium Human Methylation450 microarray.* \_Epigenetics\_, \*8\*(2), pp. 203-209.  
<a name="swan_ref"></a>[[3] Maksimovic J, Gordon L and Oshlack A (2012). *SWAN: Subset quantile Within-Array Normalization for Illumina Infinium HumanMethylation450 BeadChips.* \_Genome Biology\_, \*13\*(6), pp. R44.](http://doi.org/10.1186/gb-2012-13-6-r44)  
<a name="combat_ref"></a>[4] Johnson WE, Rabinovic A and Li C (2007). *Adjusting batch effects in microarray expression data using Empirical Bayes methods.* \_Biostatistics\_, \*8\*(1), pp. 118-127.  
<a name="sva_ref"></a>[5] Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE and Storey JD. *sva: Surrogate Variable Analysis*. R package version 3.14.0.  
<a name="pvclust_ref"></a>[[6] Suzuki R and Shimodaira H (2014). *pvclust: Hierarchical Clustering with P-Values via Multiscale Bootstrap Resampling*. R package version 1.3-2.](http://CRAN.R-project.org/package=pvclust)  
<a name="conumee_ref"></a>[[7] Hovestadt V, Zapatka M (2015). *conumee: Enhanced copy-number variation analysis using Illumina 450k methylation arrays.*, R package version 0.99.4.](http://www.bioconductor.org/packages/release/bioc/html/conumee.html)  
<a name="limma_ref"></a>[8] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies.* \_Nucleic Acids Research\_, \*43\*(7), pp. e47.  
<a name="missMethyl_ref"></a>[9] Phipson B and Oshlack A (2014). *DiffVar: A new method for detecting differential variability with application to methylation in cancer and aging.* \_Genome Biology\_, \*15\*(9), pp. 465.  
<a name="bh_ref"></a>[10] Benjamini Y and Hochberg Y (1995). *Controlling the false discovery rate: a practical and powerful approach to multiple testing.* \_Journal of the Royal Statistical Society Series B\_, \*57\*, pp. 289-300.  
<a name="10"></a>[11] Du P, Zhang X, Huang C, Jafari N, Kibbe WA, Hou L and Lin SM (2010). *Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.* \_BMC Bioinformatics\_, \*11\*, pp. 587.  
___
This pipeline has been developed by Xavier Pastor (x.pastorhostench@dkfz-heidelberg.de) at DKFZ.