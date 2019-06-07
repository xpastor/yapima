yapima
================
true
2019-06-07

  - [yapima](#yapima)
      - [Description](#description)
      - [Installation](#installation)
          - [Set up the tools using
            *conda*](#set-up-the-tools-using-conda)
          - [Update of
            *IlluminaHumanMethylationEPICmanifest*](#update-of-illuminahumanmethylationepicmanifest)
      - [Execution](#execution)
          - [Definition of variables](#definition-of-variables)
      - [Analysis steps](#analysis-steps)
          - [Preprocessing](#preprocessing)
          - [Copy Number Variation](#copy-number-variation)
          - [Differential Methylation
            Analysis](#differential-methylation-analysis)
      - [Output](#output)
          - [Basic output](#basic-output)
          - [Copy Number Variation](#copy-number-variation-1)
          - [Differential Methylation
            Analysis](#differential-methylation-analysis-1)
          - [Dependencies](#dependencies)
      - [References](#references)

# yapima

## Description

*yapima* is a pipeline developed to simplify and standardize the
preprocessing of the Illumina InfiniumHumanMethylation arrays in order
to easily produce tables with beta-values and/or M-values that can be
later used for further analysis.

For the sake of simplicity in the usage and comprehension of the
pipeline most of the options are hardcoded, and the interaction with it
almost reduced to input/output definition and processes switches.

## Installation

To install yapima just clone this repository with the following command"

    git clone https://github.com/xpastor/yapima.git

The file `environment.yml` contains the packages and software
dependencies to run *yapima*. *yapima* has been tested with the given
versions of the packages and its performance with newer versions is not
guaranteed (but I guarantee that it does not work with older ones).

### Set up the tools using *conda*

You can use *conda* to manage the software requirements and set un an
environment. Instructions on how to install and use *conda* can be found
here: <https://conda.io/projects/conda/en/latest/user-guide/index.html>

Once you have *conda* set up you can create an environment for *yapima*
running the following command:

    conda env create -f environment.yml

Be aware that if you already have a R version installed in your system,
this has priority over the one in your *conda* environment and it will
not have access to the recently created R library by *conda*. In order
to run the R instance from your *conda* environment you need to copy the
provided *env\_vars.sh* files into a specific location of your *conda*
environment. To do so, first you need to activate *yapima* environment:

    conda activate yapima

Then run the following commands:

    cp env_vars_activate.sh ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
    cp env_vars_deactivate.sh ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh

*yapima* can be run either directly using R or through a bash script
(see explanations below, in the *Execution* section). If you want to run
*yapima* directly in R, you need to activate the environment as
explained above. On the other hand, if you want to run the bash wrapper
you need to add the following line to `yapima.sh` right below the
shebang (`#!/bin/bash`):

    . activate yapima

Now you are all set up, and remember to activate the environment
everytime you run *yapima*.

Be aware that conda may not provide the last version of some
bioconductor packages. In such case, after activating the conda
environment start R and install/update the packages manually.

### Update of *IlluminaHumanMethylationEPICmanifest*

Depending on the version of R and array the execution may file due to
incompatibility caused by differences in dimensions between the manifest
provided by bioconductor and the methylation array. In such case install
the manifest provided here.

## Execution

*yapima* requires a configuration file with the specifications of the
files and steps required for the analysis.

The most basic call to run *yapima* is:

    Rscript /path/to/yapima.R /path/to/config_yapima.R

There is also a wrapper bash script that additionally places annotation
files used in the analysis, the sample sheet and a simplified script to
reproduce the results (without plots) in the results folder. In order to
run the analysis using the bash wrapper, execute the following:

    sh /path/to/yapima/yapima.sh -c /path/to/config_yapima.sh

*config\_yapima.R* and *config\_yapima.sh* can be used as templates for
the respective calls. However, be aware that the bash wrapper needs
*config\_yapima.R* as it is here delivered in order to load the values
of the bash variables into the R environment, so if you want to run R
yourself you should modify a copy of the file and not the file itself.

### Definition of variables

#### Requirements

*yapima* uses tools that rely on the standard Illumina structures to
import the data, and although it limits the flexibility of the tool,
allows the user to provide the data without much need of manipulation
since all the required files are usually provided by the core
facilities. There are three requirements regarding input data.

  - **IDAT\_DIR**: the path to the directory with all the idat files.
    All the idat files have to be located in the same directory, without
    any subdirectory structure
  - **SAMPLE\_ANNOTATION**: path to a SampleSheet file, with all the
    relevant annotation for the arrays and the samples
  - **OUTDIR**: directory to store all the files generated by *yapima*

##### A few rules for the SampleSheet

The explanations here are to understand the final structure of the
Sample Sheet, and the only manipulation of the file should be the
addition of sample information columns to the standard Illumina sample
sheet.

The basic SampleSheet is a comma delimited file that always must have
the following columns:

  - Sample\_Name: a unique identifier per sample (avoid purely numerical
    identifiers)
  - Sentrix\_ID: array barcode must correspond to its actual value for
    proper identification of the idat files
  - Sentrix\_Position: array position (6x2 disposition of the 12
    samples) must correspond to the actual position of the sample in the
    array in order to properly identify the idat files; it must follow
    the following pattern: R0\[1-8\]C0\[12\]

It is mandatory that the columns *Sample\_Name*, *Sentrix\_ID* and
*Sentrix\_Position* contain proper values. For the other columns, if no
values are provided they must be empty (i.e. no ‘NA’ values), so the
final structure of the file must be something like this:

``` 
Sample_Name,Sentrix_ID,Sentrix_Position  
Normal1,1234567890,R01C01  
Normal2,1234567890,R02C01  
Normal3,1234567890,R03C01  
Normal4,1234567890,R04C01  
Normal5,1234567890,R05C01  
Normal6,1234567890,R06C01  
Tumor1,1234567890,R01C02  
Tumor2,1234567890,R02C02  
Tumor3,1234567890,R03C02  
Tumor4,1234567890,R04C02  
Tumor5,1234567890,R05C02  
Tumor6,1234567890,R06C02  
```

Optionally, other columns of interest for the differential methylation
analysis or the correction of batch effects can be added. In such case,
if the value of a sample is unknown ‘NA’ must be provided. The values of
categorical variables can not be exclusively numerical (i.e.:
absence/presence of a mutation should be referred as something like
‘ref’/‘mut’, instead of 0/1), although numerical values may be used
if other values are in string format. However, this numerical values
will lose their numerical properties. Also, columns with numerical
variables can include only numerical values or ‘NA’.

``` 
Sample_Name,Sentrix_ID,Sentrix_Position,Group,Hospital  
Normal1,1234567890,R01C01,Normal,AAA  
Normal2,1234567890,R02C01,Normal,AAA  
Normal3,1234567890,R03C01,Normal,BBB  
Normal4,1234567890,R04C01,Normal,NA  
Normal5,1234567890,R05C01,Normal,BBB  
Normal6,1234567890,R06C01,Normal,BBB  
Tumor1,1234567890,R01C02,Tumor,BBB  
Tumor2,1234567890,R02C02,Tumor,BBB  
Tumor3,1234567890,R03C02,Tumor,AAA  
Tumor4,1234567890,R04C02,Tumor,NA  
Tumor5,1234567890,R05C02,Tumor,AAA  
Tumor6,1234567890,R06C02,Tumor,BBB  
```

#### Optional variables

  - **CONDA\_ENV**: name of the conda environment to be activated to
    execute yapima.
  - **BLACKLIST**: stores the path to a file with a set of probes that
    should be discarded from the very beginning of the analysis for
    whatever reason. The probes should be listed per line and identified
    by the Illumina identifier that generally starts with *cg*. Each
    line can have one or two columns, the second one being the sample
    where the value related to the probe should be masked. In summary, a
    probe not related to any sample will be flagged but the values will
    be left. In case a probe is related to a sample, the value for that
    probe in that sample will be replaced by ‘NA’. A probe can be
    related to multiple samples by adding a line for each sample. The
    ‘BLACKLIST’ file can be mixed with both types of probes to be
    masked, global or sample specific.
  - **RUN\_CNV**: valid values are T/TRUE or F/FALSE and toggles the
    execution of CNV analysis
  - **RUN\_DIFFERENTIAL\_METHYLATION**: valid values are T/TRUE or
    F/FALSE and toggles the execution of differential analysis
  - **REMOVE\_SNPS**: possible values are T/TRUE or F/FALSE and flags
    the polymorphic probes of a given population. By default ‘TRUE’
  - **POPULATION**: population from which the allelic frequency is used
    to flag polymorphic probes. By default ‘EUR’. Possible values are
    ’‘, ’EUR’, ‘AMR’, ‘AFR’, ‘ASN’ for 450k arrays and ’‘, ’EUR’,
    ‘AMR’, ‘AFR’, ‘EAS’, ‘SAS’ for EPIC arrays
  - **USE\_PREDICTED\_SEX**: possible values are T/TRUE of F/FALSE, to
    define if the predicted sex has to be used as confounder
  - **BATCH\_VARS**: comma separated list of columns, without spaces,
    from the sample sheet that will be used as confounders or batch
    effects. 'ArrayColumn' and 'ArrayRow', although initially not
    included in the sample sheet, can also be used as batch variables.
  - **SEED**: an integer to allow the repetition of results in random
    processes

## Analysis steps

Since the pipeline is implemented in R and the config file is just a
bash front-end for the pipeline, the values passed to the variables that
control the processess switching need valid R logical values or
expressions that can be evaluated in R context, not in bash. However, to
make it easier, it’s also possible to switch off/on the processes with
0/1.

### Preprocessing

The data is loaded into the R environment using the *minfi* (Aryee et
al. [2014](#ref-minfi)) package from bioconductor. Initially, probes are
flagged using the following values and summed accordingly, being 0 the
value for probes that pass all the filters:

  - 8: probes with low quality in at least half of the samples in a
    cohort
  - 4: probes provided in **BLACKLIST**
  - 2: crossreactive probes for the IlluminaHumanMethylation450k array
    (Chen et al. [2013](#ref-450k)), or for the
    IlluminaHumanMethylationEPIC array (McCartney et al.
    [2016](#ref-epic))
  - 1: if selected, probes with SNPs at the CpG or SBE position with an
    allelic frequency of 0.005 in the selected population for the
    IlluminaHumanMethylation450k array (Chen et al. [2013](#ref-450k))
    or for the IlluminaHumanMethylationEPIC array (McCartney et al.
    [2016](#ref-epic))

Afterwards, *Noob* background correction (Triche et al.
[2013](#ref-noob)) is performed, with an offset of 15 and dye bias
correction without a reference sample. Afterwards, SWAN normalization
(Maksimovic, Gordon, and Oshlack [2012](#ref-SWAN)) is applied, all them
from the *[minfi](https://bioconductor.org/packages/minfi)* (Aryee et
al. [2014](#ref-minfi)) Bioconductor package. Measures with low
detection are replaced by ‘NA’.

### Copy Number Variation

The CNV analysis is done using the
*[conumee](https://bioconductor.org/packages/conumee)* (Hovestadt and
Zapatka, [n.d.](#ref-conumee)) package from bioconductor. Each sample
may take around two minutes of processing, so when executing the
pipeline for exploratory purposes or in large cohorts it may e
convenient to switch off this part.

*conumee* requires a sample set without chromosomal aberrations that are
provided by
*[CopyNeutralIMA](https://bioconductor.org/packages/CopyNeutralIMA)*
(Przybilla and Pastor [2018](#ref-CopyNeutralIMA)).

This analysis may provide a good idea about chromosomal aberrations
although it should not replace an analysis at genomic level, and it
proofs to be as valid to detect sample swaps as CNV analysis using
sequencing data.

### Differential Methylation Analysis

The analysis of differentially methylated positions (DMP) is done using
the *[limma](https://bioconductor.org/packages/limma)* (Ritchie et al.
[2015](#ref-limma)) package from Bioconductor and following the
procedures detailed by the authors. It can perform only paired and
unpaired analysis and reports tables with P-values and adjusted P-values
(Benjamini and Hochberg [1995](#ref-fdr)), and takes as variables of
interest all sample information columns provided by the user in the
SampleSheet except the ones in the **BATCH\_VARS** variable, although
these are also included in the linear model to account for possible
confounding factors.

The analysis on GO and KEGG is done with the package
*[missMethyl](https://bioconductor.org/packages/missMethyl)* (Phipson,
Maksimovic, and Oshlack [2015](#ref-missMethyl2)), which takes into
account the bias due to the number of probes per gene. The probes with
an adjusted p-value below 0.05 are taken as the significant subset.

The detection of differentially methylated regions (DMR) is done using
the *[DMRcate](https://www.bioconductor.org/packages/DMRcate)* (Peters
et al. [2015](#ref-DMRcate)) package from Bioconductor. For two groups
comparisons, the t statistics from the DMP analysis are used, and the
beta log fold change is computed running the standard limma workflow on
the beta values. For comparisons with more than two groups the squared F
statistics from the DMP analysis are provided and the beta log fold
change is set to 0. All the other parameters are left as default.

## Output

### Basic output

When the pipeline is run without switching on any of the steps the
following files are produced:

  - **annotation.txt**: tab delimited file with annotation for each
    probe in InfiniumHumanMethylation450k arrays
  - **rgset.RData**: a *RGSetChannel* object with all the imported
    arrays and annotation
  - **raw\_betas.bed.gz**: compressed BED file with beta values before
    any transformation
  - **bgCorr\_meth.RData**: an object of class *MethylSet* produced
    after transforming data with *ENmix* background correction and *RCP*
    normalization
  - **normalized\_betas.bed.gz**: a compressed BED file with transformed
    beta-values, ready for analysis
  - **normalized\_M.bed.gz**: compressed BED file with transformed
    M-values, suited for statistical analysis
  - **extended\_sample\_sheet.csv**: SampleSheet with an additional
    column for the predicted gender of the samples
  - **genotyping\_betas.txt**: matrix of beta values of 65 genotyping
    probes present in the array

#### Reproducible research

The pipeline also produces two files in order to be able to reproduce
the analysis:

  - **methods\_session\_references.docx**: contains a description of the
    analysis as a paragraph that could be added to the *Methods* section
    of a paper; the session info with all the packages and versions used
    to reproduce the environment used for the analysis and the commit id
    for the version of the pipeline executed; the references to be
    mentioned in case the analysis is included in a paper
  - **script.R**: R script that reproduces a transparent analysis to
    repeat the results produced by the pipeline, with all the conditions
    evaluated by the pipeline already removed

#### QC

By default, the following QC files are also produced:

  - **qc/median\_intensities.pdf**: scatterplots with the median
    intensity of methylation signal against median intensity of
    unmethylation signal. Samples with a combined intensity below 21 are
    shown as bad quality. It is advisable to remove such probes from
    analysis and rerun *yapima*.
  - **qc/bisulfite\_conversion.pdf**: dotplots showing the success of
    the bisulfite conversion. Bad quality samples may show differing
    patterns from the majority of the samples
  - **qc/beta\_distribution.pdf**: density plots of beta values for each
    sample, before and after transformation
  - **qc/beta\_distribution\_heatmap.pdf**: density heatmaps of beta
    values, before and after transformation; these plots help better to
    identify the samples with abnormal beta values distribution
  - **qc/raw\_batch\_PCA\_**\*: PCA plots of samples before
    transformation identifying the samples by groups of
    batch/confounding variables
  - **qc/processed\_batch\_PCA\_**\*: PCA plots of samples after
    transformation identifying the samples by groups of
    batch/confounding variables
  - **qc/samples\_correlation.pdf**: heatmaps showing sample
    correlations using beta values after transformation
  - **qc/genotypes.bed**: BED file with inferred genotypes of the
    samples using the beta values of genotyping probes

Additionally, if there are columns in the SampleSheet other than the
basic ones or identified as batch effects the following plots are
produced:

  - **qc/raw\_PCA\_**\*: PCA plots of samples before transformation
    identifying the samples by variables of interest
  - **qc/processed\_PCA\_**\*: PCA plots of samples after transformation
    identifying the samples by variables of interest

### Copy Number Variation

This is mostly a QC step that produces two files per sample:

  - \*\*qc/CNV\_report/?\_CNV\_report.pdf\*\*: for each sample a pdf
    file with one plot showing all the aberrations in the genome and one
    plot per chromosome, showing the aberrations present in each
    chromosome in more detail
  - \*\*qc/CNV\_report/?\_CNV\_report.txt\*\*: for each sample a tab
    delimited file with information on the regions of the genome with
    chromosomal aberrations

### Differential Methylation Analysis

Each column in the SampleSheet that does not correspond to a basic
column or a batch variables is used to define a differential methylation
analysis, and for each of them the following files are produced:

  - **DMP/DMP\_?.bed.gz**: compressed BED file with some statistical
    parameters (*limma*) like the log fold change (‘logFC’), average
    methylation (‘AveExpr’), the P-value (‘P.value’) and the adjusted
    P-value(‘adj.P.Value’), together with array annotation

  - **DMP/DMP\_PCA\_?.png**: PCA plot using the probes with an adjusted
    P-value smaller than 0.05

  - **DMP/DMP\_cluster\_?.pdf**: hierarchical clustering of the samples
    using the probes with and adjusted P-value \<= 0.05

  - **DMP/GO\_BP\_?.txt**: tab delimited file with the name of the
    category (‘Term’) in the *Biological Process* ontology, the total
    number of probes in the category (‘N’) and the number of significant
    probes in it (‘DE’), the significance of the enrichment (‘P.DE’) and
    the significance after multiple testing correction (‘FDR’)

  - **DMP/GO\_CC\_?.txt**: same as before but for *Cellular Component*

  - **DMP/GO\_MF\_?.txt**: same as before but for *Molecular Function*

  - **DMP/KEGG\_?.txt**: same as before but for KEGG pathways

  - **DMR/DMR\_?.bed**: BED file with the identified DMRs

  - **DMR/DMR\_heatmap\_?.pdf**: heatmap of CpGs from the top DMRs that
    sum up to 500 CpGs, being the top DMRs those with a smaller
    Stouffer.

### Dependencies

The current version of yapima requires R-3.5.1 and Bioconductor 3.8. The
following packages and their dependencies must be installed:

  - biomaRt
  - circlize
  - cluster
  - ComplexHeatmap
  - conumee
  - CopyNeutralIMA
  - DMRcate
  - fpc
  - GenomicRanges
  - ggplot2
  - git2r
  - grid
  - gridExtra
  - knitr
  - limma
  - minfi
  - missMethyl

Additionally, pandoc v2.2.1 or higher is necessary in order to produce
the documentation of the execution.

## References

  - This pipeline has been developed by [Xavier
    Pastor](mailto:xpastor79@gmail.com) at DKFZ, Heidelberg, Germany.

<div id="refs" class="references">

<div id="ref-minfi">

Aryee, Martin J., Andrew E. Jaffe, Hector Corrada-Bravo, Christine
Ladd-Acosta, Andrew P. Feinberg, Kasper D. Hansen, and Rafael A.
Irizarry. 2014. “Minfi: A flexible and comprehensive Bioconductor
package for the analysis of Infinium DNA Methylation microarrays.”
*Bioinformatics* 30 (10): 1363–9.
<https://doi.org/10.1093/bioinformatics/btu049>.

</div>

<div id="ref-fdr">

Benjamini, Y, and Y Hochberg. 1995. “Controlling the False Discovery
Rate: A Practical and Powerful Approach to Multiple Testing.” *Journal
of the Royal Statistical Society Series B* 57: 289–300.

</div>

<div id="ref-450k">

Chen, Yi-an, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria
Grafodatskayak, Brent W. Zanke, Steven Gallinger, Thomas J. Hudson, and
Rosanna Weksberg. 2013. “Discovery of Cross-Reactive Probes and
Polymorphic Cpgs in the Illumina Infinium Human Methylation450
Microarray.” *Epigenetics* 8 (2): 203–9.

</div>

<div id="ref-conumee">

Hovestadt, Volker, and Marc Zapatka. n.d. *Conumee: Enhanced Copy-Number
Variation Analysis Using Illumina Dna Methylation Arrays*. Division of
Molecular Genetics, German Cancer Research Center (DKFZ), Heidelberg,
Germany. <http://bioconductor.org/packages/conumee/>.

</div>

<div id="ref-SWAN">

Maksimovic, Jovana, Lavinia Gordon, and Alicia Oshlack. 2012. “SWAN:
Subset quantile Within-Array Normalization for Illumina Infinium
HumanMethylation450 BeadChips.” *Genome Biology* 13 (6): R44.
<https://doi.org/10.1186/gb-2012-13-6-r44>.

</div>

<div id="ref-epic">

McCartney, Daniel L., Rosie M. Walker, Stewart W. Morris, Andrew M.
McIntosh, David J. Porteous, and Kathryn L. Evans. 2016. “Identification
of Polymorphic and Off-Target Probe Binding Sites on the Illumina
Infinium Methylationepic Beadchip.” *Genomics Data* 9: 22–24.

</div>

<div id="ref-DMRcate">

Peters, Timothy J, Michael J Buckley, Aaron L Statham, Ruth Pidsley,
Katherine Samaras, Reginald V Lord, Susan J Clark, and Peter L Molloy.
2015. “De Novo Identification of Differentially Methylated Regions in
the Human Genome.” *Epigenetics & Chromatin* 8: 6.
<http://www.epigeneticsandchromatin.com/content/8/1/6>.

</div>

<div id="ref-missMethyl2">

Phipson, Belinda, Jovana Maksimovic, and Alicia Oshlack. 2015.
“MissMethyl: An R Package for Analysing Methylation Data from
Illuminas Humanmethylation450 Platform.” *Bioinformatics*, btv560.

</div>

<div id="ref-CopyNeutralIMA">

Przybilla, Moritz, and Xavier Pastor. 2018. *CopyNeutralIMA: Genomic
Copy Neutral samples for Illumina Methylation arrays*.
<https://bioconductor.org/packages/CopyNeutralIMA>.

</div>

<div id="ref-limma">

Ritchie, Matthew E, Belinda Phipson, Di Wu, Yifang Hu, Charity W Law,
Wei Shi, and Gordon K Smyth. 2015. “limma Powers Differential Expression
Analyses for RNA-Sequencing and Microarray Studies.” *Nucleic Acids
Research* 43 (7): e47.

</div>

<div id="ref-noob">

Triche, Timothy J., Daniel J. Weisenberger, David Van Den Berg, Peter W.
Laird, and Kimberly D. Siegmund. 2013. “Low-Level Processing of Illumina
Infinium DNA Methylation BeadArrays.” *Nucleic Acids Research* 41 (7):
e90. <https://doi.org/10.1093/nar/gkt090>.

</div>

</div>
