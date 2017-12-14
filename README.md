# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched resequencing data
Cite and read our software note on [BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1985-1)
```
Matthias Beyens, Nele Boeckx, Guy Van Camp, Ken Op de Beeck and Geert Vandeweyer. 
pyAmpli: an amplicon-based variant filter pipeline for targeted enriched resequencing data (2017), 
18:554, DOI:10.1186/s12859-017-1985-1, BMC Bioinformatics
```
 
Questions or ideas? Feel free to contact me: ``matthias dot beyens at uantwerpen dot be``
 
## Content
 1) [Installation](https://github.com/MBeyens/pyAmpli#installation)
 2) [Getting started on example data](https://github.com/MBeyens/pyAmpli#getting-started-on-example-data)
 3) [Getting started on your own data](https://github.com/MBeyens/pyAmpli#getting-started-on-your-own-data)
 4) [What it does?](https://github.com/MBeyens/pyAmpli#what-it-does)
 5) [References](https://github.com/MBeyens/pyAmpli#references)
 
  
## Installation
#### Fully functional on Ubuntu 14.04.4 LTS platform with Python v2.7.6

```
git clone https://github.com/MBeyens/pyAmpli
cd pyAmpli
python setup.py install
 
# if permission issues: python setup.py install --user
# or run with super user

```

If running setup.py without errors, 
the required dependencies should be installed correctly (if not check dependencies.md for manual installation).

```
pyAmpli.py -h
````

If this gives problems finding the binary (rare cases), please append the pyAmpli binary to your systems PATH by doing so

```
export PATH=/ABSOLUTE/PATH/WHERE/YOU/COMPILED/PYAMPLI/pyampli:$PATH
```

 
## Getting started on example data
Create example directory 
```
mkdir example_data
```
Download example data and copy it in example_data. You do not need to provide any BAM index files, pyAmpli will do the indexing for you!
```
https://tinyurl.com/m6u7qr9
```
### Note 1
Be sure you have specified the absolute path of the hg19 genome in the configuration file (config.yaml) before running the example. You can download the full hg19 genome from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/

### Note 2
If not using the full hg19 genome sequence, please download the hg19_chr1.fa.gz from the example data repository and unzip this file. This hg19 human genome contains only chr1 sequences! Be sure to specify the absolute path in the configuration file (config.yaml).
```
gunzip hg19_chr1.fa.gz
```

Run pyAmpli example as follows
```
cd example_data
  
pyAmpli.py somatic \
-c example_config.yaml \
-bn normal_sample_chr1.bam \
-bt tumor_sample_chr1.bam \
-v somatic_variants_chr1.vcf \
-d amplicon_design_chr1.bed \
-od example_data_output/
```

## Getting started on your own data

### Data preparation

Adjust configuration file (bin/config.yaml)
``` 
REQUIRED
    reference:  hg19.fasta          # /ABSOLUTE/PATH/TO/YOUR/REFERENCE.FILE.FASTA
 
OPTIONAL
    cores: 4                        # Running on 4 cores
    change default parameters       # Optimal for somatic filtering
``` 


### Filter modes
``` 
usage: pyAMPLI [-h] {germline,somatic} ...

positional arguments:
  {germline,somatic}  commands
    germline          Input arguments for germline amplicon filter
    somatic           Input arguments for somatic amplicon filter
``` 

#### Somatic modus (default parameters can be found in configuration file)

``` 
usage: pyAMPLI somatic [-h] -bn BAM_NORMAL -bt BAM_TUMOR -v VCF -d DESIGN
                       [-f FILENAME] [-od OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -bn BAM_NORMAL, --bam_normal BAM_NORMAL
                        Binary alignment file, BAM file of the normal sample. Index located in the
                        same directory
  -bt BAM_TUMOR, --bam_tumor BAM_TUMOR
                        Binary alignment file, BAM file of the tumor sample. Index located in the
                        same directory
  -v VCF, --vcf VCF     Variant calling file, VCF file
  -d DESIGN, --design DESIGN
                        Probe/amplicon design file. Contains field, is this
                        order: Contact your manufacturer for details
  -f FILENAME, --filename FILENAME
                        Output file name
  -od OUTDIR, --outdir OUTDIR
                        Output directory
  -c CONFIG, --config CONFIG
                        Location of user-supplied configuration file
```


#### Germline modus (adjust parameters in configuration file)

``` 
usage: pyAMPLI germline [-h] -b BAM -v VCF -d DESIGN [-f FILENAME]
                        [-od OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Binary alignment file, BAM file
  -v VCF, --vcf VCF     Variant calling file, VCF file
  -d DESIGN, --design DESIGN
                        Probe/amplicon design file. Contains field, is this
                        order: Contact your manufacturer for details
  -f FILENAME, --filename FILENAME
                        Output file name
  -od OUTDIR, --outdir OUTDIR
                        Output directory
  -c CONFIG, --config CONFIG
                        Location of user-supplied configuration file
```


## What it does?
See article for more detailed explanation and tool validation.

### Variant categories
Variants are evaluated for each of the following criteria, in the order given here, and assigned to the first matching category. When all criteria are passed, a variant is classified as high quality, corresponding to the label AmpPass.
#### DepthFail | variants with low read evidence
In a first step, variants with insufficient coverage by genuine read-pairs are flagged as low read evidence variants, and not subjected to further variant filtering. FILTER field flags DepthFail and DepthFailTumor/Normal are set, respectively in germline and somatic modus. Users have the flexibility to define their own DepthFail cut-off by adjusting the min_depth_normal and/or min_depth_tumor values in the configuration file.
#### OneAmpPass | variants with panel design limitations
Variants covered by and present in a single theoretical amplicon, as a design limitation, might be more prone for systemic enrichment artefacts. As we have insufficient information to evaluate the reliability of these variants, they are flagged as OneAmpPass, and not subjected to further filtering. 
#### LowAmpFail | variants with low amount of covered amplicons
When variants are covered by multiple theoretical amplicons, we can infer variant reliability based on the number of amplicons containing the variant. Variants covered by more than two overlapping theoretical amplicons are flagged as LowAmpFail if the alternative allele is present in reads corresponding to less than three of these amplicons. Variants covered by just two amplicons are handled separately as LowAmpFail if present in reads corresponding to only one of both amplicons.
#### MatchAmpPass | variants with low amount of covered amplicons
Variants covered by just two theoretical amplicons are handled separately as MatchAmpPass if the alternative allele is present in reads from both amplicons, to indicate the limited discriminative power.
#### PositionFail | positional biases
Variants only present in the first two positions of either 3’ or 5’ read ends are flagged as PositionFail. This enrichment artefact is typically seen in Haloplex gene panels, because fragments are reproducibly generated by restriction enzymes, which cut only recognized sequences and generate non-random fragments [[Reference 1](https://github.com/MBeyens/pyAmpli#references)]. Users can adjust the min_read_pos (default 2) and min_read_pos_fraction (default 10) in the configuration file, i.e. variants will be flagged as PositionFail if more than 10% of the total reads contain the alternative allele in the first two positions of either 3’ or 5’ read ends.
#### NormalFail | low-fraction variants in normal samples
This filter is only applied in somatic mode and is more subjective to user settings. When considering paired tumor-normal samples, somatic variants are not expected to be present in the patient’s paired normal tissue sample. First, this can be indicative for a false-positive somatic variant in the tumor tissue sample, that is in fact a true-positive low-fraction germline variant in the normal sample. Secondly, it might be a systemic enrichment artefact that is more pronounced in the tumor sample and therefore called as somatic. Lastly, it could be a reliable somatic variant. This may be explained by field cancerization, which is the occurrence of genetic, epigenetic and biochemical aberrations in structurally intact cells in histologically normal tissue adjacent to cancerous lesions [[Reference 2](https://github.com/MBeyens/pyAmpli#references)]. By default, somatic variants present in more than 1% of reads from the normal sample are flagged as NormalFail. To allow the effect of field cancerization, the user can adjust the threshold (min_frac) for flagging these variants in the configuration file.
#### AmpPass | threshold-passing variants
As mentioned above, variants passing all user-defined filters are flagged as high quality variants, using the AmpPass label. 


## References

1.	Samorodnitsky E, Jewell BM, Hagopian R, Miya J, Wing MR, Lyon E, Damodaran S, Bhatt D, Reeser JW, Datta J, Roychowdhury S. Evaluation of Hybridization Capture Versus Amplicon‐Based Methods for Whole‐Exome Sequencing. Hum Mutat. 2015;36:903-914.
2.	Slaughter DP, Southwick HW, Smejkal W. Field cancerization in oral stratified squamous epithelium: clinical implications of multicentric origin. Cancer. 1953;6:963–968.
