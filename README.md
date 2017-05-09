# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched sequencing data
Cite and read our application note on pyAmpli
```
Matthias Beyens, Nele Boeckx, Guy Van Camp, Ken Op de Beeck and Geert Vandeweyer. 
pyAmpli: an amplicon-based variant filter pipeline for targeted enriched sequencing data (2017), doi: 
```

## Installation

```
git clone https://github.com/MBeyens/pyAmpli
cd pyAmpli
python setup.py install 
```

If running setup.py without errors, 
the required dependencies should be installed correctly (dependencies.md).

## Getting started on example data
Create exampe directory 
```
mkdir example_data
```
Download example data 
```
https://tinyurl.com/m6u7qr9
```
Run pyAmpli as follows
```
pyAmpli.py somatic 
-bn example_data/normal_sample_chr1.bam 
-bt example_data/tumor_sample_chr1.bam 
-v example_data/somatic_variants_chr1.vcf 
-d example_data/amplicon_design_chr1.bed 
-od example_data/
```

## Getting started on your own data

### Data preparation
Index your human genome 19 reference
``` 
samtools faidx hg19.fasta
```

Index (.BAI) you alignement (.BAM) file
``` 
samtools index your.BAM
```

Adjust config file (bin/config.yaml)
``` 
REQUIRED
    reference:  hg19.fasta          # /PATH/TO/YOUR/REFERENCE.FILE
 
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

#### Somatic modus (default)

``` 
usage: pyAMPLI somatic [-h] -bn BAM_NORMAL -bt BAM_TUMOR -v VCF -d DESIGN
                       [-f FILENAME] [-od OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -bn BAM_NORMAL, --bam_normal BAM_NORMAL
                        BAM file of the normal sample. Index located in the
                        same directory
  -bt BAM_TUMOR, --bam_tumor BAM_TUMOR
                        BAM file of the tumor sample. Index located in the
                        same directory
  -v VCF, --vcf VCF     VCF file
  -d DESIGN, --design DESIGN
                        Probe/amplicon design file. Contains field, is this
                        order: Contact your manufacturer for details
  -f FILENAME, --filename FILENAME
                        Output file name
  -od OUTDIR, --outdir OUTDIR
                        Output directory
```


#### Germline modus (adjust parameters in config file)

``` 
usage: pyAMPLI germline [-h] -b BAM -v VCF -d DESIGN [-f FILENAME]
                        [-od OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     BAM file. Index located in the same directory
  -v VCF, --vcf VCF     VCF file
  -d DESIGN, --design DESIGN
                        Probe/amplicon design file. Contains field, is this
                        order: Contact your manufacturer for details
  -f FILENAME, --filename FILENAME
                        Output file name
  -od OUTDIR, --outdir OUTDIR
                        Output directory
```


## What it does
   ######(1)	Variants present in a single theoretical amplicon are flagged as OneAmpPass, and not subjected to further variant filtering.
    VARIANT.INFO = OneAmpPass
            
   ######(2)	Variants covered by two theoretical amplicons, both having reads with the alternate allele, are flagged as MatchAmpPass. Variants are flagged as LowAmpFail, if the alternate allele is present in only one of both amplicons.
    VARIANT.INFO = MatchAmpPass
    VARIANT.INFO = LowAmpFail


   ######(3)	Variants with more than two overlapping theoretical amplicons, need the alternate allele to be present in at least three amplicons, otherwise variants are flagged as LowAmpFail.
    VARIANT.INFO = LowAmpFail

   ######(4)	Variants that only occur in the first two positions from either 3’ or 5’ ends of reads are flagged as PositionFail.
    VARIANT.INFO = PositionFail

   ######(5)	Somatic variants present in more than 1% of reads from the normal sample are flagged as NormalFail (only in somatic mode).
    VARIANT.INFO = NormalFail

   ######(6)	Variants passing all filters are flagged as AmpPass.
    VARIANT.INFO = AmpPass
