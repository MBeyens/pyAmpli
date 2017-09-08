# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched resequencing data
Cite and read our application note on pyAmpli
```
Matthias Beyens, Nele Boeckx, Guy Van Camp, Ken Op de Beeck and Geert Vandeweyer. 
pyAmpli: an amplicon-based variant filter pipeline for targeted enriched resequencing data (2017), doi: 
```

## Installation
#### Fully functional on Ubuntu 14.04.4 LTS platform with Python v2.7.6

```
git clone https://github.com/MBeyens/pyAmpli
cd pyAmpli
python setup.py install
# if permission problems: python setup.py install --user

```

If running setup.py without errors, 
the required dependencies should be installed correctly (if not check dependencies.md for manual installation).

```
pyAmpli.py -h
```

## Getting started on example data
Create example directory 
```
mkdir example_data
```
Download example data and copy it in example_data. You do not need to provide any BAM index files, pyAmpli will do the indexing via pysam.
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

Adjust config file (bin/config.yaml)
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


## What it does?
See article