# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched sequencing data
Read application note for pyAmpli rationale: <I>doi: </I>


## Installation

```
git clone https://github.com/MBeyens/pyAmpli
cd pyAmpli
python setup.py install 
```

If running setup.py without errors, 
the required dependencies should be installed correctly (dependencies.md).


## Usage

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

### Somatic modus (default)

``` 
pyAmpli.py somatic 
      -bn data/normal_sample_chr1.bam 
      -bt data/tumor_sample_chr1.bam 
      -v data/somatic_variants_chr1.vcf 
      -d data/amplicon_design_chr1.bed 
      -od data/
```


### Germline modus (adjust parameters in config file)

``` 
pyAmpli.py germline 
      -b data/normal_sample_chr1.bam 
      -v data/somatic_variants_chr1.vcf 
      -d data/amplicon_design_chr1.bed 
      -od data/
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
