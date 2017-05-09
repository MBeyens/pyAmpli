# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched sequencing data
Read application note for pyAmpli rationale: <I>doi: </I>


## Installation

```
git clone https://github.com/MBeyens/pyAmpli
cd pyAmpli
python setup.py install (--user)
```

Check dependencies at the dependencies folder

Done.

## Usage

### Data preparation
Index (.BAI) you alignement (.BAM) file
``` 
samtools index your.BAM
```

Index your human genome 19 reference
``` 
samtools faidx hg19.fasta
```

### Germline modus

``` 
./pyAmpli.py germline 
      -b ../data/normal_sample_chr1.bam 
      -v ../data/somatic_variants_chr1.vcf 
      -d ../data/amplicon_design_chr1.bed 
      -od ../data/
```


### Somatic modus

``` 
./pyAmpli.py somatic 
      -bn ../data/normal_sample_chr1.bam 
      -bt ../data/tumor_sample_chr1.bam 
      -v ../data/somatic_variants_chr1.vcf 
      -d ../data/amplicon_design_chr1.bed 
      -od ../data/
```


## What it does
    (1)	Variants present in a single theoretical amplicon are flagged as OneAmpPass, 
        and not subjected to further variant filtering. 
    (2)	Variants covered by two theoretical amplicons, both having reads with 
        the alternate allele, are flagged as MatchAmpPass. Variants are flagged 
        as LowAmpFail, if the alternate allele is present in only one of both amplicons.
    (3)	Variants with more than two overlapping theoretical amplicons, 
        need the alternate allele to be present in at least three amplicons, 
        otherwise variants are flagged as LowAmpFail.
    (4)	Variants that only occur in the first two positions from either 3’ or 5’ ends 
        of reads are flagged as PositionFail.
    (5)	Somatic variants present in more than 1% of reads from the normal sample are 
        flagged as NormalFail (only in somatic mode)
    (6)	Variants passing all filters are flagged as AmpPass.
