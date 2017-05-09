# pyAmpli: an amplicon-based variant filter pipeline for targeted enriched sequencing data
Read application note for pyAmpli rationale: <I>doi: </I>


## Installation
Download pyAmpli from github ``` git clone https://github.com/MBeyens/pyAMPLI```.

Install it ```python setup.py```

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
./amplicon_filter_v2.py germline 
      -b ../data/normal_sample_chr1.bam 
      -v ../data/somatic_variants_chr1.vcf 
      -d ../data/amplicon_design_chr1.bed 
      -od ../data/
```


### Somatic modus

``` 
./amplicon_filter_v2.py somatic 
      -bn ../data/normal_sample_chr1.bam 
      -bt ../data/tumor_sample_chr1.bam 
      -v ../data/somatic_variants_chr1.vcf 
      -d ../data/amplicon_design_chr1.bed 
      -od ../data/
```
