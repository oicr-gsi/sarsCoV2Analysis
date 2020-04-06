# sarsCoV2Analysis

Classify samples as being either SARS-CoV-2 positive or negative, identify the strain of virus, and produce statistics about the mapping.

## Overview

![Alt text](docs/SARSCoVAnalysisSummary.png?raw=true)

## Dependencies

* [bbmap 38.75](https://sourceforge.net/projects/bbmap/)
* [bowtie2 2.3.5.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zi)
* [samtools 1.9](https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
* [kraken2 2.0.8](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads)
* [ivar 1.0](https://github.com/andersen-lab/ivar)
* [bcftools 1.9](https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2)
* [vcftools 0.1.16](https://vcftools.github.io/downloads.html)
* [seqtk 1.3](https://github.com/lh3/seqtk/archive/v1.3.tar.gz)
* [bedtools 2.27](https://github.com/arq5x/bedtools2/releases/tag/v2.27.1)
* [blast 2.8.1](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz)
* [spades 3.14.0](http://cab.spbu.ru/files/release3.14.0/manual.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run sarsCoV2Analysis.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastq1`|File|Read 1 fastq file, gzipped. Can be either targeted or whole transcriptome
`fastq2`|File|Read 2 fastq file, gzipped. Can be either targeted or whole transcriptome.
`samplePrefix`|String|Prefix for output files


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`primerBed`|File?|None|
`ampliconBed`|File?|None|


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bbMap.modules`|String|"bbmap/38.75"|
`bbMap.reference`|String|"$BBMAP_ROOT/share/bbmap/resources/adapters.fa"|
`bbMap.trimq`|Int|25|
`bbMap.mem`|Int|8|
`bbMap.timeout`|Int|72|
`bowtie2HumanDepletion.modules`|String|"bowtie2/2.3.5.1 samtools/1.9 hg38-bowtie-index/2.3.5.1"|
`bowtie2HumanDepletion.reference`|String|"$HG38_BOWTIE_INDEX_ROOT/hg38_random_index"|
`bowtie2HumanDepletion.mem`|Int|12|
`bowtie2HumanDepletion.timeout`|Int|72|
`bowtie2HumanDepletion.threads`|Int|8|
`kraken2.modules`|String|"kraken2/2.0.8 kraken2-database/1"|
`kraken2.kraken2DB`|String|"$KRAKEN2_DATABASE_ROOT/"|
`kraken2.mem`|Int|8|
`kraken2.timeout`|Int|72|
`bowtie2Sensitive.modules`|String|"bowtie2/2.3.5.1 sars-covid-2-polymasked-bowtie-index/2.3.5.1 samtools/1.9"|
`bowtie2Sensitive.sarsCovidIndex`|String|"$SARS_COVID_2_POLYMASKED_BOWTIE_INDEX_ROOT/MN908947.3.mask"|
`bowtie2Sensitive.mem`|Int|8|
`bowtie2Sensitive.timeout`|Int|72|
`bowtie2Sensitive.threads`|Int|4|
`articTrimming.modules`|String|"ivar/1.0 bedtools"|
`articTrimming.mem`|Int|8|
`articTrimming.timeout`|Int|72|
`variantCalling.modules`|String|"bcftools/1.9 samtools/1.9 vcftools/0.1.16 seqtk/1.3 sars-covid-2-polymasked/mn908947.3"|
`variantCalling.sarsCovidRef`|String|"$SARS_COVID_2_POLYMASKED_ROOT/MN908947.3.mask.fasta"|
`variantCalling.mem`|Int|8|
`variantCalling.timeout`|Int|72|
`qcStats.modules`|String|"bedtools samtools/1.9"|
`qcStats.mem`|Int|8|
`qcStats.timeout`|Int|72|
`blast2ReferenceSequence.modules`|String|"blast sars-covid-2-polymasked/mn908947.3"|
`blast2ReferenceSequence.reference`|String|"$SARS_COVID_2_POLYMASKED_ROOT/MN908947.3.mask.fasta"|
`blast2ReferenceSequence.mem`|Int|8|
`blast2ReferenceSequence.timeout`|Int|72|
`spadesGenomicAssembly.modules`|String|"spades/3.14.0"|
`spadesGenomicAssembly.mem`|Int|8|
`spadesGenomicAssembly.timeout`|Int|72|


### Outputs

Output | Type | Description
---|---|---
`hostRemovedR1Fastq`|File|None
`hostRemovedR2Fastq`|File|None
`hostMappedBam`|File|None
`hostMappedBai`|File|None
`taxonomicClassification`|File|None
`bam`|File|None
`bai`|File|None
`primertrimSortedBam`|File?|None
`primertrimSortedBai`|File?|None
`cvgHist`|File?|None
`vcf`|File|None
`consensusFasta`|File|None
`variantOnlyVcf`|File|None
`bl2seqReport`|File|None
`genomecvgHist`|File|None
`genomecvgPerBase`|File|None
`hostMappedAlignmentStats`|File|None
`hostDepletedAlignmentStats`|File|None
`spades`|File|None


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
