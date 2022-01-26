#!/bin/bash
cd $1
module load samtools
ls | sed 's/.*\.//' | sort | uniq -c


find . -regex '.*\_bl2seq_report.txt' -exec md5sum {} \;
find . -regex '.*\.consensus.fasta' -exec md5sum {} \;
find . -regex '.*\.genomecvghist.txt' -exec md5sum {} \;
find . -regex '.*\.genome.cvgperbase.txt' -exec md5sum {} \;

find . -regex '.*\_bl2seq_kreport2.txt' -exec md5sum {} \;
find . -regex '.*\_qad_*.fastq.gz' -exec sh -c "echo -n "{}:";zcat {} | sort | md5sum" \;
find . -regex '.*\_host_removed_*.fastq.gz' -exec sh -c "echo -n "{}:";zcat {} | sort | md5sum" \;

find . -regex '.*\.bam' -exec sh -c "echo -n "{}:";samtools flagstat {} | md5sum" \;
find . -regex '.*\.vcf' -exec sh -c "echo -n "{}:";cat {} | grep -v '^#' | cut -f 1-4 | md5sum" \;
