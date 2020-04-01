version 1.0

workflow SARSCoV2Analysis {
  input {
    File fastqR1
    File? fastqR2
    String samplePrefix
    File bed
    Boolean trimPrimers
  }

  call bbMap {
    input:
      fastq1 = fastqR1,
      fastq2 = fastqR2,
      sample = samplePrefix
  }

  call bowtie2HumanDepletion {
    input:
      fastq1 = bbMap.out1,
      fastq2 = bbMap.out2,
      sample = samplePrefix
  }

  call kraken2 {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix
  }

  call bowtie2Sensitive {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix
  }

  if (trimPrimers) {
    call articTrimming {
      input:
        bam = bowtie2Sensitive.bamFile,
        sample = samplePrefix
    }
  }

  call variantCalling {
    input:
      sample = samplePrefix,
      bamFile = select_first([articTrimming.sortedBam, bowtie2Sensitive.bamFile])
  }

  call qcStats {
    input:
      bed = bed,
      sample = samplePrefix,
      bam = bowtie2Sensitive.bamFile,
      hostMappedBam = bowtie2HumanDepletion.hostMappedBam
  }

  call blast2ReferenceSequence {
    input:
      consensusFasta = variantCalling.consensusFasta
  }

  call spadesGenomicAssembly {
    input:
      fastq1 = bowtie2HumanDepletion.out1,
      fastq2 = bowtie2HumanDepletion.out2,
      sample = samplePrefix
  }

  output {
    File hostRemovedR1Fastq = bowtie2HumanDepletion.out1
    File hostRemovedR2Fastq = bowtie2HumanDepletion.out2
    File hostMappedBam = bowtie2HumanDepletion.hostMappedBam
    File hostMappedBai = bowtie2HumanDepletion.hostMappedBai
    File taxonomicClassification = kraken2.out
    File bam = bowtie2Sensitive.bamFile
    File bai = bowtie2Sensitive.baiFile
    File? primertrimSortedBam = articTrimming.sortedBam
    File? primertrimSortedBai = articTrimming.sortedBai
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File variantOnlyVcf = variantCalling.variantOnlyVcf
    File bl2seqReport = blast2ReferenceSequence.bl2seqReport
    File cvgHist = qcStats.cvgHist
    File genomecvgHist = qcStats.genomecvgHist
    File genomecvgPerBase = qcStats.genomecvgPerBase
    File hostMappedAlignmentStats = qcStats.hostMappedAlignmentStats
    File hostDepletedAlignmentStats = qcStats.hostDepletedAlignmentStats
    File spades = spadesGenomicAssembly.sampleSPAdes
  }
}

task bbMap {
  input {
    String modules = "bbmap/38.75"
    String bbMap = "bbmap"
    File fastq1
    File? fastq2
    String sample
    String reference = "$BBMAP_ROOT/share/bbmap/resources/adapters.fa"
    Int trimq = 25
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    bbmap bbduk in1=~{fastq1} in2=~{fastq2} \
    out1=~{sample}_qad_r1.fastq.gz out2=~{sample}_qad_r2.fastq.gz \
    ref=~{reference} \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=~{trimq}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_qad_r1.fastq.gz"
    File out2 = "~{sample}_qad_r2.fastq.gz"
  }
}

task bowtie2HumanDepletion {
  input {
    String modules = "bowtie2/2.3.5.1 samtools/1.9 hg38-bowtie-index/2.3.5.1"
    File fastq1
    File fastq2
    String reference = "$HG38_BOWTIE_INDEX_ROOT/hg38_random_index"
    String sample
    Int mem = 12
    Int timeout = 72
  }

  String hostMappedSam = "~{sample}.host.mapped.sam"
  String hostMappedBam = "~{sample}.host.mapped.bam"
  String hostMappedBai = "~{sample}.host.mapped.bai"

  command <<<
    set -euo pipefail

    bowtie2 --quiet -x ~{reference} \
    -1 ~{fastq1} -2 ~{fastq2} \
    --un-conc-gz ~{sample}_host_removed_r%.fastq.gz \
    -S ~{hostMappedSam}

    samtools view -b ~{hostMappedSam} | \
    samtools sort - -o ~{hostMappedBam}

    samtools index ~{hostMappedBam} > ~{hostMappedBai}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_host_removed_r1.fastq.gz"
    File out2 = "~{sample}_host_removed_r2.fastq.gz"
    File hostMappedBam = "~{hostMappedBam}"
    File hostMappedBai = "~{hostMappedBai}"
  }
}

task kraken2 {
  input {
    String modules = "kraken2/2.0.8 kraken2-database/1"
    File fastq1
    File fastq2
    String kraken2DB = "$KRAKEN2_DATABASE_ROOT/"
    String sample
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    kraken2 --paired ~{fastq1} ~{fastq2} \
    --db ~{kraken2DB} \
    --report ~{sample}.kreport2
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out = "~{sample}.kreport2"
  }
}

task bowtie2Sensitive {
  input {
    String modules = "bowtie2/2.3.5.1 sars-covid-2-bowtie-index/2.3.5.1 samtools/1.9"
    File fastq1
    File fastq2
    String sample
    String sarsCovidIndex = "$SARS_COVID_2_BOWTIE_INDEX_ROOT/MN908947.3"
    Int mem = 8
    Int timeout = 72
  }

  String samFile = "~{sample}.sam"
  String bamFile = "~{sample}.bam"
  String baiFile = "~{sample}.bai"

  command <<<
    set -euo pipefail

    bowtie2 --sensitive-local -p 4 \
    -x ~{sarsCovidIndex} \
    -1 ~{fastq1} -2 ~{fastq2} \
    -S ~{samFile}

    samtools view -b ~{samFile} | \
    samtools sort - -o ~{bamFile}

    samtools index ~{bamFile}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bamFile = "~{bamFile}"
    File baiFile = "~{baiFile}"
  }
}

task articTrimming {
  input {
    String modules = "ivar/1.0"
    File bam
    String sample
    String articBed = "/.mounts/labs/gsiprojects/genomics/SCTSK/analysis/bed/ARTIC-V2.bed"
    Int mem = 8
    Int timeout = 72
  }

  String bamFile = "~{sample}.bam"
  String primertrim = "~{sample}.primertrim"
  String primertrimBam = "~{sample}.primertrim.bam"
  String sortedBam = "~{sample}.primertrim.sorted.bam"
  String sortedBai = "~{sample}.primertrim.sorted.bai"

  command <<<
    set -euo pipefail

    ivar trim -i ~{bamFile} -b ~{articBed} -p ~{primertrim}

    samtools sort ~{primertrimBam} -o ~{sortedBam}

    samtools index ~{sortedBam}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? sortedBam = "~{sortedBam}"
    File? sortedBai = "~{sortedBai}"
  }
}

task variantCalling {
  input {
    String modules = "bcftools/1.9 samtools/1.9 vcftools/0.1.16 seqtk/1.3 sars-covid-2-bowtie-index/2.3.5.1 sars-covid-2/mn908947.3"
    File bamFile
    String sample
    String sarsCovidRef = "$SARS_COVID_2_ROOT/MN908947.3.fasta"
    Int mem = 8
    Int timeout = 72
  }

  String vcfName = "~{sample}.vcf"
  String fastaName = "~{sample}.consensus.fasta"
  String variantOnlyVcf = "~{sample}.v.vcf"

  command <<<
    set -euo pipefail

    samtools mpileup -aa -uf ~{sarsCovidRef} ~{bamFile} | \
    bcftools call --ploidy 1 -Mc | tee -a ~{vcfName} | \
    vcfutils.pl vcf2fq -d 10 | \
    seqtk seq -A - | sed '2~2s/[actg]/N/g' > ~{fastaName}

    bcftools mpileup -a "INFO/AD,FORMAT/DP,FORMAT/AD" \
    -d 8000 -f ~{sarsCovidRef} ~{sample}.bam | \
    tee ~{sample}.m.vcf | bcftools call --ploidy 1 -m -v > ~{variantOnlyVcf}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File vcfFile = "~{vcfName}"
    File consensusFasta = "~{fastaName}"
    File variantOnlyVcf = "~{variantOnlyVcf}"
  }
}

task qcStats {
  input {
    String modules = "bedtools samtools/1.9"
    String sample
    File bed
    File bam
    File hostMappedBam
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    bedtools coverage -hist -a ~{bed} \
    -b ~{bam} > ~{sample}.cvghist.txt

    bedtools genomecov -ibam ~{bam} > ~{sample}.genomecvghist.txt

    bedtools genomecov -d -ibam ~{bam} > ~{sample}.genome.cvgperbase.txt

    samtools stats ~{hostMappedBam} > ~{sample}.host.mapped.samstats.txt

    samtools stats ~{bam} > ~{sample}.samstats.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File cvgHist = "~{sample}.cvghist.txt"
    File genomecvgHist = "~{sample}.genomecvghist.txt"
    File genomecvgPerBase = "~{sample}.genome.cvgperbase.txt"
    File hostMappedAlignmentStats = "~{sample}.host.mapped.samstats.txt"
    File hostDepletedAlignmentStats = "~{sample}.samstats.txt"
  }
}

task blast2ReferenceSequence {
  input {
    String modules = "blast"
    String reference = "/.mounts/labs/gsiprojects/gsi/covid19/ref/MN908947.3.fasta"
    File consensusFasta
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    blastn -query ~{consensusFasta} -subject ~{reference} \
    -word_size 28 -reward 1 -penalty -2 -dust no > bl2seq_report
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bl2seqReport = "bl2seq_report"
  }
}

task spadesGenomicAssembly {
  input {
    String modules = "spades/3.14.0"
    File fastq1
    File fastq2
    String sample
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    mkdir ~{sample}.SPAdes

    rnaspades.py --pe1-1 ~{fastq1} --pe1-2 ~{fastq2} -o ~{sample}.SPAdes
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File sampleSPAdes = "~{sample}.SPAdes"
  }
}
