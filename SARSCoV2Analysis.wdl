version 1.0

workflow SARSCoV2Analysis {
  input {
    File fastqR1
    File? fastqR2
    String samplePrefix
    File bed
    String outputFileNamePrefix = "output"
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

  call variantCalling {
    input:
      sample = samplePrefix,
      bamFile = bowtie2Sensitive.bamFile
  }

  call qcStats {
    input:
      bed = bed,
      sample = samplePrefix,
      bam = bowtie2Sensitive.bamFile
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
    File taxonomicClassification = kraken2.out
    File bam = bowtie2Sensitive.bamFile
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File bl2seqReport = blast2ReferenceSequence.bl2seqReport
    File cvgHist = qcStats.cvgHist
    File genomecvgHist = qcStats.genomecvgHist
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
    String modules = "bowtie2/2.3.5.1 hg38-bowtie-index/2.3.5.1"
    File fastq1
    File fastq2
    String reference = "$HG38_BOWTIE_INDEX_ROOT/hg38_random_index"
    String sample
    Int mem = 12
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    bowtie2 --quiet -x ~{reference} \
    -1 ~{fastq1} -2 ~{fastq2} \
    --un-conc-gz ~{sample}_host_removed_r%.fastq.gz
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_host_removed_r1.fastq.gz"
    File out2 = "~{sample}_host_removed_r2.fastq.gz"
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

  String bamFile = "~{sample}.bam"

  command <<<
    set -euo pipefail

    bowtie2 --sensitive-local -p 4 \
    -x ~{sarsCovidIndex} \
    -1 ~{fastq1} -2 ~{fastq2} \
    | samtools view -b \
    | samtools sort - -o ~{bamFile}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bamFile = "~{bamFile}"
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

  command <<<
    set -euo pipefail

    samtools index ~{bamFile}

    samtools mpileup -aa -d 8000 \
    -uf ~{sarsCovidRef} ~{bamFile} | \
    bcftools call -Mc | tee -a ~{vcfName} | \
    vcfutils.pl vcf2fq -d 10 -D 100000000 | \
    seqtk seq -A - | sed '2~2s/[actg]/N/g' > ~{fastaName}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File vcfFile = "~{vcfName}"
    File consensusFasta = "~{fastaName}"
  }
}

task qcStats {
  input {
    String modules = "bedtools"
    String sample
    File bed
    File bam
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    bedtools coverage -hist -a ~{bed} \
    -b ~{bam} > ~{sample}.cvghist.txt \

    bedtools genomecov -ibam ~{bam} > ~{sample}.genomecvghist.txt \

    bedtools genomecov -d -ibam ~{bam} > ~{sample}.genome.cvgperbase.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File cvgHist = "~{sample}.cvghist.txt"
    File genomecvgHist = "~{sample}.genomecvghist.txt"
  }
}

task blast2ReferenceSequence {
  input {
    String modules = "blast"
    String bl2seq = "/.mounts/labs/gsiprojects/gsi/covid19/sw/bl2seq"
    String reference = "/.mounts/labs/gsiprojects/gsi/covid19/ref/MN908947.3.fasta"
    File consensusFasta
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

    ~{bl2seq} -i ~{consensusFasta} -j ~{reference} -p blastn \
    -W 28 -r 1 -q -2 -F F > bl2seq_report
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

    spades --pe1-1 ~{fastq1} --pe1-2 ~{fastq2} -o ~{sample}.SPAdes

    tar cf - ~{sample}.SPAdes | gzip --no-name > ~{sample}SPAdes.tar.gz
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File sampleSPAdes = "~{sample}SPAdes.tar.gz"
  }
}
