version 1.0

workflow sarsCoV2Analysis {
  input {
    File fastq1
    File fastq2
    String samplePrefix
    Boolean trimPrimers
  }

  parameter_meta {
    fastq1: "Read 1 fastq file, gzipped. Can be either targeted or whole transcriptome"
    fastq2: "Read 2 fastq file, gzipped. Optional."
    samplePrefix: "Prefix for output files"
    trimPrimers: "Whether to remove primers used for virus selection"
  }

  meta {
  	author: "Angie Mosquera"
  	email: "amosquera@oicr.on.ca"
  	description: "Classify samples as being either SARS-CoV-2 positive or negative, identify the strain of virus, and produce statistics about the mapping."
  	dependencies: [
    	{
    	    name: "bbmap/38.75",
    	    url: "https://sourceforge.net/projects/bbmap/"
    	},
      {
    	    name: "bowtie2/2.3.5.1",
    	    url: "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zi"
    	},
      {
    	    name: "samtools/1.9",
    	    url: "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2"
    	},
      {
    	    name: "kraken2/2.0.8",
    	    url: "https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads"
    	},
      {
    	    name: "ivar/1.0",
    	    url: "https://github.com/andersen-lab/ivar"
    	},
      {
    	    name: "bcftools/1.9",
    	    url: "https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2"
    	},
      {
    	    name: "vcftools/0.1.16",
    	    url: "https://vcftools.github.io/downloads.html"
    	},
      {
    	    name: "seqtk/1.3",
    	    url: "https://github.com/lh3/seqtk/archive/v1.3.tar.gz"
    	},
      {
    	    name: "bedtools/2.27",
    	    url: "https://github.com/arq5x/bedtools2/releases/tag/v2.27.1"
    	},
      {
    	    name: "blast/2.8.1",
    	    url: "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz"
    	},
      {
    	    name: "spades/3.14.0",
    	    url: "http://cab.spbu.ru/files/release3.14.0/manual.html"
    	}
  	]
  }

  call bbMap {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
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
      bam = select_first([articTrimming.sortedBam, bowtie2Sensitive.bamFile])
  }

  call qcStats {
    input:
      sample = samplePrefix,
      bam = select_first([articTrimming.sortedBam, bowtie2Sensitive.bamFile]),
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
    File? cvgHist = articTrimming.cvgHist
    File vcf = variantCalling.vcfFile
    File consensusFasta = variantCalling.consensusFasta
    File variantOnlyVcf = variantCalling.variantOnlyVcf
    File bl2seqReport = blast2ReferenceSequence.bl2seqReport
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
  String hostMappedBam_ = "~{sample}.host.mapped.bam"
  String hostMappedBai_ = "~{sample}.host.mapped.bam.bai"

  command <<<
    set -euo pipefail

    bowtie2 --quiet -x ~{reference} \
    -1 ~{fastq1} -2 ~{fastq2} \
    --un-conc-gz ~{sample}_host_removed_r%.fastq.gz \
    -S ~{hostMappedSam}

    samtools view -b ~{hostMappedSam} | \
    samtools sort - -o ~{hostMappedBam_}

    samtools index ~{hostMappedBam_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File out1 = "~{sample}_host_removed_r1.fastq.gz"
    File out2 = "~{sample}_host_removed_r2.fastq.gz"
    File hostMappedBam = "~{hostMappedBam_}"
    File hostMappedBai = "~{hostMappedBai_}"
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

  String bamFile_ = "~{sample}.bam"
  String baiFile_ = "~{sample}.bam.bai"

  command <<<
    set -euo pipefail

    bowtie2 --sensitive-local -p 4 -x ~{sarsCovidIndex} \
    -1 ~{fastq1} -2 ~{fastq2} | samtools view -b | samtools sort - -o ~{bamFile_}

    samtools index ~{bamFile_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File bamFile = "~{bamFile_}"
    File baiFile = "~{baiFile_}"
  }
}

task articTrimming {
  input {
    String modules = "ivar/1.0 bedtools"
    File bam
    String sample
    String articBed = "/.mounts/labs/gsiprojects/genomics/SCTSK/analysis/bed/ARTIC-V2.bed"
    Int mem = 8
    Int timeout = 72
  }

  String primertrim = "~{sample}.primertrim"
  String primertrimBam = "~{sample}.primertrim.bam"
  String sortedBam_ = "~{sample}.primertrim.sorted.bam"
  String sortedBai_ = "~{sample}.primertrim.sorted.bam.bai"

  command <<<
    set -euo pipefail

    ivar trim -i ~{bam} -b ~{articBed} -p ~{primertrim}

    samtools sort ~{primertrimBam} -o ~{sortedBam_}

    samtools index ~{sortedBam_}

    bedtools coverage -hist -a ~{articBed} -b ~{bam} > ~{sample}.cvghist.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? sortedBam = "~{sortedBam_}"
    File? sortedBai = "~{sortedBai_}"
    File? cvgHist = "~{sample}.cvghist.txt"
  }
}

task variantCalling {
  input {
    String modules = "bcftools/1.9 samtools/1.9 vcftools/0.1.16 seqtk/1.3 sars-covid-2-bowtie-index/2.3.5.1 sars-covid-2/mn908947.3"
    File bam
    String sample
    String sarsCovidRef = "$SARS_COVID_2_ROOT/MN908947.3.fasta"
    Int mem = 8
    Int timeout = 72
  }

  String vcfName = "~{sample}.vcf"
  String fastaName = "~{sample}.consensus.fasta"
  String variantOnlyVcf_ = "~{sample}.v.vcf"

  command <<<
    set -euo pipefail

    samtools mpileup -aa -uf ~{sarsCovidRef} ~{bam} | \
    bcftools call --ploidy 1 -Mc | tee -a ~{vcfName} | \
    vcfutils.pl vcf2fq -d 10 | \
    seqtk seq -A - | sed '2~2s/[actg]/N/g' > ~{fastaName}

    bcftools mpileup -a "INFO/AD,FORMAT/DP,FORMAT/AD" \
    -d 8000 -f ~{sarsCovidRef} ~{bam} | \
    tee ~{sample}.m.vcf | bcftools call --ploidy 1 -m -v > ~{variantOnlyVcf_}
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File vcfFile = "~{vcfName}"
    File consensusFasta = "~{fastaName}"
    File variantOnlyVcf = "~{variantOnlyVcf_}"
  }
}

task qcStats {
  input {
    String modules = "bedtools samtools/1.9"
    String sample
    File bam
    File hostMappedBam
    Int mem = 8
    Int timeout = 72
  }

  command <<<
    set -euo pipefail

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

    tar cf - ~{sample}.SPAdes | gzip --no-name > ~{sample}.SPAdes.tar.gz
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File sampleSPAdes = "~{sample}.SPAdes.tar.gz"
  }
}