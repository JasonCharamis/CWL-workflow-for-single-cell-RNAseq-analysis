cwlVersion: v1.0
class: CommandLineTool
inputs:
  genome_dir:
    type: Directory
  sample:
    type: string
  r1_trimmed:
    type: File
  r2_trimmed:
    type: File
outputs:
  gene_expression:
    type: File
baseCommand: STAR
arguments:
  - --runThreadN
  - 20
  - --genomeDir
  - $(inputs.genome_dir)
  - --readFilesIn
  - $(inputs.sample)_1.trimmed.fastq
  - $(inputs.sample)_2.trimmed.fastq
  - --outSAMtype
  - BAM SortedByCoordinate
  - --outFileNamePrefix
  - $(inputs.sample.path).bam
  - --soloMultiMappers
  - EM
  - --soloType
  - Droplet
  - --soloCBwhitelist
  - whitelist
  - --clipAdapterType
  - CellRanger4
  - --outFilterScoreMin
  - 30
  - --soloCBmatchWLtype
  - 1MM_multi_Nbase_pseudocounts
  - --soloUMIfiltering
  - MultiGeneUMI_CR
  - --soloUMIdedup
  - 1MM_CR

stdout: index.log