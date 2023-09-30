cwlVersion: v1.0
class: Workflow
inputs:
  genome_fasta:
    type: File
  gtf:
    type: File
  samples:
    type: record
    fields:
      r1:
        type: File
      r2:
        type: File
      sample:
        type: string
outputs:
  gene_expression_files:
    type: File[]
steps:
  trim_reads:
    run: trim_reads.cwl
    in:
      r1: samples.r1
      r2: samples.r2
    out:
      r1_trimmed: trimmed/r1_trimmed.fastq.gz
      r2_trimmed: trimmed/r2_trimmed.fastq.gz
      r1_garbage: trimmed/r1_garbage.fastq.gz
      r2_garbage: trimmed/r2_garbage.fastq.gz

  build_index:
    run: build_star_index.cwl
    in:
      genome_fasta: genome/$(inputs.genome_fasta)
      gtf: genome/$(inputs.gtf)
    out: [star_index]

  map_and_quantify:
    run: mapping_and_quantification.cwl
    in:
      genome_fasta: genome/$(inputs.genome_fasta)
      gtf: genome/$(inputs.gtf)
      samples:
        - r1: trim_reads/r1_trimmed.fastq.gz
          r2: trim_reads/r2_trimmed.fastq.gz
          sample: samples.sample
    out: [gene_expression_files]
