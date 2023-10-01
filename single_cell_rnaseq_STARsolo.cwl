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
  counts_matrix:
    type: File
  pca_plot:
    type: File
  tsne_plot:
    type: File
  umap_plot:
    type: File
  marker_genes_detected:
    type: File
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
    out: "counts_matrix"

  sce_analysis:
    run: sce_analysis.cwl
    in:
      r_script: single_cell_rnaseq.R
      counts_matrix: map_and_quantify/counts_matrix
    out:
      pca_plot: "PCA.svg"
      tsne_plot: "tSNE.svg"
      umap_plot: "UMAP.svg"
      marker_genes_detected: "Marker_genes_detected.txt"
