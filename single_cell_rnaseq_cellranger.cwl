cwlVersion: v1.0
class: Workflow
inputs:
  transcriptome:
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

  cellranger_quantify:
    run: cellranger_quantification.cwl
    in:
      transcriptome: $(inputs.genome_fasta)
      samples:
        - r1: trim_reads/r1_trimmed.fastq.gz
          r2: trim_reads/r2_trimmed.fastq.gz
          sample: samples.sample
    out:
      count_matrix: $(inputs.sample).matrix.mtx

  cellranger_aggregate:
    run: cellranger_aggregate.cwl
    in:
      id: $(inputs.sample)
    out:
      count_matrix: $(inputs.sample).spreadsheet.csv

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
