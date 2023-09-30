cwlVersion: v1.0
class: CommandLineTool
inputs:
  genome_fasta:
    type: File
  gtf:
    type: File
outputs:
  star_index:
    type: Directory
baseCommand: STAR
arguments:
  - --runThreadN
  - $(runtime.cores)
  - --runMode
  - genomeGenerate
  - --genomeDir
  - $(runtime.outdir)/star_index/
  - --genomeFastaFiles
  - $(inputs.genome_fasta.path)
  - --sjdbGTFfile
  - $(inputs.gtf.path)
  - --sjdbOverhang
  - 149
stdout: index.log
