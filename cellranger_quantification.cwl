cwlVersion: v1.0
class: CommandLineTool
inputs:
   sample_id:
      type: string
   transcriptome:
      type: File
   fastq_files:
      type: File
baseCommand: cellranger
arguments:
  - count
  - --id
  - $(inputs.sample_id)
  - --transcriptome
  - $(inputs.transcriptome)
  - --fastqs
  - $(inputs.sample_id).fastq
  - --localcores
  - 8
  - --localmem
  - 64
outputs:
  - counts_matrix: $(inputs.sample_id).matrix.mtx    

stdout: cellranger.$(inputs.sample_id).log
stderr: cellranger.$(inputs.sample_id).err

