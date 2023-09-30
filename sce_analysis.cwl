cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["Rscript"]
inputs:
   r_script:
       type: File
       inputBinding:
         position: 1
   count_matrix:
       type: File
       inputBinding:
         position: 2
outputs:
    output_file:
       type: File
       outputBinding:
         glob: "Marker_genes_detected.txt"
