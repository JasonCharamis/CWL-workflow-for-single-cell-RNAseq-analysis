cwlVersion: v1.0
class: CommandLineTool
inputs:
  r1:
    type: File
  r2:
    type: File
outputs:
  r1_trimmed:
    type: File
  r2_trimmed:
    type: File
  r1_garbage:
    type: File
  r2_garbage:
    type: File
baseCommand: java -jar Trimmomatic-0.39/trimmomatic-0.39.jar
arguments:
  - PE
  - -threads
  - $(runtime.cores)
  - $(inputs.r1.path)
  - $(inputs.r2.path)
  - $(outputs.r1_trimmed.path)
  - $(outputs.r1_garbage.path)
  - $(outputs.r2_trimmed.path)
  - $(outputs.r2_garbage.path)
  - ILLUMINACLIP:/home/Programs/trimmomatic_0.39/adapters/TruSeq3-PE.fa:2:30:10
  - LEADING:3
  - TRAILING:3
  - SLIDINGWINDOW:4:15
  - MINLEN:50
stdout: trimmed.log
