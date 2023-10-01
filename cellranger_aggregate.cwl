cwlVersion: v1.0
class: CommandLineTool
inputs:
   sample_id:
      type: File
outputs:
   csv:
      type: File
baseCommand: cellranger
arguments:
   - aggr
   - --id
   - $(inputs.sample_id)
   - --csv
   - $(outputs.csv)
stdout: cellrange_aggr.log