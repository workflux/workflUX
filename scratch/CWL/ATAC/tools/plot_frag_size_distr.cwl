cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 5000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/tidyverse:1.2.1
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["Rscript", "-e"]
arguments:
  - valueFrom: |
      ${
        var r_cmd = "library(tidyverse); " +
                  "library(gridExtra); " +
                  "size_freqs <- read.table('" + 
                      inputs.fragment_sizes_tsv.path +
                        "', header=F)[,1] %>% as.integer() %>% table(); " +
                  "size_freq_tb <- tibble(frag_size=1:max(as.integer(names(size_freqs)))) %>% " +
                    "full_join( tibble( frag_size=as.integer(names(size_freqs)), freq=size_freqs), by='frag_size') %>% " +
                    "mutate(freq=sapply(freq,function(f) ifelse(is.na(f), 0, f))); " +
                  "plot_ <- ggplot(size_freq_tb) + " +
                    "geom_line(mapping=aes(x=frag_size,y=freq)) + " +
                    "ggtitle('Fragment Size Distribution') + " +
                    "xlab('') + ylab('frequency'); " +
                  "plot_log <- ggplot(size_freq_tb) + " +
                    "geom_line(mapping=aes(x=frag_size,y=freq)) + " +
                    "scale_y_log10() + " +
                    "xlab('fragment size [bp]') + ylab('frequency'); " +
                  "png( '" + inputs.output_basename + "_frag_size_distr.png', " +
                      "width = 850, height = 850); " +
                  "grid.arrange(plot_, plot_log, nrow=2);" +
                  "dev.off();"+
                  "mqc_file_headers <- c("+
                    "'# id: \\\"Fragment Size Distribution\\\"', " +
                    "'# description: \\\"Distribution of fragment sizes (<1000) which typically shows the DNA pitch and multimers of nucleosomes. " + 
                                    "For a better visualization, which also includes higher fragment sizes, please see the files ending with _frag_size_distr.png\\\"', " +
                    "'# pconfig:', " +
                    "'#    xmax: 1000'" +
                    ");" +
                  "mqc_file <- file('" + inputs.output_basename + "_fragm_sizes_mqc.tsv');" +
                  "writeLines(mqc_file_headers, mqc_file);" +
                  "close(mqc_file);" +
                  "write.table(size_freq_tb, file=\'" + inputs.output_basename + "_fragm_sizes_mqc.tsv\', " +
                    "sep='\\t', row.names=F, col.names=F, append=T);" +
                  "print('done')"
        return r_cmd
      }

### INPUT PART:
##################################################
inputs:
  fragment_sizes_tsv:
    type: File
  output_basename:
    type: string
 
### OUTPUT PART:
##################################################
outputs:
  frag_size_distr_plot:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_frag_size_distr.png")
  frag_size_distr_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_fragm_sizes_mqc.tsv")
      
    
    