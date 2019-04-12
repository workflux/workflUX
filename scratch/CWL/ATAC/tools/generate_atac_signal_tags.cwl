cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
  
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bash", "-c"]
arguments:
  - valueFrom: |
      ${
        var cmd_line = "touch \"" + inputs.output_basename + "_irreg_mappings.bedpe\" " +
            "\"" + inputs.output_basename + "_tn5_bind_region_unsorted.bed\" " +
            "\"" + inputs.output_basename + "_fragments_tn5_incl_tags_unsorted.bed\" " +
            "\"" + inputs.output_basename + "_tn5_center_1bp_unsorted.bed\" " +
            "\"" + inputs.output_basename + "_pot_nucl_bound_tags_unsorted.bed\" " +
            "\"" + inputs.output_basename + "_nucl_free_tags_unsorted.bed\" " +
            "\"" + inputs.output_basename + "_irreg_mappings.bedpe\"; " +
            " awk \'" +
                    "BEGIN {" +
                        "OFS=\"\\t\";" +
                        "chrM_read_count=0;" +
                        "interchrom_map_read_count=0;" +
                        "regular_read_count=0;" +
                        "irregular_read_count=0;" +
                        "too_small_fragment_count=0;" +
                        "nucl_free_fragment_count=0;" +
                        "nucl_bound_fragment_count=0;" +
                        "wrong_strand_orient_count=0;" +
                    "}" +
                    "{" +
                    "	if ( $1==\"chrM\" || $4==\"chrM\") {" +
                            "irregular_read_count += 2;" +
                            "if ( $1==$4 ) {" +
                                "chrM_read_count += 2;" +
                                "print $0, \"chrM\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                            "}" +
                            "else {" +
                                "chrM_read_count += 1;" +
                                "interchrom_map_read_count += 2;" +
                                "print $0, \"interchrom_map\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                            "}" +
                        "}" +
                        "else if ( $1!=$4 ) {" +
                            "irregular_read_count += 2;" +
                            "interchrom_map_read_count += 2;" +
                            "print $0, \"interchrom_map\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                        "}" +
                        "else if ( $9==$10 ) {" +
                            "irregular_read_count += 2;" +
                            "wrong_strand_orient_count += 2;" +
                            "print $0, \"wrong_read_pair_orientation\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                        "}" +
                        "else {" +
                            "right_orientation=0;" +
                            "if ( $9==\"+\" && $2<=$5 && $3<=$6 ){" +
                                "right_orientation=1;" +
                                "fragment_size=$6-$2;" +
                                "start=$2;" +
                                "end=$6;" +
                                "print $1, $2+4, $2+5, $7, $8 > \"" + inputs.output_basename + "_tn5_center_1bp_unsorted.bed\";" +
                                "print $1, $6-5, $6-4, ($7 \"(mate)\"), $8 > \"" + inputs.output_basename + "_tn5_center_1bp_unsorted.bed\";" +
                            "}" +
                            "else if ( $9==\"-\" && $2>=$5 && $3>=$6 ){" +
                                "right_orientation=1;" +
                                "fragment_size=$3-$5;" +
                                "start=$5;" +
                                "end=$3;" +
                                "print $1, $5+4, $5+5, $7, $8 > \"" + inputs.output_basename + "_tn5_center_1bp_unsorted.bed\";" +
                                "print $1, $3-5, $3-4, ($7 \"(mate)\"), $8 > \"" + inputs.output_basename + "_tn5_center_1bp_unsorted.bed\";" +
                            "}" +
                            "else {" +
                                "wrong_strand_orient_count += 2;" +
                                "print $0, \"wrong_read_pair_orientation\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                            "}" +
                            "print fragment_size > \"" + inputs.output_basename + "_fragment_sizes.txt\";" +
                            "if ( fragment_size>=38 && right_orientation ) {" +
                                "regular_read_count += 2;" +
                                "tag_start = start-10;" +
                                "tag_end = start+19;" +
                                "if ( tag_start < 0) { tag_start = 0 }" +
                                "print $1, tag_start, tag_end, $7, $8 > \"" + inputs.output_basename + "_tn5_bind_region_unsorted.bed\";" +
                                "tag_start = end-19;" +
                                "tag_end = end+10;" +
                                "if ( tag_start < 0) { tag_start = 0 }" +
                                "print $1, tag_start, tag_end, ($7 \"(mate)\"), $8 > \"" + inputs.output_basename + "_tn5_bind_region_unsorted.bed\";" +
                                "tag_start = start-10;" +
                                "tag_end = end+10;" +
                                "if ( tag_start < 0) { tag_start = 0 }" +
                                "print $1, tag_start, tag_end, $7, $8 > \"" + inputs.output_basename + "_fragments_tn5_incl_tags_unsorted.bed\";" +
                                "if (fragment_size>=185) { " +
                                    "nucl_bound_fragment_count += 1;" +
                                    "tag_start = start+19;" +
                                    "tag_end = end-19;" +
                                    "if ( tag_end < 0) { tag_end = 0 }" +
                                    "print $1, tag_start, tag_end, $7, $8 > \"" + inputs.output_basename + "_pot_nucl_bound_tags_unsorted.bed\";" +
                                    "tag_start = start-10;" +
                                    "tag_end = start+19;" +
                                    "if ( tag_start < 0) { tag_start = 0 }" +
                                    "print $1, tag_start, tag_end, $7, $8 > \"" + inputs.output_basename + "_nucl_free_tags_unsorted.bed\";" +
                                    "tag_start = end-19;" +
                                    "tag_end = end+10;" +
                                    "if ( tag_start < 0) { tag_start = 0 }" +
                                    "print $1, tag_start, tag_end, ($7 \"(mate)\"), $8 > \"" + inputs.output_basename + "_nucl_free_tags_unsorted.bed\";" +
                                "}" +
                                "else {" +
                                    "nucl_free_fragment_count += 1;" +
                                    "tag_start = start-10;" +
                                    "tag_end = end+10;" +
                                    "if ( tag_start < 0) { tag_start = 0 }" +
                                    "print $1, tag_start, tag_end, $7, $8 > \"" + inputs.output_basename + "_nucl_free_tags_unsorted.bed\";" +
                                "}" +
                            "}" +
                            "else {" +
                                "if ( right_orientation ){" +
                                    "irregular_read_count += 2;" +
                                    "too_small_fragment_count += 1;" +
                                    "print $0, \"too_small_frag_size\" >\"" + inputs.output_basename + "_irreg_mappings.bedpe\";" +
                                "}" +
                            "}" +
                    "	}" +
                    "}" +
                    "END {" +
                        "print \"\# id: \\\"Filtering Statistics\\\"\" > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"\# description: \\\"- This section shows statistics on read filtering: (1) reads pairs that are mapping to different chromosomes as well as " + 
                            "(2) read pairs located on ChrM are filtered out; (3) read pairs that have a wrong orientation towards each other " +
                            "(e.g. both reads on same strand, or reads pointing to different direction) are removed, too. " +
                            "(4) Only the remaining regular reads are used for the fragment size analysis and the generation of atac signal tracks. " +
                            "In addition to the filtering shown here, reads were also selected for high mapping quality and to be mapped in a proper pair." +
                            "\\\"\" > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"\# plot_type: \\\"bargraph\\\"\" > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"chrM_reads\", chrM_read_count > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"interchrom_map_reads\", interchrom_map_read_count > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"wrong_read_pair_orientation\", wrong_strand_orient_count > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"regular_reads\", regular_read_count > \"" + inputs.output_basename + "_filtering_stats_mqc.tsv\";" +
                        "print \"\# id: \\\"Fragment Length Classification\\\"\" > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                        "print \"\# description: \\\"Fragments are classified by their size: (1) nucleosome free fragements are smaller than a typical a nucleosome binding region while "+
                            "(2) potentially nucleosome bound fragments are larger;(3) fragments that are classified as too small are shorter than expected for a Tn5 digestion. " +
                            "For these calculations, the DNA span that is covered by the Tn5 enzyme during the trasposition is taken into account." + 
                            "\\\"\" > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                        "print \"\# plot_type: \\\"bargraph\\\"\" > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                        "print \"too_small_fragments\", too_small_fragment_count > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                        "print \"nucl_free_fragments\", nucl_free_fragment_count > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                        "print \"pot_nucl_bound_fragments\", nucl_bound_fragment_count > \"" + inputs.output_basename + "_frag_size_classification_mqc.tsv\";" +
                    "}" +
                    "\' " + inputs.bedpe_alignm.path;
        cmd_line += ";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + inputs.output_basename + "_tn5_bind_region_unsorted.bed > " + inputs.output_basename + "_tn5_bind_region.bed " +
                    ";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + inputs.output_basename + "_tn5_center_1bp_unsorted.bed > " + inputs.output_basename + "_tn5_center_1bp.bed " +
                    ";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + inputs.output_basename + "_nucl_free_tags_unsorted.bed > " + inputs.output_basename + "_nucl_free_tags.bed " +
                    ";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + inputs.output_basename + "_pot_nucl_bound_tags_unsorted.bed > " + inputs.output_basename + "_pot_nucl_bound_tags.bed " +
                    ";LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n " + inputs.output_basename + "_fragments_tn5_incl_tags_unsorted.bed > " + inputs.output_basename + "_fragments_tn5_incl_tags.bed " ;

        return cmd_line;
      }
        
### INPUT PART:
##################################################
inputs:
  output_basename:
    type: string
  bedpe_alignm:
    type: File
 
### OUTPUT PART:
##################################################
outputs:
  bed_tn5_bind_region_signal:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_bind_region.bed")
  bed_tn5_center_1bp_signal:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_1bp.bed")
  bed_nucl_free_signal:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_nucl_free_tags.bed")
  bed_nucl_bound_signal:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_pot_nucl_bound_tags.bed")
  bed_fragments_tn5_incl_signal:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_fragments_tn5_incl_tags.bed")
  fragment_sizes_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_fragment_sizes.txt")
  filtering_stats_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_filtering_stats_mqc.tsv")
  frag_size_stats_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_frag_size_classification_mqc.tsv")
  irreg_mappings_bedpe:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_irreg_mappings.bedpe")