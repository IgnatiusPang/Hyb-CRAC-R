## ----include = FALSE----------------------------------------------------------
## Load libraries 
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(stringr)
p_load(here)


## ----include = FALSE----------------------------------------------------------
## Directories Management
base_dir <-  here::here() 
data_dir <- file.path( base_dir, "Data")
results_dir <- file.path(base_dir, "Results")
source_dir <- file.path(base_dir, "Source")



## -----------------------------------------------------------------------------
# pws_file <- file.path(results_dir,
#                       "PWS", "DM_MiSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.pws")
# 
# hyb_file <- file.path(results_dir,
#                       "Hybrids_Annot",
#                       "DM_MiSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.compressed.pyCRAC_comp_Sa_JKD6009_hybrids_annot.hyb")
# 
# output_file <- file.path(results_dir,
#                       "Hybrids_Annot",
#                       "DM_MiSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.hyb.max_reads.txt")


# pws_file <- file.path(results_dir,
#                       "PWS", "DM_NextSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.pws")
# 
# hyb_file <- file.path(results_dir,
#                       "Hybrids_Annot",
#                       "DM_NextSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.compressed.pyCRAC_comp_Sa_JKD6009_hybrids_annot.hyb")
# 
# output_file <- file.path(results_dir,
#                       "Hybrids_Annot",
#                       "DM_NextSeq",
#                       "all_NNNCACTAGC_L5Bc_rnc-HTF_1.hyb.max_reads.txt")

output_file <- snakemake@output[["MAX_READS"]]
hyb_file <- snakemake@input[["HYB"]]
pws_file <- snakemake@input[["PWS"]]




## -----------------------------------------------------------------------------
pws_tbl <- vroom::vroom( pws_file, 
                         col_names = FALSE) %>%
            dplyr::rename( position = "X2",
                           strand = "X4",
                           reads_depth = "X5") %>%
  dplyr::select( position, strand, reads_depth)


pos_max_pos <- pws_tbl %>%
  dplyr::filter( strand == "+"  ) %>%
  summarise( max_pos = max(position)) %>%
  pull(max_pos)

neg_max_pos <- pws_tbl %>%
  dplyr::filter( strand == "-"  ) %>%
  summarise( max_pos = max(position)) %>%
  pull(max_pos)

pos_max_pos + neg_max_pos == nrow(pws_tbl)

pos_reads_depth <- pws_tbl %>%
  dplyr::filter( strand == "+"  ) %>%
  pull( reads_depth)

pos_positions <- pws_tbl %>%
  dplyr::filter( strand == "+"  ) %>%
  pull( position)


pws_pos_array <- rep(0, pos_max_pos)
pws_pos_array[pos_positions] <- pos_reads_depth

sum(pws_pos_array) ==  sum(pos_reads_depth)

neg_reads_depth <- pws_tbl %>%
  dplyr::filter( strand == "-"  ) %>%
  pull( reads_depth)

neg_positions <- pws_tbl %>%
  dplyr::filter( strand == "-"  ) %>%
  pull( position)


pws_neg_array <- rep(0, neg_max_pos)
pws_neg_array[neg_positions] <- neg_reads_depth

sum(pws_neg_array ) == sum(neg_reads_depth)


## -----------------------------------------------------------------------------
  
hyb_annot_tbl <- vroom::vroom( hyb_file, 
                         col_names = FALSE ) %>%
  dplyr::rename(  hyb_id="X1",
    hyb_seq="X2",
    left_binding_energy="X3",
    left_seq_id="X4",
    left_transcript_start="X5",
    left_transcript_stop="X6",
    left_genomic_start= "X7",
    left_genomic_stop = "X8",
    right_binding_energy="X9",
    right_seq_id="X10",
    right_transcript_start="X11",
    right_transcript_stop="X12",
    right_genomic_start = "X13",
    right_genomic_stop = "X14",
    mapping_score="X15",
    counts_stats="X16") 

hyb_annot_tbl_cln <- hyb_annot_tbl %>%
  dplyr::mutate( hyb_count =  str_match_all(counts_stats, "count_total=(\\d+)")  %>% purrr::map_chr(~as.integer(.[2])) ) %>%
  dplyr::mutate( two_way_merged =  str_match_all(counts_stats, "two_way_merged=(\\d+)")  %>% purrr::map_chr(~as.integer(.[2])) ) %>%
  dplyr::mutate( left_rna_class =   purrr::map_chr(left_seq_id, ~str_split_fixed(., "\\|", 4)[4] ) )  %>%
  dplyr::mutate( right_rna_class =  purrr::map_chr(right_seq_id, ~str_split_fixed(., "\\|", 4)[4] ) ) %>%
  dplyr::mutate( left_chromo =  purrr::map_chr(left_seq_id, ~str_split_fixed(., "\\|", 4)[1] ) ) %>%
  dplyr::mutate( right_chromo =  purrr::map_chr(right_seq_id, ~str_split_fixed(., "\\|", 4)[1] ) ) %>%
  dplyr::mutate( left_name =  purrr::map_chr(left_seq_id, ~str_split_fixed(., "\\|", 4)[2] ) ) %>%
  dplyr::mutate( right_name = purrr::map_chr(right_seq_id, ~str_split_fixed(., "\\|", 4)[2] ) ) %>%
  dplyr::mutate( left_strand = case_when( left_genomic_start < left_genomic_stop ~ "+",
                                          left_genomic_start > left_genomic_stop ~ "-",
                                          TRUE ~ NA_character_    )) %>%
  dplyr::mutate( right_strand = case_when( right_genomic_start < right_genomic_stop ~ "+",
                                           right_genomic_start > right_genomic_stop ~ "-",
                                           TRUE ~ NA_character_    )) 



## -----------------------------------------------------------------------------

get_max_read_depth <- function( start, stop, pws_pos_array, pws_neg_array) {
    if ( start > stop) {
      return( max( pws_neg_array[start:stop])  ) 
    } else if ( start < stop ) {
      return ( max( pws_pos_array[stop:start] ))
    } else {
      stop( "get_max_read_depth: Start position and stop position of hybrid-half is the same position (e.g. hybrid-half of length 1).")
    }

} 




## -----------------------------------------------------------------------------

max_reads_tbl_temp1 <- hyb_annot_tbl_cln %>%
  dplyr::mutate( left_max_reads = purrr::map2_dbl(left_genomic_start, left_genomic_stop, ~get_max_read_depth(.x, .y,  pws_pos_array, pws_neg_array))  ) %>%
  dplyr::mutate( right_max_reads = purrr::map2_dbl(right_genomic_start, right_genomic_stop, ~get_max_read_depth(.x, .y,  pws_pos_array, pws_neg_array))  ) %>%
  dplyr::mutate( smaller_max_depth = ifelse( left_max_reads < right_max_reads, left_max_reads, right_max_reads),
          larger_max_depth = ifelse( left_max_reads > right_max_reads, left_max_reads, right_max_reads) ) %>%
  dplyr::mutate( left_hyb_id = purrr::map_chr( hyb_id, ~as.character(.))) %>%
  dplyr::mutate( left_hyb_count = hyb_count )

colnames( max_reads_tbl_temp1)


max_reads_tbl <- max_reads_tbl_temp1   %>%      
  dplyr::rename( right_hyb_id = "hyb_id" ) %>% 
  dplyr::mutate( right_hyb_count = hyb_count ) %>%
  dplyr::select( left_hyb_id, left_max_reads, left_hyb_count, left_rna_class, left_chromo, 
                 left_name, left_genomic_start, left_genomic_stop, left_strand,
                 right_hyb_id, right_max_reads, right_hyb_count, right_rna_class, right_chromo, 
                 right_name, right_genomic_start, right_genomic_stop, right_strand,
                 smaller_max_depth, larger_max_depth, 
                 two_way_merged,
                 left_binding_energy,
                 left_transcript_start,
                 left_transcript_stop,   
                 right_binding_energy,
                 right_transcript_start,
                 right_transcript_stop,
                 hyb_seq)

# print( paste("output_file =", output_file) )
# 
# write.table(max_reads_tbl, 
#                   output_file, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


if ( is.na( output_file) ) {
   stop( "output_file is not in the input arguments list.")
}

# saveRDS( max_reads_tbl, path=paste0(output_file, ".rds"  ))

vroom::vroom_write( max_reads_tbl, 
                    path=output_file)


