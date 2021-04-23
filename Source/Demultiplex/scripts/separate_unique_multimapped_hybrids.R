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
# 
# p_values_file <- file.path( results_dir, "Hybrids_P_Values/DM_NextSeq/all_NNNCACTAGC_L5Bc_rnc-HTF_1.hyb.max_reads.txt.p_values.txt" ) 
# unique_hybrids_file <- file.path( results_dir, "Hybrids_P_Values/DM_NextSeq/all_NNNCACTAGC_L5Bc_rnc-HTF_1.unique.hyb.p_values.txt" ) 
# multimapped_hybrids_file <- file.path( results_dir, "Hybrids_P_Values/DM_NextSeq/all_NNNCACTAGC_L5Bc_rnc-HTF_1.multimapped.hyb.p_values.txt" )  
#   
p_values_file <- snakemake@input[["p_values_file"]]
unique_hybrids_file <- snakemake@output[["unique_hybrids_file"]]
multimapped_hybrids_file <- snakemake@output[["multimapped_hybrids_file"]]



## -----------------------------------------------------------------------------

p_values_tbl <- vroom::vroom( p_values_file)




## -----------------------------------------------------------------------------

hybrid_repeat_counts  <- p_values_tbl %>%
  group_by( left_hyb_id ) %>%
  summarise( repeat_counts= n()) %>%
  ungroup()

unique_hybrid_ids <- hybrid_repeat_counts %>%
    dplyr::filter( repeat_counts == 1)

multimapped_hybrid_ids <- hybrid_repeat_counts %>%
    dplyr::filter( repeat_counts > 1)



## -----------------------------------------------------------------------------

unique_hybrids_tbl <- p_values_tbl %>%
  inner_join( unique_hybrid_ids, by=c( "left_hyb_id")) %>%
  dplyr::select( -repeat_counts )

vroom::vroom_write( unique_hybrids_tbl, 
                    path=unique_hybrids_file)


multimapped_hybrids_tbl <- p_values_tbl %>%
  inner_join( multimapped_hybrid_ids, by=c( "left_hyb_id")) %>%
  dplyr::select( -repeat_counts )


vroom::vroom_write( multimapped_hybrids_tbl, 
                    path=multimapped_hybrids_file)


