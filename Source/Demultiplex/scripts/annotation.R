## ----include = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
## Load libraries 
if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(stringr)
p_load(here)
p_load(igraph)
p_load(magrittr)

p_load(Biostrings)
p_load(rtracklayer )
p_load(furrr)



## ----include = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
## Directories Management
base_dir <-  here::here() 
data_dir <- file.path( base_dir, "Data")
results_dir <- file.path(base_dir, "Results")
source_dir <- file.path(base_dir, "Source")
hyb_p_val_dir <- file.path( results_dir, "Hybrids_P_Values")

dataset_type <-    "JT_SM"   #  "SM_Only" # "JT_Only" #

merged_hyb_p_val_dir <- file.path( results_dir, "Merged_Hybrids_P_Values", dataset_type)

source( file.path( source_dir, "Common", "common_functions.R") ) 

create_dir_if_not_exists(merged_hyb_p_val_dir )

num_groups_for_calc <- 30  ## split gff file into groups to aid quicker computation
 



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# hyb_file <- "/home/ignatius/PostDoc/2020/StaphCLASH2020/Results/Hybrids_Merged/DM_MiSeq/all_NNNCACTAGC_L5Bc_rnc-HTF_1.compressed.pyCRAC_comp_Sa_JKD6009_hybrids.merged.hyb"
# 
# gtf_file <- "/home/ignatius/PostDoc/2020/StaphCLASH2020/Data/Genomic/JKD6009/Annotation/Hyb_pyCRAC/compiled_edited_features.gtf"
# 
# rna_type_file <- "/home/ignatius/PostDoc/2020/StaphCLASH2020/Data/Genomic/JKD6009/Annotation/Hyb_pyCRAC/rna_type_priority.tab"

hyb_file <- snakemake@input[["HYB_FILE"]]
gtf_file <- snakemake@input[["GTF_FILE"]]
rna_type_file <- snakemake@input[["RNA_TYPE_FILE"]]
chormosome_feature_file <- snakemake@input[["CHROMOSOME_FEATURE_FILE"]]
output_file <- snakemake@output[["hyb_annot_output"]]

# hyb_file <- "/media/ignatius/Dinosaur/Winton_CLASH/DM_NextSeq/Processed_Data/Hyb_Merge/all_NNNTCTCTAGC_L5Bd_rnc-HTF_2.compressed.pyCRAC_comp_Sa_JKD6009_hybrids.merged.hyb"
# gtf_file <- "/home/ignatius/PostDoc/2020/StaphCLASH2020/Data/Genomic/JKD6009/Annotation/Hyb_pyCRACv/1_JKD6009_genomic_merge_features_simple_annot_checked_edited.gtf"
# rna_type_file <- "/home/ignatius/PostDoc/2020/StaphCLASH2020/Data/Genomic/JKD6009/Annotation/Hyb_pyCRACv/rna_type_priority.tab"
# output_file <- "/media/ignatius/Dinosaur/Winton_CLASH/DM_NextSeq/Processed_Data/Hyb_Annot/all_NNNTCTCTAGC_L5Bd_rnc-HTF_2.compressed.pyCRAC_comp_Sa_JKD6009_hybrids_annot.hyb"



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------


hyb_tbl <- vroom::vroom( hyb_file, col_names = FALSE) %>%
  mutate(row_id = row_number())

hyb_colnames <- c(
"left_hyb_id",
"hyb_seq",
"binding_energy",
"left_chromo",	"left_read_start", "left_read_stop",	"left_genomic_start",	"left_genomic_stop",	"left_strand"	,


"right_chromo",	"right_read_start", "right_read_stop",	"right_genomic_start",	"right_genomic_stop",	"right_strand",
"optional_text", "row_id"
)

# Column 1, unique sequence identifier.
# Column 2, read sequence (truncated here for clarity).
# Column 3, predicted binding energy in kcal/mol.
# Columns 4–9, mapping information for first fragment of read: name of matched transcript, coordinates in read, coordinates in transcript, mapping score.
# Columns 10–15, mapping information for second fragment of read.
# Column 16 (optional, not shown here), list of annotations in the format: ‘‘feature1=value1; feature2=value2;..."

colnames( hyb_tbl ) <- hyb_colnames

# hyb_tbl



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
gff_obj <- import.gff2(gtf_file)

gff_tbl <- as.data.frame( gff_obj )

# "_" in hyb file converted to "-"

chromosome_features_tbl <- vroom::vroom( chormosome_feature_file, colnames=c("gtf_chromosome_name", "hyb_chromosome_name") )

rna_coordinates <- gff_tbl %>%
  left_join(chormosome_feature_file
            , by=c("seqnames"="gtf_chromosome_name") )
  mutate( seqnames = as.character(hyb_chromosome_name),
          type = as.character( type)) %>%
  mutate( feature_name = pmap_chr( list( seqnames, gene_id,  transcript_id, type ), ~paste0(c(..1, ..2, ..3, ..4), collapse="|")      )) %>%
  dplyr::select ( seqnames, start, end, width, strand, feature_name, type, rna_type) %>%
  dplyr::rename( feature_chromosome = "seqnames",
                 feature_start = "start", 
                 feature_end = "end",
                 feature_strand = "strand",
                 feature_type  = "type" ) %>%
  dplyr::filter( feature_type != "gene")

## substitute into column 4 and column 10


# 
# gff_tbl %>%
#   distinct( type)

group_size <- ceiling( max( rna_coordinates$feature_end ) / num_groups_for_calc ) 

rna_coordinates <- rna_coordinates %>%
  mutate( group_id = (feature_end %/% group_size)   + 1 ) 



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
hyb_tbl <- hyb_tbl %>%
  mutate( left_group_id_stop = (left_genomic_stop %/% group_size)   + 1 )  %>%
  mutate( left_group_id_start = (left_genomic_start %/% group_size)   + 1 )  %>%
  mutate( right_group_id_stop = (right_genomic_stop %/% group_size)   + 1 ) %>%
  mutate( right_group_id_start = (right_genomic_start %/% group_size)   + 1 )


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rna_type_tbl <- vroom::vroom( rna_type_file )


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
left_names_step_1.a <- hyb_tbl  %>%
  dplyr::filter( left_group_id_stop == left_group_id_start  ) %>%
  dplyr::select( row_id, left_chromo, left_genomic_start, left_genomic_stop, left_group_id_stop) %>%
  dplyr::left_join( rna_coordinates, by = c( "left_chromo" = "feature_chromosome",
                                              "left_group_id_stop" = "group_id") )  %>%
    dplyr::filter ( ( feature_start <= left_genomic_start & left_genomic_start <= feature_end ) |
                    ( feature_start <= left_genomic_stop & left_genomic_stop <= feature_start ) ) %>%
    left_join( rna_type_tbl, by=c("rna_type" = "rna_type"))  %>%
    dplyr::select( row_id, feature_name, rank, rna_type)

left_names_step_1.b <-  hyb_tbl %>%
  dplyr::filter( left_group_id_stop != left_group_id_start  )  %>%
  dplyr::select( row_id, left_chromo, left_genomic_start, left_genomic_stop) %>%
  dplyr::left_join( rna_coordinates, by = c( "left_chromo" = "feature_chromosome" ) )  %>%
    dplyr::filter ( ( feature_start <= left_genomic_start & left_genomic_start <= feature_end ) |
                    ( feature_start <= left_genomic_stop & left_genomic_stop <= feature_start ) ) %>%
    left_join( rna_type_tbl, by=c("rna_type" = "rna_type"))  %>%
    dplyr::select( row_id, feature_name, rank, rna_type)

left_names_step_1 <- left_names_step_1.a %>%
  bind_rows( left_names_step_1.b)

rm( left_names_step_1.a, left_names_step_1.b)
gc()

left_min_rank <- left_names_step_1 %>%
      dplyr::select(-feature_name) %>%
      group_by( row_id) %>%
      arrange( row_id, rank) %>%
      mutate( rank_row_id = row_number()) %>%
      ungroup() %>%
      filter( rank_row_id == 1)  

left_names_step_2 <- left_names_step_1 %>%
  dplyr::inner_join( left_min_rank, by =c("rank" = "rank", 
                                          "row_id" = "row_id",
                                          "rna_type" = "rna_type"))  %>%
  dplyr::select(-rank_row_id, -rank)  %>%
  dplyr::mutate(left_feature = str_split( feature_name, "\\|")  %>% 
                  purrr::map2_chr ( rna_type, 
                                    ~{ paste( c(.x[1:3], .y)  , collapse="|")}))
    
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
right_names_step_1.a <- hyb_tbl  %>%
  dplyr::filter( right_group_id_stop == right_group_id_start  ) %>%
  dplyr::select( row_id, right_chromo, right_genomic_start, right_genomic_stop, right_group_id_stop) %>%
  dplyr::left_join( rna_coordinates, by = c( "right_chromo" = "feature_chromosome",
                                              "right_group_id_stop" = "group_id") )  %>%
    dplyr::filter ( ( feature_start <= right_genomic_start & right_genomic_start <= feature_end ) |
                    ( feature_start <= right_genomic_stop & right_genomic_stop <= feature_start ) ) %>%
    left_join( rna_type_tbl, by=c("rna_type" = "rna_type"))  %>%
    dplyr::select( row_id, feature_name, rank, rna_type)

right_names_step_1.b <-  hyb_tbl %>%
  dplyr::filter( right_group_id_stop != right_group_id_start  )  %>%
  dplyr::select( row_id, right_chromo, right_genomic_start, right_genomic_stop) %>%
  dplyr::left_join( rna_coordinates, by = c( "right_chromo" = "feature_chromosome" ) )  %>%
    dplyr::filter ( ( feature_start <= right_genomic_start & right_genomic_start <= feature_end ) |
                    ( feature_start <= right_genomic_stop & right_genomic_stop <= feature_start ) ) %>%
    left_join( rna_type_tbl, by=c("rna_type" = "rna_type"))  %>%
    dplyr::select( row_id, feature_name, rank, rna_type)

right_names_step_1 <- right_names_step_1.a %>%
  bind_rows( right_names_step_1.b)

rm( right_names_step_1.a, right_names_step_1.b)
gc()



right_min_rank <- right_names_step_1 %>%
      dplyr::select(-feature_name) %>%
      group_by( row_id) %>%
      arrange( row_id, rank) %>%
      mutate( rank_row_id = row_number()) %>%
      ungroup() %>%
      filter( rank_row_id == 1)  

right_names_step_2 <- right_names_step_1 %>%
  dplyr::inner_join( right_min_rank, by =c("rank" = "rank", 
                                           "row_id" = "row_id", 
                                           "rna_type" = "rna_type"))  %>%
  dplyr::select(-rank_row_id, -rank) %>%
  dplyr::mutate(right_feature = str_split( feature_name, "\\|")  %>% 
                  purrr::map2_chr ( rna_type, 
                                    ~{ paste( c(.x[1:3], .y)  , collapse="|")}))



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
 hyb_tbl_merged <-  hyb_tbl %>%
  left_join( left_names_step_2, by=c( "row_id"))  %>% 
  left_join( right_names_step_2, by=c( "row_id")) %>%
  dplyr::select(-row_id)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# feature_overlap_partial <- partial( feature_overlap, feature_table = rna_coordinates, rna_type_priority=rna_type_tbl)
# 
# hyb_tbl_merged <-  hyb_tbl %>% 
#   mutate( left_feature = future_pmap_chr( list( left_chromo, left_genomic_start, left_genomic_stop), ~ feature_overlap_partial(..1, ..2, ..3)   ) ) %>%
#   mutate( right_feature = future_pmap_chr( list( right_chromo, right_genomic_start, right_genomic_stop), ~ feature_overlap_partial(..1, ..2, ..3)   ) ) 


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
hyb_tbl_final <- hyb_tbl_merged %>%
  dplyr::select(one_of( c(
"left_hyb_id",
"hyb_seq",
"binding_energy",
"left_feature",
"left_read_start", "left_read_stop",	"left_genomic_start",	"left_genomic_stop",	"left_strand"	,
"right_feature",
"right_read_start", "right_read_stop",	"right_genomic_start",	"right_genomic_stop",	"right_strand",
"optional_text" )))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
vroom::vroom_write( hyb_tbl_final, col_names=FALSE, file=output_file ) 

