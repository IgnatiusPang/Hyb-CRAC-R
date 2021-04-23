## calculate_hybrid_p_values.R
## 
## copyright 2017 Chi Nam Ignatius Pang and Jai J. Tree 
## GNU Lesser General Public License
##
## usage: 
##	 	calculate_generic_hybrid_p_values.R MAX_READS.txt X
## where: 
## 		MAX_READS.txt 	= the output from max_reads.sh (maximum CRAC reads at hybrid half)
## 		X 				= the total number of mapped reads in CRAC experiment

## NB: must use R/3.2.2 or higher for dplyr.

if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(tidyverse)
p_load(vroom)

## test if there is at least one argument: if not, return an error
args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("Two arguments must be supplied, MAX_READS.txt and the total number of reads in dataset", call.=FALSE)
}

print ('input file:')
print (args[1])
print ('number of reads:')
print (args[2])

## import dataset from command line and set constants

results_directory <- getwd()
replicate_1 <- vroom::vroom(args[1])
replicate_1_total_reads <- as.numeric(args[2])

# replicate_1 <- vroom::vroom("/home/ignatius/PostDoc/2020/StaphCLASH2020/Results/Hybrids_Annot/DM_MiSeq/all_NNNCACTAGC_L5Bc_rnc-HTF_1.hyb.max_reads.txt")
# replicate_1_total_reads <- 897

p_value_threshold <- 0.05

## Rename left_max_reads to g_x, right_max_reads to g_y
## g_x = Reads Per Million (RPM) for g_x
## g_y = RPM for g_y

if( "left_max_reads" %in% colnames(replicate_1)) {
	replicate_1 <- dplyr::rename( replicate_1,  g_x= "left_max_reads")
}

if( "right_max_reads" %in% colnames(replicate_1)) {
	replicate_1 <- dplyr::rename( replicate_1,   g_y= "right_max_reads" )
}

##############################
### 	Functions          ###
##############################

### Calculate the p-value
## Adds the columns p_x_y, pdf_x_y and bh_adj_p_value to the table

calculate_adjusted_p_value <- function(input_table, total_reads) { 

	input_table <- dplyr::mutate( input_table, g_x=as.numeric(g_x), g_y=as.numeric(g_y))
	
	## Calculate a proxy for the probability P(g_x, g_y)  =  2 * (g_x RPM) * (g_y RPM)
	input_table <- dplyr::mutate( input_table,  p_x_y=2*g_x*g_y )
	
	## Calculate the normalized value 
	total_for_normalization <- sum(input_table[, "p_x_y"] )
	
	input_table <- input_table %>% 
	                dplyr::mutate( pdf_x_y = p_x_y/total_for_normalization ) %>%
	
	                ### Calculate the binomial probability
                  dplyr::mutate(  p_value = dbinom( left_hyb_count, size= total_reads, pdf_x_y )) %>%
	  
	                ### Calculate the Benjamini-Hochberg adjusted p-value
                  dplyr::mutate(  bh_adj_p_value = p.adjust( p_value, method="BH" ) )
	
	return(input_table)
} 

### Calculate the connection score as described by Chang et al. in the PARIS paper
calculate_connection_score <- function ( input_table, hyb_count_column, g_x_column, g_y_column) {

  connection_score_vector <- as.vector(input_table[,hyb_count_column] /
              sqrt( input_table[,g_x_column] * input_table[,g_y_column] ) )  

	updated_column <- data.frame( connection_score_vector )
	
 colnames( updated_column) <- c("connection_score")

	output_table <- input_table %>% bind_cols( updated_column)

	return( output_table)

}

##############################


print ("file formatted ok")
replicate_1_edited <- calculate_adjusted_p_value(replicate_1, replicate_1_total_reads) 
print ("adjusted p value colculated")
replicate_1_edited_v2 <- calculate_connection_score(replicate_1_edited, "left_hyb_count", "g_x", "g_y")
print ("connection score calculated")

# add fields for total hybrid count, cumulative p-value, number of experiments

# colnames( replicate_1_edited_v2)


updated_replicate_1_edited <- replicate_1_edited_v2 %>% 
            dplyr::mutate( total_hybrids= purrr::map_dbl( left_hyb_count, ~.) ) %>% 
            dplyr::mutate( cumulative_p_value=   purrr::map_dbl( bh_adj_p_value, ~.)  ) %>% 
            dplyr::mutate( cumulative_connection_score=  purrr::map_dbl( connection_score, ~.) ) %>% 
            dplyr::mutate( number_expts=1) %>% 
            dplyr::select( left_hyb_id, g_x, left_hyb_count, left_rna_class, left_chromo, 
                           left_name, left_genomic_start, left_genomic_stop, left_strand,
                           right_hyb_id, g_y, right_hyb_count, right_rna_class, right_chromo, 
                           right_name, right_genomic_start, right_genomic_stop, right_strand,
                           smaller_max_depth, larger_max_depth, 
                           p_x_y,	pdf_x_y,	p_value,	bh_adj_p_value,	connection_score,	total_hybrids,	
                           cumulative_p_value,	cumulative_connection_score,	number_expts,
                           two_way_merged,
                           left_binding_energy,
                           left_transcript_start,
                           left_transcript_stop,   
                           right_binding_energy,
                           right_transcript_start,
                           right_transcript_stop,
                           hyb_seq )  
  
## write the hyb file with p-values to file
vroom::vroom_write (updated_replicate_1_edited, path=file.path( paste(args[1],".p_values.txt", sep=""))) 									 	
									 	