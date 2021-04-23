## merge_hybrid_p_values.R
##
## copyright 2017 Chi Nam Ignatius Pang and Jai J. Tree 
## GNU Lesser General Public License
##
## merges hybrids that overlap between datasets and updates the p-value
## appends a field containing the number of experiments recovered
## updates the connection score with the read and hybrid information from both experments
## updates the hybrid count by adding overlapping expts together
## updates the max_read with the highest max_reads count.
##
## takes the output from calculate_generic_hybrid_p_values.R as inputs (two only) 
## usage: 
## 		merge_hybrid_p_values.R INPUT1.p_values.txt INPUT2.p_values.txt
## NB: must use R/3.2.2 or higher for dplyr.

library(sqldf)
library(dplyr)

## test if there is at least one argument: if not, return an error
args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("Two arguments must be supplied, both p_value tables", call.=FALSE)
}

## import dataset from command line and set constants

results_directory <- getwd()
hybA <- read.csv.sql(args[1], sep="\t")
hybB <- read.csv.sql(args[2], sep="\t")

p_value_threshold <- 0.05

print ("input files are: ")
print (args[1])
print (args[2])

#########################################################################################
## rename column headers in file1 and file2 - these must all be unique

hybA <- dplyr::rename( hybA,
							ID_a_r1=ID,
							rna_class_a_r1=rna_class,
							chromo_a_r1=chromo,
							name_a_r1=name,
							strand_a_r1=strand,
							ID_b_r1=ID.1,
							rna_class_b_r1=rna_class.1,
							chromo_b_r1=chromo.1,
							name_b_r1=name.1,
							strand_b_r1=strand.1,
							hyb_count_a_r1=hyb_count,
							hyb_count_b_r1=hyb_count.1,
							p_x_y_r1= p_x_y, 
							pdf_x_y_r1=pdf_x_y, 
							p_value_r1 = p_value,
							bh_adj_p_value_r1 = bh_adj_p_value,
							g_x_r1 = g_x,
							g_y_r1 = g_y,
							connection_score_r1= connection_score, 
							start_a_r1=start,
							end_a_r1=end,
							start_b_r1=start.1,
							end_b_r1=end.1,
							min_r1=min,
							max_r1=max,
							total_hybrids_r1=total_hybrids,
							cumulative_p_value_r1=cumulative_p_value,
							number_expts_r1=number_expts,
							cumulative_connection_score_r1=cumulative_connection_score)
									 								 	
hybB <- dplyr::rename( hybB,
							ID_a_r2=ID,
							rna_class_a_r2=rna_class,
							chromo_a_r2=chromo,
							name_a_r2=name,
							strand_a_r2=strand,
							ID_b_r2=ID.1,
							rna_class_b_r2=rna_class.1,
							chromo_b_r2=chromo.1,
							name_b_r2=name.1,
							strand_b_r2=strand.1,
							hyb_count_a_r2=hyb_count,
							hyb_count_b_r2=hyb_count.1,
							p_x_y_r2= p_x_y, 
							pdf_x_y_r2=pdf_x_y, 
		 					p_value_r2 = p_value,
							bh_adj_p_value_r2 = bh_adj_p_value,
							g_x_r2 = g_x,
							g_y_r2 = g_y,
							connection_score_r2= connection_score,									 	
							start_a_r2=start,
							end_a_r2=end,
							start_b_r2=start.1,
							end_b_r2=end.1,
							min_r2=min,
							max_r2=max,
							total_hybrids_r2=total_hybrids,
							cumulative_p_value_r2=cumulative_p_value,
							number_expts_r2=number_expts,
							cumulative_connection_score_r2=cumulative_connection_score)

print ("...table headers renamed...")

#########################################################################################
### find overlapping hybrids between datasets

## make indexes for each table
### Match, forward direction (Replicate 1, gene 1) == (Replicate 2, gene 1), (Replicate 1, gene 2) == (Replicate 2, gene 2)
### Match, reverse direction (Replicate 1, gene 1) == (Replicate 2, gene 2), (Replicate 1, gene 2) == (Replicate 2, gene 1)

merged_table <- sqldf( c( 
			'create index idx_start_a_r1 on hybA("start_a_r1")',
			'create index idx_end_a_r1   on hybA("end_a_r1")',
			'create index idx_start_b_r1 on hybA("start_b_r1")',
			'create index idx_end_b_r1   on hybA("end_b_r1")',	
			'create index idx_start_a_r2 on hybB("start_a_r2")',
			'create index idx_end_a_r2   on hybB("end_a_r2")',
			'create index idx_start_b_r2 on hybB("start_b_r2")',
			'create index idx_end_b_r2   on hybB("end_b_r2")', 
			'select * from hybA left outer join hybB on (
												  		(
															( start_a_r1 <= start_a_r2 and start_a_r2 <= end_a_r1 )
															or( start_a_r1 <= end_a_r2 and end_a_r2 <= end_a_r1 ) 
												  		) and (
															(start_b_r1 <= start_b_r2 and start_b_r2 <= end_b_r1 )
															or( start_b_r1 <= end_b_r2 and end_b_r2 <= end_b_r1 ) 
														)
											    		) or ( 
														(	
															(start_a_r1 <= start_b_r2 and start_b_r2 <= end_a_r1 )  
															or( start_a_r1 <= end_b_r2 and end_b_r2 <= end_a_r1 ) 
														) and (
															( start_b_r1 <= start_a_r2 and start_a_r2 <= end_b_r1 )
															or( start_b_r1 <= end_a_r2 and end_a_r2 <= end_b_r1 ) 
														)
													) union 
			select * from hybB left outer join hybA on (
												  		(
															( start_a_r1 <= start_a_r2 and start_a_r2 <= end_a_r1 )
															or( start_a_r1 <= end_a_r2 and end_a_r2 <= end_a_r1 ) 
												  		) and (
															(start_b_r1 <= start_b_r2 and start_b_r2 <= end_b_r1 )
															or( start_b_r1 <= end_b_r2 and end_b_r2 <= end_b_r1 ) 
														)
											    		) or ( 
														(	
															(start_a_r1 <= start_b_r2 and start_b_r2 <= end_a_r1 )  
															or( start_a_r1 <= end_b_r2 and end_b_r2 <= end_a_r1 ) 
														) and (
															( start_b_r1 <= start_a_r2 and start_a_r2 <= end_b_r1 )
															or( start_b_r1 <= end_a_r2 and end_a_r2 <= end_b_r1 ) 
														)
													)'
		) )

print ("...forward overlaps found...")
print ("...reverse overlaps found...")
							  	
## merged_table <- merged_table %>% mutate_if(is.factor, as.character) 
## merged_table_rev <- merged_table_rev %>% mutate_if(is.factor, as.character) 

##########################################################################################

## output tables
### Merge forward and reverse direction

print ("...merging forward and reverse tables...")

merged_table_edited <- dplyr::mutate ( merged_table, combined_bh_adj_p_value=1-pchisq(-2*(log(bh_adj_p_value_r1) + log(bh_adj_p_value_r2)), 4) ) # degree of freedom = 2 * number of samples = 4 in this case
write.table(merged_table_edited, paste(args[2], "merged_output_table.txt", sep="_"), sep="\t", quote=FALSE)
print (paste( args[2], "merged_output_table.txt", sep="_"))




