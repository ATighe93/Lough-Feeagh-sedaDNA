library("tidyverse")

DCC_table_combo <- read.csv("Fwd_abundance_table.csv", header = T, sep=",")

DCC_table_combo <- DCC_table_combo %>% mutate(seq_id=row_number())
blast <- read.csv("ASV_fwd_blast_top_hit.csv", header = FALSE, sep="\t")
blast=blast %>% dplyr::select(V1,V2,V3,V4,V5)
colnames(blast)=c("seq_id","Blast_hit","Percent_ID","Align_length","e-value")


My_table_combo <- full_join(DCC_table_combo, blast, by="seq_id")
write.csv(My_table_combo,"Fwd_combined_results.csv")
