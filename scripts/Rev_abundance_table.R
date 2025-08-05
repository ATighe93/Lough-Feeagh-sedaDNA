library("phyloseq")
library("tidyverse")

### LIB1 ###
sd <- read.csv("track_rev_orient.csv", header = T, sep=",")
sd <- column_to_rownames(sd, "X")
sd <- phyloseq(sample_data(sd))
otu1 <- readRDS("seqtab_nochim2.rds")
otu1 <- phyloseq(otu_table(otu1, taxa_are_rows = F))
ps1 <- merge_phyloseq(otu1, sd)
saveRDS(ps1, "ps_d1.rds")

## MAKE SURE ALL OF THE FOLLOWING MATCH (ie, they all relate to the same seq run and orientation)
rdp_output="RDP_classified_rev.txt" #CHANGE ME to the path to your RDP output
dada_phyloseq="ps_d1.rds" #CHANGE ME to the path to your phyloseq object from the DADA2 run (combined metadata and seqtab_nochim file: R script: make_ps_object_for_parsing)
new_phyloseq="ps_d1_parsed.rds" #CHANGE ME to the path for your output file, call it whatever you want :)
threshold=0.6 #CHANGE ME to the confidence threshold you want to use for determining a successful assignment at a given taxonomic level

#load table and select the columns for the different taxonomic levels and their confidences
rdp=read.table(rdp_output,header=F,sep="\t")
rdp=rdp %>% dplyr::select(V1,V6,V8,V9,V11,V12,V14,V15,V17,V18,V20,V21,V23,V24,V26,V27,V29)
colnames(rdp)=c("ASVs","superkingdom","conf_superkingdom","kingdom","conf_kingdom","phylum","conf_phylum","class","conf_class","order","conf_order","family","conf_family","genus","conf_genus","species","conf_species")

# read in phyloseq object
physeq=readRDS(dada_phyloseq)

#parse the taxonomy to the lowest correctly assigned level
#note this is contained in an if loop to control that the rdp output and the phyloseq object are already in the correct order
#this script assumes that the names of the ASVs in the original phyloseq object match the names of the ASVs that were classified using the RDP classifier

if((sum(taxa_names(physeq)==rdp$ASVs)==nrow(rdp))=="TRUE"){
  taxtemp= rdp%>%dplyr::select(kingdom,phylum,class,order,family,genus,species)#make a temp object with just the taxonomic info that will be altered according to the thresholds
  taxtemp=sapply(taxtemp,function(x)as.character(x))
  for (i in (1:length(taxtemp[,1]))){
    if(rdp$conf_kingdom[i]<threshold){#labels things unassigned at kingdom level as unassigned
      taxtemp[,1][i]="NA"
      taxtemp[,2][i]="NA"
      taxtemp[,3][i]="NA"
      taxtemp[,4][i]="NA"
      taxtemp[,5][i]="NA"
      taxtemp[,6][i]="NA"
      taxtemp[,7][i]="Unassigned"
    }else{
      if(rdp$conf_kingdom[i]>threshold&rdp$conf_phylum[i]<threshold){ #labels kingdom sp.
        taxtemp[,2][i]="NA"
        taxtemp[,3][i]="NA"
        taxtemp[,4][i]="NA"
        taxtemp[,5][i]="NA"
        taxtemp[,6][i]="NA"
        taxtemp[,7][i]=paste0(rdp$kingdom[i]," sp.")
      }else{
        if(rdp$conf_phylum[i]>threshold-0.01&rdp$conf_class[i]<threshold){ #labels phylum sp.
          taxtemp[,3][i]="NA"
          taxtemp[,4][i]="NA"
          taxtemp[,5][i]="NA"
          taxtemp[,6][i]="NA"
          taxtemp[,7][i]=paste0(rdp$phylum[i]," sp.")
        }else{
          if(rdp$conf_class[i]>threshold-0.01&rdp$conf_order[i]<threshold){#labels class sp.
            taxtemp[,4][i]="NA"
            taxtemp[,5][i]="NA"
            taxtemp[,6][i]="NA"
            taxtemp[,7][i]=paste0(rdp$class[i]," sp.")
          }else{
            if(rdp$conf_order[i]>threshold-0.01&rdp$conf_family[i]<threshold){#labels order sp.
              taxtemp[,5][i]="NA"
              taxtemp[,6][i]="NA"
              taxtemp[,7][i]=paste0(rdp$order[i]," sp.")
            }else{
              if(rdp$conf_genus[i]<threshold & rdp$conf_family[i]>threshold-0.01){#labels family sp.
                taxtemp[,6][i]="NA"
                taxtemp[,7][i]=paste0(rdp$family[i]," sp.")
              }else{
                if(rdp$conf_species[i]<threshold & rdp$conf_genus[i]>threshold-0.01 & rdp$conf_family[i]>threshold-0.01){#labels genus sp.
                  taxtemp[,7][i]=paste0(rdp$genus[i], " sp.")
                }}}}}}}}}

rownames(taxtemp)=rdp$ASVs#label the rows in the new taxonomy table

if((sum(taxa_names(physeq)==rownames(taxtemp))==nrow(taxtemp))=="TRUE"){#make the new phyloseq object and save it
  tax_table(physeq)=taxtemp
  saveRDS(physeq,file=new_phyloseq)
} else {
  print("rownames don't match between phyloseq object and taxtemp table!")
}


ps <- readRDS("ps_d1_parsed.rds")
subset


g.ps <- subset_taxa(ps, !phylum=="NA")

tax <- as.data.frame(tax_table(g.ps))
tax$sums <- taxa_sums(g.ps)

DCC_otu <- as.data.frame(t(otu_table(g.ps)))

tax <- rownames_to_column(tax,"asv")
DCC_otu <- rownames_to_column(DCC_otu,"asv")
DCC_table_combo <- full_join(tax, DCC_otu)
write.csv(DCC_table_combo,"Rev_abundance_table.csv")

