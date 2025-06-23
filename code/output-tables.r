export LD_LIBRARY_PATH='/usr/lib/jvm/java-17-openjdk-17.0.15.0.6-2.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
R 


library(tidyverse)
library(qiime2R)
library(taxize)
require(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(rJava)
library(xlsx)
library(Biostrings)
require(rBLAST)
require(RFLPtools)
require(phyloseq)

options <- c(
  'mifish/results/Chocorua_table.qza',
  'mifish/results/Chocorua_hybrid_taxonomy_10.qza',
  'mifish/results/Chocorua_rep-seqs.qza',
  'Chocorua_output',
  'mifish/results/Chocorua_metadata.tsv',
  'mifish/results/Chocorua_rooted-tree.qza')

threads <- 12

source("/home/users/jtm1171/code/qiimeras/R-functions.r")

blasthits <- read.table('/home/users/jtm1171/Chocorua/mifish/results/12S-derep-tophits-blast/blast6.tsv', sep = '\t')
vtax <- read_qza('mifish/results/Chocorua_vsearch_taxonomy_all-accepts_90perc_12S-derep.qza')$data
fish <- read_qza('/home/users/jtm1171/refdbs/mifish/june_2025/DAR+12S-tax-derep.tax.qza')$data 
############################################################################################
############################################################################################
############################################################################################
# attach metadata
# metadata_out <- read.table(options[5], sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = "@")

############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
## plot sequences per sample
physeq <- qza_to_phyloseq(
  features = options[1], 
  taxonomy = options[2],
  tree = options[6],
  metadata = options[5]
  )
tax_table(physeq)[,'Kingdom'] <- gsub("tax=",'',tax_table(physeq)[,'Kingdom'])
tax_table(physeq)[,'Kingdom']  <- gsub("d__",'',tax_table(physeq)[,'Kingdom'])

metadata_out <- meta(physeq)
############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
## Feature table
feat_table <- read_qza(options[1])$data
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
## Rep seqs / ASVs
reps <- read_qza(options[3])$data
reps <- reps[rownames(feat_table)]
ASVs <- names(reps)
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
#### Hybrid Taxonomy
tax <- read_qza(options[2])$data
tax$Method <- paste0('hybrid_',tax$Method)
rownames(tax) <- tax$Feature.ID
tax <- tax[which(rownames(tax)%in% rownames(feat_table) ),]
tax[,'Taxon'] <- gsub("d__",'k__',tax[,'Taxon'])
ptax <- tax %>% parse_taxonomy()

### Qiime Vesearch Taxonomy
## vtax <- read_qza('mifish/results/Chocorua_vsearch_taxonomy_all-accepts_90perc_12S-derep.qza')$data
vtax$Method <- 'VSEARCH_MF'
rownames(vtax) <- vtax$Feature.ID
vtax <- vtax[which(rownames(vtax)%in% rownames(feat_table) ),]
vtax[,'Taxon'] <- gsub("d__",'k__',vtax[,'Taxon'])
vtax <- vtax[which(!vtax$Taxon == 'Unassigned'),]
pvtax <- vtax %>% parse_taxonomy()
colnames(pvtax) <- paste0('VS_',colnames(pvtax))
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
### VSEARCH BLAST HITS
### This subsets down to fish hits- not what we needed ### 
## blasthits <- read.table('/home/users/jtm1171/Chocorua/mifish/results/12S-derep-tophits-blast/blast6.tsv', sep = '\t')
 keep <- blasthits$V1[ which(blasthits$V3 > 0) ]
 blasthits <- blasthits[blasthits$V1 %in% keep,]
 physeq_fish <- prune_taxa(rownames(tax_table(physeq)) %in% keep, physeq) 
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################

## Count levels for each
pvna <- rowSums(!is.na(pvtax))
pna <- rowSums(!is.na(ptax))

level_count <- sapply(ASVs,function(ASV){
 which.max( c(pna[ASV],pvna[ASV]) )
})

use_pna <- names(which(level_count == 1))
use_pvna <- names(which(level_count == 2))

### Use update_ptax to edit the assignments
### Use vsearch or sklearn, whichever has the higher level of taxonomy
update_ptax <- pvtax[use_pvna,]
colnames(update_ptax) <- colnames(ptax)
update_ptax$Method <- 'VSEARCH_MF'
update_ptax <- data.frame(rbind(update_ptax, data.frame(ptax[use_pna, ],Method=tax[use_pna,'Method' ])))

#sum(ASVs %in% rownames(ptax))
#sum(ASVs %in% rownames(pvtax))
#ptax_no_id <- ASVs[!ASVs %in% rownames(ptax)[which(!is.na(ptax$Species))]]
#pvtax_yes_id <- ASVs[ASVs %in% rownames(pvtax)[which(!is.na(pvtax$VS_Species))]]
#update_ptax <- ptax_no_id[which(ptax_no_id %in% pvtax_yes_id )]
#update_ptax <- pvtax[update_ptax,]
#colnames(update_ptax) <- colnames(ptax)
#update_ptax$Method <- 'VSEARCH_MF'

### need to add in a search for ASVs that were partially id'd by VSEARCH (not just species)
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
## generate blastout
## blast full database
Sys.setenv(BLASTDB = '/home/share/databases/ncbi_nt' )
db <- blast('/home/share/databases/ncbi_nt/nt')

blast_outs <- predict(db, reps, BLAST_args=paste('-max_target_seqs 5 -num_threads',threads), custom_format='qseqid evalue bitscore pident sskingdoms sscinames scomnames',verbose = TRUE, keep_tmp = TRUE)
blast_outs <- fmtBlast(bo_in = blast_outs)
ASVs_blasted <- unique(blast_outs$qseqid) 
euk_asvs <- ASVs_blasted[!ASVs_blasted %in% unique(blast_outs[which(blast_outs$sskingdoms == "Bacteria" ),'qseqid'])]

blast_outs_20 <- predict(db, reps[euk_asvs], BLAST_args=paste('-max_target_seqs 20 -num_threads',threads), custom_format='qseqid evalue bitscore pident sskingdoms sscinames scomnames',verbose = TRUE, keep_tmp = TRUE)
blast_outs_20 <- fmtBlast(bo_in = blast_outs_20)
############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
## List for converting other sci to common from blast database
sci <- unique(blast_outs_20$sscinames)
com <- lapply(sci,function(X, blast_res = blast_outs_20 ){
    #unique(blast_res[which(blast_res$sscinames == X),'scomnames'])
    #blast_res[which(blast_res$sscinames == X),'scomnames'][1]
    unique(blast_res[which(blast_res$sscinames == X),'scomnames'])
})
com <- unlist(com)
names(com) <- sci
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
### Add common name to mitofish database
## fish <- read_qza('/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR.tax.qza')$data 
fish <- fish[-1,]
fish[,'Taxon']  <- gsub("^d__",'k__',fish[,'Taxon']) 
fish <- parse_taxonomy(fish)
uniq_reftax <- fish[!duplicated(fish$Species),] ##unique taxonomies 
rownames(uniq_reftax) <- uniq_reftax$Species
uniq_reftax$Common <- com[uniq_reftax$Species]
fish <- fish$Species
## add common names to current taxonomy update tables
ptax$Common <- com[ptax$Species]
pvtax$VS_Common <- com[pvtax$VS_Species]
update_ptax$Common <- com[update_ptax$Species]
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
### Get the top evalue for each asv
blast5 <- lapply(ASVs_blasted, bhit_sum, input = blast_outs)
blast5 <- data.frame(do.call(rbind,blast5))
rownames(blast5) <- ASVs_blasted

blast20 <- lapply(euk_asvs, bhit_sum, input = blast_outs_20 )
blast20 <- data.frame(do.call(rbind,blast20))
rownames(blast20) <- euk_asvs

indx <- which(!rownames(blast5) %in% rownames(blast20))
blast20 <- rbind(blast20, blast5[indx,])
############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
###### Taxonomy out
tax_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    row.names = ASVs
)
rownames(tax_out) <- ASVs

tax_out <- tax_out[rownames(tax_out)%in% rownames(feat_table),]

## Take the confidence and consensus from the qiime taxonomy
tax_out[rownames(tax), c('Confidence','Consensus')] <- tax[, c('Confidence','Consensus')] 
## Update tax_out with the updated taxonomy
tax_out[rownames(update_ptax), colnames(update_ptax) ] <- update_ptax

## tax_out has 
############################################################################################
############################################################################################
############################################################################################
## ASVs out with no blast corrections
ASVs_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    Total_reads = NA,
    Samples_Pos = NA,
    #Qiime_in_GBIF = NA,
    #BLAST_in_GBIF = NA,
    #BLAST_GBIF_GUESS = NA,
    #BLAST_GBIF_GUESS_in_refdb = NA,
    feat_table[ASVs,],
    seq = as.character(reps[ASVs]),
    check.names = FALSE
    )
rownames(ASVs_out) <- ASVs
ASVs_out[rownames(update_ptax), colnames(update_ptax) ] <- update_ptax
ASVs_out[rownames(tax), c('Confidence','Consensus')] <- tax[, c('Confidence','Consensus')] 
ASVs_out[, 'Total_reads'] <- rowSums(feat_table)[rownames(ASVs_out)]
ASVs_out[, 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(ASVs_out)]
############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################
# Species table
## Collapse taxa (should not have sequences- multiple sequences per taxa)
otus <- tax_glom(physeq, taxrank=rank_names(physeq)[7], NArm=F, bad_empty=c(NA, "", " ", "\t"))

OTUs_out <- data.frame(
  ASV_ID = rownames(tax_table(otus)),
  tax_table(otus),
  Common = com[tax_table(otus)[,'Species']],
  num_samples_pos = rowSums(otu_table(otus)>0),
  total_reads = rowSums(otu_table(otus)),
  otu_table(otus),
  check.names = FALSE
  )
OTUs_out <- OTUs_out[,-1]
############################################################################################
############################################################################################
############################################################################################

### HERE 

############################################################################################
############################################################################################
############################################################################################
Updates_table <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    Total_reads = NA,
    Samples_Pos = NA,
    ## Qiime_in_GBIF = NA,
    ## BLAST_in_GBIF = NA,
    ## BLAST_GBIF_GUESS = NA,
    ## BLAST_GBIF_GUESS_in_refdb = NA,
    blast20[ASVs,],
    seq = as.character(reps[ASVs]),
    check.names = FALSE,
    row.names = ASVs
    )

### ADD qiime taxonomy before updates
indx <- rownames(Updates_table)[rownames(Updates_table) %in% rownames(ptax)]
Updates_table[indx, colnames(ptax)] <- ptax[indx,]

indx <- rownames(Updates_table)[rownames(Updates_table) %in% rownames(tax)]
Updates_table[indx,c('Confidence','Consensus','Method') ] <- tax[indx,c('Confidence','Consensus','Method')]

Updates_table[rownames(Updates_table), 'Total_reads'] <- rowSums(feat_table)[rownames(Updates_table)]
Updates_table[rownames(Updates_table), 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(Updates_table)]



## DONT USE CHOCURA

############################################################################################
############################################################################################
############################################################################################
## Blast Output with all the info for corrections
blast_out <- data.frame(
    ASVID = ASVs,
    Confidence = NA, 
    Consensus = NA,  
    Method = NA,
    Kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = NA,
    Common = NA,
    species_king = NA,
    species_com = NA, 
    species_sci = NA,
    evalue = NA,       
    pident = NA,
    is.fish = NA,
    Total_reads = NA,
    Samples_Pos = NA,
    seq = as.character(reps[ASVs])
    )

### Uncorrected taxonomy from sklearn and vsearch (update_ptax is before gbif corrections)
blast_out[rownames(update_ptax), colnames(update_ptax)] <- update_ptax
blast_out[rownames(tax_out),c('Confidence','Consensus','Method') ] <- tax_out[,c('Confidence','Consensus','Method')]
blast_out[names(rowSums(feat_table)), 'Total_reads'] <- rowSums(feat_table)
blast_out[names(rowSums(feat_table)), 'Samples_Pos'] <- rowSums(feat_table > 0)

tax_blast <- blast_out

blast_out[ rownames(blast20),colnames(blast20)] <- blast20
blast_out <- blast_out[rownames(blast_out) %in% rownames(feat_table), ] 
blast_out <- data.frame(cbind(blast_out,feat_table[blast_out$ASVID,]),check.names = FALSE)







##### ID fish
bhits_fish <- rownames(blast_out)[grep('YES',blast_out$is.fish )]
blast_fish <- unique(c(rownames(ptax),rownames(pvtax),bhits_fish))

## non-fish blast out
non_fish <- blast_out[!rownames(blast_out) %in% blast_fish, ]
non_fish <- non_fish[!is.na(non_fish$pident),]

getcol <- c("ASVID",colnames(non_fish)[13:length(non_fish[1,])])
non_fish <- non_fish[,getcol]
non_fish <- non_fish[order(non_fish$species_king,non_fish$species_sci),]

## only fish blast out
blast_out <- blast_out[rownames(blast_out) %in% blast_fish, ]

### redefining tax_out here (why?)
tax_out2 <- data.frame(cbind(tax_blast, pvtax[tax_blast$ASVID,], blast20[tax_blast$ASVID,]), check.names = FALSE)

## only keep ptax, pvtax, and is.fish results in final blast table
############################################################################################
############################################################################################
############################################################################################



############################################################################################
############################################################################################
############################################################################################


############################################################################################
############################################################################################
############################################################################################


##Updated_ASVs

metadata_in = meta(physeq)
ASVs_in = ASVs_out
nfish = non_fish
blast = blast_out

## metadata_out[grep('TB',samps_in,value=T),'sample.blank'] <- 'BLANK'

all_in <- rownames(metadata_in)
asvs_out_all <- ASVs_in[rowSums(ASVs_in[,all_in],na.rm = T)>0,]
otus_out_all <- add_tax(Updated_ASV = ASVs_in, samps = all_in)

############################################################################################
############################################################################################
############################################################################################




############################################################################################
############################################################################################
############################################################################################
### OUTPUT XLSX #####
#####################
# out <- list(ASVs_out, otus_out,summary_out, tree_out, tip_out_hclust, tip_out_agnes, blast_out)
out <- list(
    asvs_out_all,
    otus_out_all
    )

tax_cols <- colnames(out[[1]])[5:11]
out <- lapply(out, function(df, cols = tax_cols){  
    orders <- order(apply(df[,cols],1,paste0,collapse=""))

    new_out <- df[orders,]
    #print(head(new_out)[,1:20])
    return(new_out)
})

out <- c(out,list(non_fish, metadata_out))

xlsx_out <- paste0("Chocura",'_BLAST.xlsx')

names(out) <- c("ASVs_sklearn","OTUs_sklearn","non_fish","metadata")
openxlsx::write.xlsx(out, file = xlsx_out, rowNames=FALSE)





### OUTPUT XLSX #####
#####################
############################################################################################
############################################################################################
############################################################################################



## ADD 