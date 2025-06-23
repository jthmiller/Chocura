#!/home/users/jtm1171/.conda/envs/qiime2-amplicon-2024.5/bin/Rscript --vanilla


sklearn_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_hybrid-taxonomy-all-95.qza'
blast_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425.blastout'
vsearch_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_vsearch-taxonomy-all-95.qza'
target <- 'Actinopteri'
tlevel <- 'Class'
feat_table <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_table.qza'
reps <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_rep-seqs.qza'
blast_names <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/blast-sci-com-names.tsv'
xlsx_out <- 'test.xlsx'


sklearn_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_hybrid-taxonomy-all-95.qza'
blast_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425.blastout'
vsearch_in <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_vsearch-taxonomy-all-95.qza'
target <- 'Actinopteri'
tlevel <- 'Class'
feat_table <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_table.qza'
reps <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_rep-seqs.qza'
blast_names <- '/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/blast-sci-com-names.tsv'
xlsx_out <- 'test.xlsx'

-rw-r--r--. 1 jtm1171 domain users  77K Jun 20 10:10 Chocorua_table.qza
-rw-r--r--. 1 jtm1171 domain users  74K Jun 20 10:10 Chocorua_rep-seqs.qza
-rw-r--r--. 1 jtm1171 domain users 6.6M Jun 20 10:22 mifish-2024.5-classifier.qza
-rw-r--r--. 1 jtm1171 domain users 109K Jun 20 10:24 Chocorua_hybrid_taxonomy.qza

### asv-output.r
library(tidyverse)
require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)
require(Biostrings)

source('/home/users/jtm1171/code/qiimeras/R-functions.r')

#sklearn_in <- options[1]  
#blast_in <- options[2]  
#vsearch_in <- options[3]  
#target <- options[4]  
#tlevel <- options[5]  
#feat_table <- options[6]  
#reps <- options[7]  


refs <- read_qza('/home/users/jtm1171/refdbs/mifish/june_2025/12S-seqs-derep-uniq.qza')$data

add_tax <- names(refs)[which(lengths(refs) < 250)]

write_lines(
  add_tax,
  file= '/home/users/jtm1171/refdbs/mifish/june_2025/add_tax.txt',
  sep = "\n",
  na = "NA",
  append = FALSE
)



feats <- read_qza('Chocorua_table.qza')$data
seqs <- read_qza('Chocorua_rep-seqs.qza')$data
tax <- read_qza('Chocorua_hybrid_taxonomy.qza')$data

tax <- read_qza('Chocorua_hybrid_taxonomy_old.qza')$data

rownames(tax) <- tax$Feature.ID

sort(feats[,'Ostrich-overnight'])

as.character(seqs['7aa794d91e82fea3d5d238cefc098359'])
tax['7aa794d91e82fea3d5d238cefc098359',]

blast_in <- read.csv(blast_in,row.names=1)
vsearch_in <- read_qza(vsearch_in)$data
sklearn_in <- read_qza(sklearn_in)$data
feat_table <- read_qza(feat_table)$data
reps <- read_qza(reps)$data

blast_names <- read.table(blast_names, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(blast_names) <- blast_names$sci

sklearn_in$Method <- 'SKLEARN_MF'
rownames(sklearn_in) <- sklearn_in$Feature.ID
sklearn_in <- sklearn_in[which(rownames(sklearn_in)%in% rownames(feat_table) ),]
sklearn_in[,'Taxon'] <- gsub("d__",'k__',sklearn_in[,'Taxon'])

ptax <- sklearn_in %>% parse_taxonomy()
sklearn_targets <- rownames(ptax)[which(ptax[,tlevel] == target)]

### Qiime Vesearch Taxonomy

vsearch_in$Method <- 'VSEARCH_MF'
rownames(vsearch_in) <- vsearch_in$Feature.ID
vsearch_in <- vsearch_in[which(rownames(vsearch_in)%in% rownames(feat_table) ),]
vsearch_in[,'Taxon'] <- gsub("d__",'k__',vsearch_in[,'Taxon'])
vsearch_in <- vsearch_in[which(!vsearch_in$Taxon == 'Unassigned'),]
pvtax <- vsearch_in %>% parse_taxonomy()
vsearch_targets <- rownames(pvtax)[which(pvtax[,tlevel] == target)]

#colnames(pvtax) <- paste0('VS_',colnames(pvtax))

############################################################################################
############################################################################################
############################################################################################


## Count levels for each
pvna <- rowSums(!is.na(pvtax))
pna <- rowSums(!is.na(ptax))

ASVs <- rownames(feat_table)
level_count <- sapply(ASVs,function(ASV){
 which.max( c(pna[ASV],pvna[ASV]) )
})
names(level_count) <- ASVs

use_pna <- names(which(level_count == 1))
use_pvna <- names(which(level_count == 2))

### Use update_ptax to edit the assignments
### Use vsearch or sklearn, whichever has the higher level of taxonomy
update_ptax <- pvtax[use_pvna,]
colnames(update_ptax) <- colnames(ptax)
update_ptax$Method <- 'VSEARCH_MF'
update_ptax <- data.frame(rbind(update_ptax, data.frame(ptax[use_pna, ],Method=sklearn_in[use_pna,'Method' ])))
update_ptax$Common <- blast_names[update_ptax$Species,'com']
##################################################################################################
target_ASVs <- unique(c(vsearch_targets,sklearn_targets))
offtarget_ASVs <- rownames(feat_table)[which(!rownames(feat_table) %in% target_ASVs)]
##################################################################################################


##################################################################################################
ASVs_targets_out <- data.frame(
    ASVID = target_ASVs,
    seq = as.character(reps[target_ASVs]),
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
    Blast_species_king = NA,
    Blast_species_com = NA,
    Blast_species_sci = NA,
    Blast_evalue = NA, 
    Blast_pident = NA,
    feat_table[target_ASVs,],
    check.names = FALSE
    )
rownames(ASVs_targets_out) <- target_ASVs
ASVs_targets_out[rownames(update_ptax), colnames(update_ptax) ] <- update_ptax

ASVs_targets_out[rownames(sklearn_in), c('Confidence','Consensus')] <- sklearn_in[, c('Confidence','Consensus')] 
ASVs_targets_out[, 'Total_reads'] <- rowSums(feat_table)[rownames(ASVs_targets_out)]
ASVs_targets_out[, 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(ASVs_targets_out)]
ASVs_targets_out[rownames(ASVs_targets_out), c('Blast_species_king','Blast_species_com','Blast_species_sci','Blast_evalue','Blast_pident')] <- blast_in[rownames(ASVs_targets_out), c('species_king','species_com','species_sci','evalue','pident')]


OTUs_targets_out <- add_otu_taxa(Updated_ASV = target_ASVs, tax = update_ptax, samps = colnames(feat_table), com = blast_names)



############################################################################################
############################################################################################


##################################################################################################
ASVs_offtarget_out <- data.frame(
    ASVID = offtarget_ASVs ,
    seq = as.character(reps[offtarget_ASVs ]),
    Total_reads = NA,
    Samples_Pos = NA,
    Blast_species_king = NA,
    Blast_species_com = NA,
    Blast_species_sci = NA,
    Blast_evalue = NA, 
    Blast_pident = NA,
    feat_table[offtarget_ASVs ,],
    check.names = FALSE
    )
rownames(ASVs_offtarget_out) <- offtarget_ASVs 
ASVs_offtarget_out[, 'Total_reads'] <- rowSums(feat_table)[rownames(ASVs_offtarget_out)]
ASVs_offtarget_out[, 'Samples_Pos'] <- rowSums(feat_table > 0)[rownames(ASVs_offtarget_out)]
ASVs_offtarget_out[rownames(ASVs_offtarget_out), c('Blast_species_king','Blast_species_com','Blast_species_sci','Blast_evalue','Blast_pident')] <- blast_in[rownames(ASVs_offtarget_out), c('species_king','species_com','species_sci','evalue','pident')]

cols <- c('Blast_species_king','Blast_species_com','Blast_species_sci')
orders <- order(apply(ASVs_offtarget_out[,cols],1,paste0,collapse=""))
ASVs_offtarget_out <- ASVs_offtarget_out[orders,]
##################################################################################################
##################################################################################################

############################################################################################
############################################################################################
############################################################################################
### OUTPUT XLSX #####
#####################
# out <- list(ASVs_out, otus_out,summary_out, tree_out, tip_out_hclust, tip_out_agnes, blast_out)
out <- list(
    ASVs_targets_out,
    OTUs_targets_out
    )
tax_cols <- colnames(out[[1]])[6:12]
out <- lapply(out, function(df, cols = tax_cols){  
    orders <- order(apply(df[,cols],1,paste0,collapse=""))
    new_out <- df[orders,]
    #print(head(new_out)[,1:20])
    return(new_out)
})
out <- c(out,list(ASVs_offtarget_out))

## xlsx_out <- paste0("MBON_JAN28",'BLAST_GBIF_CORRECTED.xlsx')
names(out) <- c("ASVs","OTUs","off-target")
openxlsx::write.xlsx(out, file = xlsx_out, rowNames=FALSE)
#openxlsx::write.xlsx(out, file = "", rowNames=FALSE)





### OUTPUT XLSX #####
#####################
############################################################################################
############################################################################################
############################################################################################
