#!/home/unhAW/jtmiller/.conda/envs/qiime2R/bin/Rscript --vanilla

export LD_LIBRARY_PATH='/usr/lib/jvm/java-11-openjdk-11.0.25.0.9-3.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
R 



## BLAST6
## conda activate qiime2R
#  

# options <- c('/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-seqs-derep-uniq_export/blast/mifishdb',
#     '/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-tax-derep-uniq.qza',
#     'NERRS-MF-Q1-Q5.rep-seqs.qza',
#     'NERRS_MF_Q1-Q5_vsearch-taxonomy-all-98.qza',
#     'NERRS_MF_Q1-Q5_BLAST6-tophit-all-98.qza',
#     'NERRS-MF-Q1-Q5.table.qza',
#     '/home/unhAW/jtmiller/watts/data/NERR/mifish/results',
#     'Q1-Q5',
#     'NERRS-MF-Q1-Q5_sklearn-taxonomy.qza',
#     'NERRS_SPECIES_LIST.txt',
#     130,
#     210,
#     'NERRS_MF_Q1-Q5_BLAST6-tophit-all-98_export/blast6.tsv',
#     'NERRS_SPECIES_LIST.txt')

#options <- c("/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-seqs-derep-uniq_export/blast/mifishdb",
#    "/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-tax-derep-uniq.qza",
#    "results/MBON_Sep21_rep-seqs.qza",
#    "results/MBON_Sep21_vsearch-taxonomy-all-95.qza",
#    "results/MBON_Sep21_BLAST6-tophit-all-95.qza",
#    "results/MBON_Sep21_table.qza",
#    "/home/unhAW/jtmiller/watts/data/MBON/mifish/results",
#    "MBON_Sep21",
#    "MBON_Sep21_sklearn-taxonomy.qza",
#    "results/N_Amer_Fish_Species.txt",
#    130,
#    210,
#    'results/MBON_Sep21_BLAST6-tophit-all-95_export/blast6.tsv')

options <- c("/home-wd/home/unhAW/jtmiller/watts/ref-database/MiFish/MitoFish/july2023/12S-seqs-derep-uniq_export/blast/mifishdb",
    "/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.tax.qza",
    "AW-Test-Chocorua-MFNX012425_rep-seqs.qza",
    "AW-Test-Chocorua-MFNX012425_vsearch-taxonomy-all-95.qza",
    "AW-Test-Chocorua-MFNX012425_BLAST6-tophit-all-95.qza",
    "AW-Test-Chocorua-MFNX012425_table.qza",
    "/home/users/jtm1171/Chocorua/mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out",
    "Chocorua",
    "classify-sklearn.qza",
    "N_Amer_Fish_Species.txt",
    "130",
    "210", 
    "AW-Test-Chocorua-MFNX012425_BLAST6-tophit-all-95_export/blast6.tsv")
################ BEGIN #######
options <- commandArgs(trailingOnly = TRUE)

db <- options[1]
reftax <- options[2]
seq <- options[3]
tax <- options[4]
top <- options[5]
stable <- options[6] 
out_dir <- options[7]
tag <- options[8]
sklearn <- options[9]
gbif <- options[10]
seq_min <- as.numeric(options[11])
seq_max <- as.numeric(options[12])
blast6 <- options[13]


require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)
#require(Biostrings)

source('/home/users/jtm1171/code/qiime2R_custom_functions.R')

seq <- read_qza(seq)$data
ASVs <- names(seq)

reftax <- read_qza(reftax)$data
rownames(reftax) <- reftax$Feature.ID

blast6 <- read.table(blast6, sep='\t')
colnames(blast6) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
blast6 <- blast6[which(blast6$qseqid %in% ASVs),]

gbif <- read.table(gbif, row.names = NULL, sep = '\t', header=F)
colnames(gbif) <- 'species'


ASVs <- unique(blast6$qseqid)
### top 5 blast6 hits
blast6 <- lapply(ASVs, function(id, input = blast6 ){
        indx <- which(input$qseqid == id)
        hits <- input[indx,][which( input[indx,'pident'] == max(input[indx,'pident'] ) ), ]
        
        species <- unique(parseMiFishTaxNames(reftax[hits$sseqid,'Taxon']))
        taxa <- paste(species, collapse = ';')
        pos_gbif <- paste(species[species %in% gbif$species], collapse = ';')

        seqid <- paste(hits$sseqid, collapse = ';')
        pident <- paste(unique(hits$pident), collapse = ';')

        return(data.frame(seqid,taxa,pident,pos_gbif))

})

blast6 <- data.frame(do.call(rbind,blast6))
rownames(blast6) <- ASVs 
write.table(blast6 , file = file.path(out_dir,paste0(tag,'_qiime_blast6_hits.txt')), sep = '\t', quote = F, row.names = T, col.names = T)
