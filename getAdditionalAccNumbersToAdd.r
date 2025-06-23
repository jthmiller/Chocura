## Need to add the mitohelper entries that do not have primer sequence. Pull out by lenght of ASV?
#### R 
library(tidyverse)
require(qiime2R)
require(taxize)
require(rBLAST)
require(RFLPtools)
require(Biostrings)

source('/home/users/jtm1171/code/qiimeras/R-functions.r')


#refs <- read_qza('/home/users/jtm1171/refdbs/mifish/june_2025/12S-seqs-derep-uniq.qza')$data
refs <- read_qza('/home/users/jtm1171/refdbs/mifish/june_2025/12S-16S-18S-seqs.qza')$data
refs <- names(refs)[which(lengths(refs) < 250)]


tls <- read.table('/home/users/jtm1171/refdbs/mifish/june_2025/12S-16S-18S_extracted.tax.list')$V1
#tls <- read.table('/home/users/jtm1171/refdbs/mifish/june_2025/extracted.tax.list')$V1
refs <- refs[!refs %in% tls] 

## Write IDs in refs that is not in the extracted list.  
write_lines(
  refs,
  file= '/home/users/jtm1171/refdbs/mifish/june_2025/12S-16S-18S_add_length_tax.txt',
  sep = "\n",
  na = "NA",
  append = FALSE
)