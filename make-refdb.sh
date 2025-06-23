qiime feature-classifier extract-reads \
  --i-sequences 12S-seqs-derep-uniq.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-identity 0.65 \
  --o-reads 12S-seqs-derep-uniq_mifish_extracted_reads.qza \
  --p-n-jobs 32

qiime tools export --input-path 12S-seqs-derep-uniq_mifish_extracted_reads.qza --output-path 12S-seqs-derep-uniq_mifish_extracted_reads_export
## head ../*/add-dar-seqs.fasta
## Tax list of extracted reads and DAR
##
grep '^>' 12S-seqs-derep-uniq_mifish_extracted_reads_export/dna-sequences.fasta | sed 's/>//g' > extracted.tax.list
## cat tax.list add_tax.txt | sort | uniq > tax.list.sorted


### Full Mitoheper db, get tax for extracted reads 
qiime tools export --input-path 12S-tax-derep-uniq.qza --output-path 12S-tax-derep-uniq

while read acc ; do 
grep "^${acc}\s" 12S-tax-derep-uniq/taxonomy.tsv
done < extracted.tax.list > extracted.tax

## extracted fasta: 12S-seqs-derep-uniq_mifish_extracted_reads_export/dna-sequences.fasta
## extracted tax: extracted.tax

## DAR tax: DAR.tax
## DAR fasta: DAR.fasta


### Run R script to get lenght added ACC
while read acc ; do 
grep -A1 "^>${acc}$" 12S-seqs-derep-uniq/dna-sequences.fasta
done < add_length_tax.txt > lenght-seqs.fasta

while read acc ; do 
grep "^>${acc}$" 12S-tax-derep-uniq/taxonomy.tsv
done < add_length_tax.txt > lenght-seqs.tax

## Length fasta: lenght-seqs.fasta
## Length tax: lenght-seqs.tax

### DAR, Mitohelper, and mitohelper length seqs
cat \
  12S-seqs-derep-uniq_mifish_extracted_reads_export/dna-sequences.fasta  \
  DAR.fasta \
  lenght-seqs.fasta \
  > mitofish+DAR-65match.fasta

qiime tools import \
  --input-path mitofish+DAR-65match.fasta  \
  --output-path mitofish+DAR-65match.fasta \
  --type 'FeatureData[Sequence]'

awk -v FS='\t' -v OFS='\t' '{print $1,$2}' extracted.tax > mitofish+DAR-65.tax.ref
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' DAR.tax >> mitofish+DAR-65.tax.ref
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' lenght-seqs.tax >> mitofish+DAR-65.tax.ref

qiime tools import \
    --type "FeatureData[Taxonomy]" \
    --input-format "HeaderlessTSVTaxonomyFormat" \
    --input-path mitofish+DAR-65.tax.ref \
    --output-path mitofish+DAR-65.tax

# train without the weights
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads mitofish+DAR-65match.fasta.qza \
  --i-reference-taxonomy mitofish+DAR-65.tax.qza \
  --o-classifier mitofish+DAR-65-lenght-250_classifier.qza



awk -v FS='\t' -v OFS='\t' '{print $1,$2}' 12S-tax-derep-uniq/taxonomy.tsv > DAR+12S-tax-derep.tax
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' DAR.tax >> DAR+12S-tax-derep.tax

qiime tools import \
    --type "FeatureData[Taxonomy]" \
    --input-format "HeaderlessTSVTaxonomyFormat" \
    --input-path DAR+12S-tax-derep.tax \
    --output-path DAR+12S-tax-derep.tax



########################################################################################
########################################################################################
### 12s-18s-16s mifish
wget https://zenodo.org/records/15028392/files/12S-16S-18S-seqs.qza?download=1 -O 12S-16S-18S-seqs.qza
wget https://zenodo.org/records/15028392/files/12S-16S-18S-tax.qza?download=1 -O 12S-16S-18S-tax.qza

########################################################################################
########################################################################################

qiime feature-classifier extract-reads \
  --i-sequences 12S-16S-18S-seqs.qza \
  --p-f-primer GTCGGTAAAACTCGTGCCAGC \
  --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
  --p-identity 0.65 \
  --o-reads 12S-16S-18S-seqs_mifish_extracted-reads.qza \
  --p-n-jobs 12

grep '^>' 12S-16S-18S-seqs_mifish_extracted-reads/dna-sequences.fasta | sed 's/>//g' > 12S-16S-18S_extracted.tax.list

qiime tools export --input-path 12S-16S-18S-seqs_mifish_extracted-reads.qza --output-path 12S-16S-18S-seqs_mifish_extracted-reads
qiime tools export --input-path 12S-16S-18S-tax.qza --output-path 12S-16S-18S-tax
qiime tools export --input-path 12S-16S-18S-seqs.qza --output-path 12S-16S-18S-seqs

while read acc ; do 
grep "^${acc}\s" 12S-16S-18S-tax/taxonomy.tsv
done < 12S-16S-18S_extracted.tax.list > 12S-16S-18S_extracted.tax

## extracted fasta: 12S-16S-18S-seqs_mifish_extracted-reads/dna-sequences.fasta
## extracted tax: 12S-16S-18S_extracted.tax

############################################
### Run R script to get lenght added ACC ##
############################################
## output: 12S-16S-18S_add_length_tax.txt

while read acc ; do 
grep -A1 "^>${acc}$" 12S-16S-18S-seqs/dna-sequences.fasta
done < 12S-16S-18S_add_length_tax.txt > 12S-16S-18S_lenght-seqs.fasta

while read acc ; do 
grep "^${acc}\s" 12S-16S-18S-tax/taxonomy.tsv
done < 12S-16S-18S_add_length_tax.txt > 12S-16S-18S_lenght-seqs.tax

## testing
while read acc ; do 
grep "^${acc}\s" 12S-16S-18S-tax/taxonomy.tsv
done <<< "$(head 12S-16S-18S_add_length_tax.txt)" > 12S-16S-18S_lenght-seqs.tax


## Length fasta: 12S-16S-18S_lenght-seqs.fasta
## Length tax: 12S-16S-18S_lenght-seqs.tax

### DAR, Mitohelper, and mitohelper length seqs
cat \
  12S-16S-18S-seqs_mifish_extracted-reads/dna-sequences.fasta  \
  DAR.fasta \
  12S-16S-18S_lenght-seqs.fasta \
  > 12S-16S-18S+DAR-65match.fasta

qiime tools import \
  --input-path 12S-16S-18S+DAR-65match.fasta  \
  --output-path 12S-16S-18S+DAR-65match.fasta \
  --type 'FeatureData[Sequence]'

awk -v FS='\t' -v OFS='\t' '{print $1,$2}' 12S-16S-18S_extracted.tax > 12S-16S-18S+DAR-65.tax.ref
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' DAR.tax >> 12S-16S-18S+DAR-65.tax.ref
awk -v FS='\t' -v OFS='\t' '{print $1,$2}' 12S-16S-18S_lenght-seqs.tax >> 12S-16S-18S+DAR-65.tax.ref

qiime tools import \
    --type "FeatureData[Taxonomy]" \
    --input-format "HeaderlessTSVTaxonomyFormat" \
    --input-path 12S-16S-18S+DAR-65.tax.ref \
    --output-path 12S-16S-18S+DAR-65.tax






