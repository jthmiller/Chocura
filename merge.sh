
## merge table
qiime feature-table merge \
  --i-tables mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_table.qza \
  --i-tables mifish/runs/Chocorua-MFNX052825-90/qiime_out/Chocorua-MFNX052825-90_table.qza \
  --i-tables mifish/runs/Chocorua-MFNX052825-ON/qiime_out/Chocorua-MFNX052825-ON_table.qza \
  --o-merged-table mifish/results/Chocorua_table.qza

## merge asvs
qiime feature-table merge-seqs \
  --i-data mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_rep-seqs.qza \
  --i-data mifish/runs/Chocorua-MFNX052825-90/qiime_out/Chocorua-MFNX052825-90_rep-seqs.qza \
  --i-data mifish/runs/Chocorua-MFNX052825-ON/qiime_out/Chocorua-MFNX052825-ON_rep-seqs.qza \
  --o-merged-data mifish/results/Chocorua_rep-seqs.qza



    ## taxonomy
    maxaccepts=all
    query_cov=0.90 
    perc_identity=0.90
    weak_id=0.80 
    tophit_perc_identity=0.90

    refreads=${refreads:-/home/users/jtm1171/refdbs/mifish/july2023/12S-seqs-derep-uniq.qza}
    reftax=${reftax:-/home/users/jtm1171/refdbs/mifish/july2023/12S-tax-derep-uniq.qza}
    blastdb=${blastdb:-/home/share/databases/ncbi_nt/nt}
    sklearn=${sklearn:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish-classifier.qza}


    echo "running mifsh"
    echo "refreads: $refreads"
    echo "reftax: $reftax"
    echo "sklearn: $sklearn"
    echo "blastdb: $blastdb"

    min=99
    max=251



refreads=${refreads:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.fasta.qza}
reftax=${reftax:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.tax.qza}
sklearn=${sklearn:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed-classifier.qza}

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $refreads \
  --i-reference-taxonomy $reftax \
  --o-classifier mifish/results/mifish-2024.5-classifier.qza

### sklearn hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query mifish/results/Chocorua_rep-seqs.qza \
  --i-classifier mifish/results/mifish-2024.5-classifier.qza \
  --i-reference-reads /home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.fasta.qza \
  --i-reference-taxonomy /home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0 \
  --o-classification mifish/results/Chocorua_hybrid_taxonomy


### Generate taxonomy tables
conda activate qiime2R

head -n2 mifish/runs/AW-Test-Chocorua-MFNX012425/qiime_out/AW-Test-Chocorua-MFNX012425_dns_export/metadata.tsv > mifish/results/Chocorua_metadata.tsv
find . -name "metadata.tsv" -exec cat {} \; | grep -v "sample-id" | grep -v "#q2:types" >> mifish/results/Chocorua_metadata.tsv

qiime2_output_tables.r mifish/results/Chocorua_table.qza mifish/results/Chocorua_hybrid_taxonomy.qza mifish/results/Chocorua_rep-seqs.qza mifish/results/Chocorua_ASVs.csv mifish/results/Chocorua_metadata.tsv


refreads=${refreads:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.fasta.qza}
reftax=${reftax:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed.tax.qza}
sklearn=${sklearn:-/home/users/jtm1171/refdbs/mifish/july2023/mitofish+DAR+relaxed-classifier.qza}

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /home/users/jtm1171/refdbs/mifish/july2023/12S-seqs-derep-uniq.qza \
  --i-reference-taxonomy /home/users/jtm1171/refdbs/mifish/july2023/12S-tax-derep-uniq.qza \
  --o-classifier mifish/results/mifish-2024.5-classifier-old.qza



### sklearn hybrid taxonomy
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query mifish/results/Chocorua_rep-seqs.qza \
  --i-classifier /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65-lenght-250_classifier.qza \
  --i-reference-reads /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65match.fasta.qza \
  --i-reference-taxonomy /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65.tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts 10 \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0 \
  --o-classification mifish/results/Chocorua_hybrid_taxonomy_10


export LD_LIBRARY_PATH='/usr/lib/jvm/java-17-openjdk-17.0.15.0.6-2.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
export PATH="$PATH:~/code"
qiime2_output_tables.r mifish/results/Chocorua_table.qza mifish/results/Chocorua_hybrid_taxonomy_10.qza mifish/results/Chocorua_rep-seqs.qza mifish/results/Chocorua_ASVs_10 mifish/results/Chocorua_metadata.tsv


qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query mifish/results/Chocorua_rep-seqs.qza \
  --i-classifier /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65-lenght-250_classifier.qza \
  --i-reference-reads /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65match.fasta.qza \
  --i-reference-taxonomy /home/users/jtm1171/refdbs/mifish/june_2025/mitofish+DAR-65.tax.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
  --p-confidence 0 \
  --o-classification mifish/results/Chocorua_hybrid_taxonomy_all

export LD_LIBRARY_PATH='/usr/lib/jvm/java-17-openjdk-17.0.15.0.6-2.el9.x86_64/lib/server:$LD_LIBRARY_PATH'
export PATH="$PATH:~/code"
qiime2_output_tables.r mifish/results/Chocorua_table.qza mifish/results/Chocorua_hybrid_taxonomy_all.qza mifish/results/Chocorua_rep-seqs.qza mifish/results/Chocorua_ASVs_all mifish/results/Chocorua_metadata.tsv







qiime feature-classifier classify-consensus-vsearch \
    --i-query mifish/results/Chocorua_rep-seqs.qza \
    --i-reference-reads /home/users/jtm1171/refdbs/mifish/june_2025/12S-seqs-derep-uniq.qza \
    --i-reference-taxonomy /home/users/jtm1171/refdbs/mifish/june_2025/12S-tax-derep-uniq.qza \
  --p-threads 8 \
  --p-query-cov 0.95 \
  --p-perc-identity 0.90 \
  --p-maxrejects all \
  --p-maxaccepts all \
  --p-maxhits all \
  --p-min-consensus 0.51 \
    --o-classification mifish/results/Chocorua_vsearch_taxonomy_all-accepts_90perc_12S-derep \
    --o-search-results mifish/results/Chocorua_vsearch_taxonomy_all-accepts_90perc_12S-derep-tophits


qiime phylogeny align-to-tree-mafft-fasttree \
   --i-sequences mifish/results/Chocorua_rep-seqs.qza \
   --o-alignment mifish/results/Chocorua_aligned-rep-seqs \
   --o-masked-alignment mifish/results/Chocorua_masked-aligned-rep-seqs.qza \
   --o-tree mifish/results/Chocorua_unrooted-tree.qza \
   --o-rooted-tree mifish/results/Chocorua_rooted-tree.qza \
   --p-n-threads 24