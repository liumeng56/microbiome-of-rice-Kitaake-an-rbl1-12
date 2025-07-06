qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-n-threads 40 --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats stats.qza
qiime feature-classifier classify-sklearn --i-classifier /home/liumeng/rbl12-microbiome/silva-138-99-nb-classifier.qza --i-reads rep-seqs.qza --p-n-jobs 40 --o-classification taxonomy.qza
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast --p-include d__Bacteria --o-filtered-table feature-table-filt-contam.qza
qiime feature-table filter-features --i-table feature-table-filt-contam.qza --p-min-frequency 50 --p-min-samples 3 --o-filtered-table feature-table-final.qza
qiime feature-table filter-seqs --i-table feature-table-final.qza --i-data rep-seqs.qza --o-filtered-data rep-seqs-final.qza
qiime feature-table summarize --i-table feature-table-final.qza --o-visualization feature-table.qzv
qiime feature-table rarefy --i-table feature-table-final.qza --p-sampling-depth 32048 --o-rarefied-table table-rarefied.qza
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast --p-include d__Bacteria --o-filtered-table feature-table-filt-contam.qza
qiime feature-table filter-features --i-table feature-table-filt-contam.qza --p-min-frequency 50 --p-min-samples 3 --o-filtered-table feature-table-final.qza
qiime feature-table filter-seqs --i-table feature-table-final.qza --i-data rep-seqs.qza --o-filtered-data repset-seqs-final.qza
qiime tools export --input-path table-rarefied.qza --output-path feature-table
biom convert -i feature-table/feature-table.biom -o feature-table/feature-table.txt --to-tsv 
qiime tools export --input-path rep-seqs-final.qza --output-path rep-seqs
qiime tools export --input-path taxonomy.qza --output-path taxonomy
grep -v ";\_" taxonomy-level-6.txt > filtered_6_table.txt
sed 's/.*g__/g__/' filtered_6_table.txt > simplified_6_table.txt
awk '!/g__uncultured/' simplified_6_table.txt > final_6_table.txt
