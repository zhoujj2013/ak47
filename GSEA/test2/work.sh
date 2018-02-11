Enrichment analysis with user defined geneset.

less -S GSM1908737_cuffdiff.gene.level_copy2.txt | cut -f 1,10 | grep -v "^#" | sort -k2gr | awk '$2<100 && $2 > -100' > gene.rnk
