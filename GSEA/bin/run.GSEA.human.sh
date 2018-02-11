#!/bin/sh

if [ $# -lt 2 ];then
	echo ""
	echo "	Pipeline for running GSEA"
	echo "	Author: zhoujj2013@gmail.com";
	echo "	usage: sh $0 <expr.foldchange.rnk> <prefix>";
	echo ""
	exit;
fi

gsea=$(dirname $0)/gsea-3.0.jar
expr=$1
prefix=$2

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.H -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.H -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.mir.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.mir -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.mir.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.mir -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.tft.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.tf -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c3.tft.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.tf -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.C5 -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.all.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.C5 -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.bp -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.bp -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

echo java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.kegg -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false

java -Xmx512m -cp $gsea xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $expr -scoring_scheme weighted -rpt_label $prefix.kegg -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ./ -gui false
