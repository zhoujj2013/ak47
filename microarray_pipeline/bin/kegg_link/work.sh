grep GOTERM up.input.tableReport.final.txt |  cut -f 4,12 | awk -F"\t" '{print $2"\t"$1}' | head -10  > GO.txt

/x400ifs-accel/zhoujj/software/R-3.2.0/bin/Rscript anno.barplot.r
convert -density 150 test.pdf -quality 90 test.png

less -S human.genenames.txt | cut -f 2,9 | grep -v "Approved Symbol" > human.geneid.lst

less -S MGI_Gene_Model_Coord.rpt | cut -f 3,6 | grep -v "3. marker symbol" > mouse.geneid.lst

perl getPATHWAYmap.pl up.input.tableReport.final.txt mouse.geneid.lst > kegg.link

perl getPATHWAYmap.pl up.input.tableReport.final.txt down.input.tableReport.final.txt human.geneid.lst > kegg.link.txt

