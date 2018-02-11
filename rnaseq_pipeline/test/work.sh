cp ../data/expr.table ../data/normal.lst ../data/if.lst .
perl ../bin/run.pl GSE343141 expr.table normal.lst if.lst

ls |grep -v work | while read line; do rm -r $line;done;

