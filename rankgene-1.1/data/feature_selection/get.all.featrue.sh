for ((i=1;i<=$n;i++)); do gene_name=`sed -n "$i"p svm.map | cut -f2`; cat ../../../GPL570-13270.txt | grep "$gene_name" | cut -f1 >> temp; echo "\n" >>temp; done
