cd ./step7/cmap
mkdir ./filter_t1/
cat output_t1/*txt | awk -F "\t" '{if($6<0.05){print $0}}' > filter_t1/filter_pvalue.txt
le filter_t1/filter_pvalue.txt | awk -F "\t" '{print $NF}' | sort | uniq > filter_t1/drug.txt
sed -i "s/'//g" filter_t1/drug.txt
sed -i "s/\[//g" filter_t1/drug.txt
sed -i "s/\]//g" filter_t1/drug.txt
cat filter_t1/drug.txt | while read line; do a=`grep $line filter_t3/filter_pvalue.txt |wc -l`; echo -e $a"\t"$line >> filter_t1/drug_num.txt; done
grep "^2" filter_t1/drug_num.txt | awk -F "\t" '{print $2}' > filter_t1/drug_num_2.txt #将up和down的p值都小于0.05的药物筛选出来
cat filter_t1/drug_num_2.txt | while read line; do grep $line filter_t1/filter_pvalue.txt >> filter_t1/filter_pvalue_2.txt; done

mkdir ./filter_t3/
cat output_t3/*txt | awk -F "\t" '{if($6<0.05){print $0}}' > filter_t3/filter_pvalue.txt
le filter_t3/filter_pvalue.txt | awk -F "\t" '{print $NF}' | sort | uniq > filter_t3/drug.txt
sed -i "s/'//g" filter_t3/drug.txt
sed -i "s/\[//g" filter_t3/drug.txt
sed -i "s/\]//g" filter_t3/drug.txt
cat filter_t3/drug.txt | while read line; do a=`grep $line filter_t3/filter_pvalue.txt |wc -l`; echo -e $a"\t"$line >> filter_t3/drug_num.txt; done
grep "^2" filter_t3/drug_num.txt | awk -F "\t" '{print $2}' > filter_t3/drug_num_2.txt #将up和down的p值都小于0.05的药物筛选出来
cat filter_t3/drug_num_2.txt | while read line; do grep $line filter_t3/filter_pvalue.txt >> filter_t3/filter_pvalue_2.txt; done

#le compoundinfo_beta.txt | awk -F "\t" '{print $1"\t"$2"\t"$7}' > compoundinfo_beta.txt.used
