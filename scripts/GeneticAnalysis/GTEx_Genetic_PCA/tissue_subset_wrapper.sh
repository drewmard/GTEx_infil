while read -r TISSUE
do
	echo "${TISSUE}"
	qsub /home/anm2868/scripts/GTEx_Genetic_PCA/PCA_calc.sh $TISSUE
done < "/home/anm2868/etc/TISSUES.txt"
