#!/bin/sh

# Written by Manik Uppal (2018-2019). Contact: mdu2002@med.cornell.edu


#this function is to calculate the proper factor to subsample each bam file; it will take 3 arguments - 1) the normal tissue library size, 2) sample (normal or immune to keep it general) library size, 3) sample percentage as specified in the input arguments
#--> General Formula: DesiredImm% x NormLibSize = ImmLibSize x NecImm%
# alias bamtools=/Users/ManikUppal/bamtools/bamtools

generateFactor () {
bc -l <<< "($1/$2)*$3"
}
calc ()
{
bc -l <<< "$@"
}


IFS=' ' read normal nmlpct <<< $1
#echo $normal
#echo $nmlpct
NumNormalReads=$(samtools idxstats $normal | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}')
#echo $NumNormalReads

finalOutputString=""
for var in "$@"
do
    IFS=' ' read var1 var2 <<< "$var"
    immune=$var1
    immunepct=$var2
    NumImmuneReads=$(samtools idxstats $immune | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}')
#echo $NumImmuneReads
    FACTOR=$(generateFactor $NumNormalReads $NumImmuneReads $immunepct)
#echo $FACTOR
    halfName=$(echo $immune | sed 's/_.*//')
    fullName=$(echo ${halfName//*\/})
    ImmName=$(calc $immunepct*100)
    ImmNameInt=${ImmName%.*}
    intOutputName="${fullName}_${ImmName}pct.intermediate.bam"
    finalOutputString+="${fullName}_${ImmName}pct."
#echo $outputName
    echo "Subsampling ${fullName} to ${ImmName} percent..."
    samtools view -s $FACTOR -b $immune -o $intOutputName
done
finalOutputName="${finalOutputString}bam"
finalOutputSortedName="${finalOutputString}sorted.bam"
echo $finalOutputName

ls *.intermediate.bam > files.bamlist
echo "Merging intermediate bam files..."
#bamtools used here because of proper handling of headers in BAM/SAM files for merging
bamtools merge -list files.bamlist -out $finalOutputName
# samtools merge -list files.bamlist -out $finalOutputName
rm *intermediate.bam
samtools sort $finalOutputName -o $finalOutputSortedName
rm $finalOutputName
samtools index $finalOutputSortedName

#include rm *intermediate.bam and rm $finalOutputName

# move synthetic mix file
mv ${finalOutputString}* /athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/SynMix
