

##################################  Process genotypes in Vierstra et al ##################################################

VIERSTRA_DIR=/n/groups/price/kushal/ENCODE/data/Vierstra/Version1
cd $VIERSTRA_DIR
bcftools query -l genotypes.vcf.gz > sample_names.txt
cat sample_names.txt | wc -l

for line in `cat sample_names.txt | awk '{print $1}' | sort | uniq`;
do
samp_name=`echo $line | awk '{print $1}'`
cmd="bcftools query -s $samp_name -f '%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT\t%ARD\t%RD]\n' genotypes.vcf.gz > $samp_name.refalt.txt"
sbatch --time=30:00 --mem=20000 --output=vierstra.out --error=vierstra.err -p short -c 1 --wrap="$cmd"
echo "Finished processing for sample" $samp_name "\n"
done

##################################  Process genotypes in Vierstra et al ##################################################


VIERSTRA_DIR=/n/groups/price/kushal/ENCODE/data/Vierstra/Version2
cd $VIERSTRA_DIR
bcftools query -l allele_counts.vcf.gz > sample_names.txt
cat sample_names.txt | wc -l

for line in `cat sample_names.txt | awk '{print $1}' | sort | uniq`;
do
samp_name=`echo $line | awk '{print $1}'`
cmd="bcftools query -s $samp_name -f '%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%GT\t%ARD\t%RD]\n' allele_counts.vcf.gz > $samp_name.refalt.txt"
sbatch --time=30:00 --mem=20000 --output=vierstra.out --error=vierstra.err -p short -c 1 --wrap="$cmd"
echo "Finished processing for sample" $samp_name "\n"
done
