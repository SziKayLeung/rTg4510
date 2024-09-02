
cd /lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/2_cutadapt_merge/DN
for i in *_merged_combined.fasta*; do
  count=$(grep -c description "$i")
  echo -e "$i\t$count" >> read_numbers.txt
done

cd /lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/2_cutadapt_merge/NeuN
for i in *_merged_combined.fasta*; do
  count=$(grep -c description "$i")
  echo -e "$i\t$count" >> read_numbers.txt
done