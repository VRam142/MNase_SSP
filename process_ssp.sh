#Wrapper script for processing MNase-SSP data
#VARS
r1=$1
r2=$2
ref=$3
outdir=$4
#MAIN
mkdir $outdir
#Adaptor clip raw reads
SeqPrep -L 15 -A AGATCGGAAGAGCAC  -B AGATCGGAAGAGCGT -f $r1 -r $r2 -1 $r1.clipped.gz -2 $r2.clipped.gz 2>> seq_prep_out
#Align reads
bowtie2 -x $ref -1 $r1.clipped.gz -2 $r2.clipped.gz -p 8 | samtools view -bS - > $outdir/$outdir.bam 2> $outdir/$outdir.alignment_stats  
#Deduplicate aligned reads
cd $outdir
samtools sort $outdir.bam $outdir.sorted
~mkircher/bin/pipeline2.0/FilterUniqueBAM.py $outdir.sorted.bam  --outprefix=$outdir.deduped
#Resort deduped stuff and generate BAM index
samtools sort $outdir.deduped.bam $outdir.deduped.sorted
samtools index $outdir.deduped.sorted.bam
python ~/python/git/vramani/generate_length_histogram.py /net/shendure/vol8/projects/HiC.TCC.DHC.project/nobackup/HiC_resources/mm10.chrom.sizes $outdir.deduped.sorted.bam False > ${outdir}.length 
