echo $1 #fasta
echo $2 #gtf
echo $3 #Mt /mt / MT
echo $4 #genome_name

mkdir -p  ./$4/Annotation/Genes/ ./$4/Sequence/WholeGenomeFasta
if [[ $1 == *.gz ]] && [[ $2 == *.gz ]]
then
   gunzip -c $1 > genome.fa
   gunzip -c $2 > genes.gtf
   mv ./genome.fa ./$4/Sequence/WholeGenomeFasta/
   mv ./genes.gtf ./$4/Annotation/Genes/
   grep "^$3" ./$4/Annotation/Genes/genes.gtf | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > mt_genes.txt
   grep rRNA ./$4/Annotation/Genes/genes.gtf | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > rRNA_genes.txt
else
   cp $1 ./$4/Sequence/WholeGenomeFasta/genome.fa
   cp $2 ./$4/Annotation/Genes/genes.gtf
   grep "^$3" $2 | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > mt_genes.txt
   grep rRNA $2 | grep -oP "gene_name \"\K[^\"]+" | sort | uniq > rRNA_genes.txt
fi
mv mt_genes.txt rRNA_genes.txt ./$4/Annotation/Genes/
outdir=./$4/Sequence/
output_fasta=./$4//Sequence/WholeGenomeFasta/genome.fa
output_gtf=./$4/Annotation/Genes/genes.gtf
mkdir ${outdir}/STARIndex
   
STAR --runMode genomeGenerate --genomeDir ${outdir}/STARIndex --genomeFastaFiles $output_fasta --sjdbGTFfile $output_gtf --sjdbOverhang 50 --runThreadN 10 --outTmpDir ${outdir}/STARtmp
mv  ./Log.out ./$4/Sequence/STARIndex/

