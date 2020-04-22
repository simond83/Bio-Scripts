# Bio-Scripts
Bioinformatics scripts used for data analysis in Mox study in Qiime 1

##Bash loop for panda seq for all files= 

for FILE in *_L001_R1_001.fastq
do FILE2=`echo $FILE | sed 's/_L001_R1_001.fastq/_L001_R2_001.fastq/g'`;
PRE=`echo $FILE | sed 's/_L001_R1_001.fastq//g'`;
pandaseq -f $FILE -r $FILE2 -u $PRE.unalined > $PRE.panda.fa;
done;

grep -c "^>" *.panda.fa;
grep -c "^>" *.unaligned;

for FILE in *.panda.fa;
do
grep -v "^>" $FILE | perl -ne 'print "'"$FILE"'". "\t" . length($_) . "\n";';
done >> length_all_samples.txt;

R

x<-read.table(file="length_all_samples.txt", header=F);
head(x)
names(x)<-c("Sample","Length")
boxplot(x$Length~x$Sample,ylim=c(0,500))

for FILE in *_L001_R1_001.fastq
do FILE2=`echo $FILE | sed 's/_L001_R1_001.fastq/_L001_R2_001.fastq/g'`;
PRE=`echo $FILE | sed 's/_L001_R1_001.fastq//g'`;
pandaseq -f $FILE -r $FILE2 -l 380 -L 420 -u $PRE.unalined > $PRE.panda.fa;
done;

mv *.panda.fa pandaseq_cutoff (this moves pandaseq output files)
mv *.unalined pandaseq_cutoff (again moves output files)

grep -c "^>" *.panda.fa;
grep -c "^>" *.unaligned;


for FILE in *.panda.fa;
do
grep -v "^>" $FILE | perl -ne 'print "'"$FILE"'". "\t" . length($_) . "\n";';
done >> length_all_samples.txt;

R

x<-read.table(file="length_all_samples.txt", header=F);
head(x)
names(x)<-c("Sample","Length")
boxplot(x$Length~x$Sample,ylim=c(0,500))

pico

add_qiime_labels.py -i pandaseq_alligned -m map3.txt -c InputFileName -n 1 -o combined_fasta

pick_otus.py -m usearch61 -i combined_seqs.fna -o usearch_qf_results --db_filepath=/usr/local/share/qiime/gg_otus-12_10-release/rep_set/97_otus.fasta;

pick_rep_set.py -m most_abundant -i combined_seqs_otus.txt -o rep_set.fna -f ../combined_seqs.fna

assign_taxonomy.py -m rdp -i rep_set.fna -o rdp_assigned_taxonomy --rdp_max_memory=100000;

align_seqs.py -i rep_set.fna -o pynast_aligned_seqs --template_fp=/usr/local/share/qiime/core_set_aligned.fasta.imputed --alignment_method pynast --pairwise_alignement_method uclust --min_percent_id 75.0 --min_length 150

make_otu_table.py -i combined_seqs_otus.txt -t rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o otu_table.biom -e pynast_aligned_seqs/rep_set_failures.fasta;

filter_samples_from_otu_table.py -i otu_table.biom -o otu_tabl_samples_to_keep.biom --sample_id_fp ../../samples_to_keep.txt

 summarize_taxa_through_plots.py -o taxa_summary -i otu_table.biom -m ../../map3.txt

firefox bar_charts.html

summarize_taxa_through_plots.py -o taxa_summary_groups -i otu_table.biom -m ../../map3.txt -c 

make_phylogony.py -i pynast_aligned_seqs/rep_set_aligned.fasta -o rep_set.tree -t fasttree -r tree_method_default;

alpha_refraction.py -i otu_table.biom -o arare_max100 -t rep_set.tree -m ../../map3.txt
multiple_rarefactions.py -i otu_table.biom -m 25000 -x 400000 -s 25000 -n 10 -o rare_2500-400000
alpha_diversity.py -i rare_25000-400000 -o alpha_rare_400000 -t rep_set.tree -m observed_species,chao1,PD_whole_tree

beta_diversity_through_plots.py -i otu_table.biom -o bdiv_even100 -t rep_set.tree -m ../../map3.txt

compare_alpha_diversity.py -i chao1.txt -o alpha_chao1_stats -m map.txt -t nonparametric -c Treatment

compare_alpha_diversity.py -i arareS/alpha_div_collated/chao1.txt -o alpha_chao1_stats -m map_t.txt -t nonparametric -c type_treatment

compare_categories.py --method adonis -i weighted_unifrac_dm.txt -o bdiv_stats_adonis_weighted -m map_treat.txt -c Treatment 


##Metabolomics in house matlab scripts## 

plot(x’)
a
b
x(:,a:b)=[];
ppm(:,a:b)=[];
a=doAlignment(x,ppm,median(x),0)
n=normalise(a,’total’)

PCA
mpca=pcamodel(x,0.5)
mpca.addClass('class name',find(y==0));
mpca.addClass('class name',find(y==1));
mpca.twoFirstComponents
mpca.addComponent
mpca.scorePlot(3,4,’classes’)
sx=std(x);
colorplot(ppm,mpca.P(:,1)’.*sx,mpca.P(:,1).^2)
mpca.twoFirstComponents
text(mpca.T(:,1),mpca.T(:,2),label)


OPLS-DA
mjro2pls=mjrMainO2pls(x,y,1,1,0,7,'nfold','mc','uv',[],'re','y','standard');
mjrO2plsSummaryPlot(mjro2pls,x,y,ppm);
[pv,Q2p, R2p,m] = JTPpermutate(x,y,1000,1,0,7,’uv’)










