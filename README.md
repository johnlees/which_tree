Evaluation of phylogenetic reconstruction methods
=================================================
This repository contains code used to run the simulations and tree
comparisons in this paper:

[Evaluation of phylogenetic reconstruction methods using bacterial whole genomes: a simulation based study](https://wellcomeopenresearch.org/articles/3-33/v1)

Simulations
-----------
Simulations were run in steps as follows:

Using [ALF](http://alfsim.org/#index) for genes with the
parameters in `alf-tree-testing-sim.drw`. The required starting genome and tree
are included in this repository.
```
alfsim alf-tree-testing-sim.drw
```

Using [dawg](https://github.com/reedacartwright/dawg) for intergenic
regions, with parameter file `intergenic.dawg`:
```
dawg intergenic.dawg > intergenic_tree_test.mfa
```

These are then combined with the following commands:
```
./dawg_to_alf.pl speciesMapping.txt intergenic_coordinates.txt dawg_output.mfa
./alf_db_to_fasta_splice_intergenic.pl SE0xx_dna.fa dawg/MSA > SE0xx_genome.fa
```
Running the second line over all 96 sequences.

Finally, the resulting alignments (DB files) are converted into reads
and assembled:
```
pirs simulate -i DB/SE0xx_genome.fa -l 100 -x 50 -m 250 -v 80 -o reads/SE0xx
VelvetOptimiser.pl --s 37 --e 81 --x 4 --t 1 -f '-shortPaired -fastq.gz -separate reads/SE0xx_100_250_1.fq.gz reads/SE0xx_100_250_2.fq.gz' -p SE0xx
```

Reads are mapped with bwa and called with samtools to produce and
alignment:
```
bwa mem TIGR4_ref.fa SE0xx_1.fq.gz SE0xx_2.fq.gz | samtools fixmate -O bam - SE0xx.bam
samtools mpileup -C 50 -m 2 -F 0.0005 -d 1000 -t DP,SP -L 1000 -g -u -f TIGR4_ref.fa SE0xx.bam | bcftools call -vm -S samples.txt -O z -o SE0xx.vcf.gz
```

Associated data can be downloaded [here](https://doi.org/10.6084/m9.figshare.5483461).

Tree construction
-----------------
Software and versions are documented in the paper. Distance matrices are
converted to trees in the next step.

Command lines used, where appropriate:

andi:
```
andi -j -t 1 SE*.fa
```

bigsdb:
```
./bigsdb_distances.pl roary/pan_genome_sequences 96 > dist_mat.txt
```

raxml binary data:
```
raxmlHPC-PTHREADS-SSE3 -m BINGAMMAI -p 12345 -s gene_presence_absence.fa -n gene_presence_absence
```

raxml other data:
```
raxmlHPC-PTHREADS-SSE3 -T 4 -f d -s input.phy -m GTRGAMMA -p 12345 -n ml_mlst_alignment
```
with ascertainment bias correction:
```
raxmlHPC-PTHREADS-SSE3 -T 4 -f d -s input.phy -m ASC_GTRGAMMA -p 12345 -n ml_mlst_alignment --asc-corr=lewis
```

iqtree:
```
iqtree -s ../TIGR4_ref_bwa_snps.aln -nt 1 -pre gtr_g -m GTR+ASC+G
```

fasttree:
```
FastTree -nt -pseudo -fastest < TIGR4_ref_bwa.aln > tigr4_fasttree_fast.tree
FastTree -nt -gttr < TIGR4_ref_bwa.aln > tigr4_fasttree_slow.tree

mash:
```
mash sketch -o reference.msh *.fa
mash dist reference.msh reference.msh > mash_distances.txt
```

mrbayes:
```
mpirun mb < mrbayes_long.txt
```

ncd:
```
ls SE*.fa > assemblies.txt
./ncd_distance.pl assemblies.txt > ncd_distances.txt
```

parsnp:
```
parsnp -g Streptococcus_pneumoniae_TIGR4_v3.gbk -d ./assemblies -p 1 -c -x
```

progressiveCactus:
```
runProgressiveCactus.sh assemblies.txt ./ ./simulation_alignment.hal --database kyoto_tycoon --maxThreads 16
hal2maf ./simulation_alignment.hal ./simulation_alignment.maf
./maf2aln.pl --exclude Spn23F,TIGR4 -o simulation_alignment_fixed.aln
```

roary:
```
roary -p 4 -r *.gff
```

Comparing trees
---------------
All code is in `tree_compare.R`. The input trees can be downloaded from
[figshare](10.6084/m9.figshare.5483464), or load the `.rds` objects in the repo.

Comparison between the trees for real data is in `mass_trees.R`.

Interactive plots can be knitted from the `.rmd` files, or download the
[resulting html from figshare](https://doi.org/10.6084/m9.figshare.5923300).

