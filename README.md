# plotRNAedits
plot RNA edits from vcf
This package utilizes MutationalPatterns to plot 192 tricnucleotide profile RNA edits from vcf files. 
The 192 profile is the full trinucleotide profile for nearest neighbour RNA edit context. 
The assumption is that the data come from stranded RNASeq experiments. This way the exact edit is presented. 
No condensing the reverse complements is performed as is the standard approach in somatic DNA mutation profiles.
