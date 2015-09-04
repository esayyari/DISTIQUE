Random scripts that do phylogenetic plumbing ,,,

* [summarize.triplets.py](summarize.triplets.py): this is called from [a shell script](../../shell/triplet.freq.sh) to summarize triplet frequencies across many trees

* [gc-stats.py](gc-stats.py): outputs the GC content statistics (and more) for a given alignment. 
  * The output is very elaborate, and has many fields. Here are the fields:
    * `"SEQUENCE","TAXON","A_C","C_C","G_C","T_C","A_R","C_R","G_R","T_R","c1_A_C","c1_C_C","c1_G_C","c1_T_C","c1_A_R","c1_C_R","c1_G_R","c1_T_R","c2_A_C","c2_C_C","c2_G_C","c2_T_C","c2_A_R","c2_C_R","c2_G_R","c2_T_R","c3_A_C","c3_C_C","c3_G_C","c3_T_C","c3_A_R","c3_C_R","c3_G_R","c3_T_R"`
    * The first element is the sequence name, the second field is the everything in the sequence name before `_`
    * Then there are counts of A, C, G, T, and then the ratio of A, C, G, and T. 
    * Then comes the counts and ratios for each of the three codon positions.
