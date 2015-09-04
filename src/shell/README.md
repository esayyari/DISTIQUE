Here are descriptions of some (maybe someday all) of the scripts belonging to this directory.
Each script would ideally give its usage (again someday!):

On MAC, you need to install the standard coreutils and then you need GNU command line tools for many of these. 
I found [this](https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/) link useful. 

* [setup.sh](setup.sh): this file is (should be) internally called by other scripts to setup some env variables globally used. 
  For it to work, you need to have `$WS_HOME` defined as explained [here](../../README).

* [annotate-with-genes.sh](annotate-with-genes.sh): reads a bunch of source trees and a destination tree, and annotates the branches of the destination tree 
  with the number (or frequency) of source trees that have that branch with high support

* [bootstrap_summary.sh](bootstrap_summary.sh): Summarizes the bootstrap support of branches. Outputs min, max, and average bootstrap support, plus the number of highly supported branches. Input is newick tree, and the script assumes bootstrap support values are given as node labels (so for example, `((a,b)76,(c,d)45)`).

* [contract_low_support.sh](contract_low_support.sh): contracts low support branches in a newick tree. 
  * this is just a wrapper for a perl script that I obtained from Kevin Liu many years ago. The perl script is packaged. Thanks Kevin! 

* [convert_to_phylip.sh](convert_to_phylip.sh): converts fasta alignments into phylip alignments. 

* [compareTrees.missingBranch](compareTrees.missingBranch): compare an estimated tree versus a reference tree, and outputs the 
  missing branch rate (FN). 
  * The output format is:

    `[number of branches in reference] [number of branches in the reference tree not found in the estiamte tree] [fraction of branches in refrence not found in estimated]`
 
  * This is just a wrapper around a perl script by Morgan Price, and available [here](http://www.microbesonline.org/fasttree/treecmp.html). Thanks Morgan!
  * No extra steps for installation should be required, except you need to have perl installed. 
  * The estimated tree can be multiple lines; in that case, the comparison is between the reference tree (only first line) and all estimated trees. 
  * use `-simplify` if the two trees are not on the same taxa but you want to force them to be
  * If you get a criptic error message, run [compareTrees](compareTrees) instead, and you should get more useful error messages

* [compareTrees](compareTrees): similar to [compareTrees.missingBranch](compareTrees.missingBranch), but instead of printing just the FN, it prints an extensive branch by branch comparison of the two trees. 
  Refer to Price's [website](http://www.microbesonline.org/fasttree/treecmp.html) for more info. 
  * This is useful for branch length and bootstrap support comparison. 

* [compareTrees.compatibility](compareTrees.compatibility): compare trees in terms of compatibility, and thus, the two trees need not be on the same exact set of taxa, 
 or fully binary. However, this code is just a wrapper around an old hacked [phylonet](http://bioinfo.cs.rice.edu/phylonet) code, with some bug fixes I added.
 This old phylonet code is packaged with the repository. Some limitation of the code are:
   * Internal node labels are problematic (if repeated). Usually, you need to remove internal node labels (e.g., bootstrap support). 
   * I think this old phylonet dislikes funky characters in taxon names; this includes quotation.

* [compareTrees.triplets](compareTrees.triplets): finds the triplet distance between two trees. 
  * This is just a wrapper around an excellent program from a [MS thesis by Jens Johansen](http://jensjohansen.com/thesis/). Thanks Jens!
  * The binaries of Jen's code are *not* included in this package. Please download the code from the link above, build it, and put the binaries under `$WS_HOME/bin` and name the binary file `triplets.soda2103`.
  * Note that the script can be easily adopted to do quartets instead of triplets (TODO). 

* [create_1stAnd2ndcodon_alignment.sh](create_1stAnd2ndcodon_alignment.sh): creates an alignment where the first half are the first codon positions of the input alignment and the second half is the second codon position. It also outputs a file that gives the boundaries between the two codons. 

* [distance.sh](distance.sh): calculates a bunch of statistics on a given fasta alignment file. 
  These include things like mean and max p-distance, gappiness, length, etc.
  * This is just a wrapper around a java code I obtained from Kevin Liu and packaged here. 

* [reformat_bootstrap_support.sh](reformat_bootstrap_support.sh): bootstrap support is sometimes shown in brackets (e.g., `(a,b)[98])`.
  To convert that to the format that most tools here understand (e.g., `(a,b)98`), you can use this script. 

* [remove_taxon_from_fasta.sh](remove_taxon_from_fasta.sh): removes a set of sequences from a fasta file.
  * The set of sequence names should be | delimited; so to remove human, gorilla, and chimp, the first argument should be `human|gorilla|chimp`.
  * You should be able to use more fancy regular expressions in the first argument (sequence name).
  * Use `-rev` as the third argument to keep sequences in the first argument and remove everything else

* [remove_empty_sequences_from_fasta.sh](remove_empty_sequences_from_fasta.sh): removes empty sequence headers from a fasta file

* [simplifyfasta.sh](simplifyfasta.sh): a one-liner that standardizes fasta alignment files so that each sequence is only one line. Simple and sweet with no dependencies. 

* [subsample.fasta.sh](subsample.fasta.sh): restricts a fasta file to a subset of columns. Does not remove all gap columns.

* [trim-fasta.sh](trim-fasta.sh): strips down a fasta file to its first N sites

* [triplet.freq.sh](triplet.freq.sh): given a file with a bunch of newick trees, computes the frequencies of the three topologies for all n choose 3 triplets and outputs results to standard output. 
  * This is just a wrapper around an excellent program from a [MS thesis by Jens Johansen](http://jensjohansen.com/thesis/). Thanks Jens!
  * The binaries of Jen's code are *not* included in this package. Please download the code from the link above, build it, and put the binaries under `$WS_HOME/bin` and name the binary file `trip_bin`.
  * Note that the script can be easily adopted to do quartets instead of triplets (TODO). 
  * python > 2.6 is a requirement (but python >3.0 would probably not work). 


