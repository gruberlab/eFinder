# E-Finder: a tool to find multigene elements in assembled sequences using profile HMMs.

e-Finder is a generic tool for detection and extraction of multigene elements from assembled genomes using profile HMMs. e-Finder executes hmmsearch program (HMMER package) to run similarity searches using profile as queries against translated sequences of the assembled genomes. Any region containing a cluster of genes can be detected, such as transposons, CRISPR-Cas systems, prophages, operons, etc.

# Instalation

E-Finder does not need to be installed. The user should only download the e-finder.pl file.

# Requirements

- transeq program: EMBOSS package - http://emboss.sourceforge.net/
- hmmsearch program: HMMER3 package - http://hmmer.org/
- tblastn program: BLAST package – http://blast.ncbi/nlm/nih/gov

# Usage
```
perl e-finder -df <file> -i <file> -s|-e <decimal>  <optional parameters>
```
### Mandatory parameters
```
-df|dataset_file  <file>        : Dataset (FASTA file)
-dd|dataset_dir <directory>	: Dataset directory.
-i|input_file                	: Input file (single or multiple profile HMMs).

```
### Optional parameters
```
-ce|circular<yes|no>		: Assume that the element is originally derived from a circular element (e.g. a prophage derived
				  from a circular phage genome). Default: no.
-conf <file>                    : Use a configuration file that lists all parameters for execution, overriding any parameter of
				  the command line.
-cpu|cpu_threads <integer>      : Number of threads to be used by hmmsearch. If not specified, e_Finder determines the number of
				  threads available in the multiprocessor server and uses by default half of this value.
-ed|element_distance <integer>	: Minimum distance between elements. Default = 5000.
-ex|extension <fna|faa|fasta>   : Extension for input files. When using PATRIC database data, we suggest using fna, since this is
				  the extension of the contig nucleotide sequence files provided in this resource. Default: fna.
-e|e-value | -s|score <decimal> : E-value (-e) or score (-s) threshold value. E-value (-e) or score (-s) threshold value. Report
				  hmmsearch hits that present values equal to or lower than the E-value or equal to or larger than
				  the score. Only one of these parameters and the respective value shall be provided. Default = -e 10.
-fs|flanking_size <integer>	: Size (in bp) of the 5’ and 3’ flanking regions that will be excised together with the multigene
				  element. If a 0 (zero) value is used, e-Finder will extract the sequence comprised between the start
				  codon of the first gene and the stop codon of the last gene. Default = 5000.
-gc|genetic_code  <integer>     : Genetic code to define start codons and perform conceptual translation of the genes. e-Finder uses
				  numbering codes defined on the NCBI Genetic Codes page and implemented on transeq (EMBOSS package):
				  0 (Standard); 1 (Standard with alternative initiation codons); 2 (Vertebrate Mitochondrial);
				  3 (Yeast Mitochondrial); 4 (Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma);
				  5 (Invertebrate Mitochondrial); 6 (Ciliate Macronuclear and Dasycladacean); 9 (Echinoderm Mitochondrial);
				  10 (Euplotid Nuclear); 11 (Bacterial); 12 (Alternative Yeast Nuclear); 13 (Ascidian Mitochondrial);
				  14 (Flatworm Mitochondrial); 15 (Blepharisma Macronuclear); 16 (Chlorophycean Mitochondrial);
				  21 (Trematode Mitochondrial); 22 (Scenedesmus obliquus); 23 (Thraustochytrium Mitochondrial).
				  Default: 0.
-h|help             		: Display help screen.
-ic|ignore_cutoff <yes|no>      : Ignore cutoff scores in the profile HMMs and use a custom value defined by parameters -e or -s
                                  for all input models. If -r no is used, e-finder will use the cutoff scores specified in the
				  respective CUTOFF SCORE tag of each profile HMM. For models not containing cutoff values, e-finder
				  will use the cutoff value specified by the parameter -e ou -s. If none of these parameters has been
				  specified, the program will then use hmmsearch's default cutoff value (-E 10). Default = yes.
-gc | ignore_genomes		: List of genomes previously analysed. These genomes will be ignored by this execution.
-id|intergenic_dist <integer>	: Maximum distance (in bp) between intergenic regions. This value is adopted for all intergenic distances
				  between any pair of genes. If the -synteny parameter is used, its values take precedence over those set
				  for -id. Default = 5000.
-mg|min_gene <integer>          : Minimum number of genes. This parameter specifies the minimum number of genes that must be found for a
				  sequence to be considered positive. For instance, if profile HMMs from four different proteins are used
				  and the user specifies -mg 2, any sequence containing at least two of the four markers will be considered
				  for downstream analysis of the remaining criteria. Default = 2.
-o|output <string>      	: This is the directory where e-Finder will store all output directories/file. All similarity search results
				  are stored in the all_results subdirectory. If the user specifies an output directory that already exists,
				  e-Finder inspects the all_results subdirectory and uses the hmmsearch results from the previous run. This
				  feature saves processing time since it skips the relatively slow similarity search step. For each run,
				  e-Finder creates a run_# (e.g. run_1, run_2, run_3, etc.) subdirectory where all output files are stored.
				  Default = output_dir.
-ol|overlap <integer>		: Maximum allowed overlap distance (in bp) between open reading frames in the same coding strand. If a 0 (zero)
				  value is used (default), no overlap is allowed.
-pl|patric_list	<string>	: Input file in PATRIC-like (https://www.patricbrc.org) tabular format. This two-column file lists accession codes
				  and organism names, respectively, and provides information for e-Finder to generate a final CSV file reporting all
				  found multigene regions associated with the respective organism names.
-sf|size_filter <integer> 	: Minimum size (bp) of the excised element, not including user-defined flanking regions (parameter -fs). Default: 1000.
-synteny <string>               : Define gene order and maximum allowed distance (kb) between genes. Each marker is defined by a letter and the distances
				  can be specified by decimals. Intergenic distances defined for parameter -synteny take precedence over those set for -id.
				  If parameter -synteny is not specified, e-Finder will accept any sequences presenting the minimum number of genes (specified
				  in parameter -mg), with any intergenic distances (up the maximum value defined by parameter -id).
				  Example:
					-Circular element: a,2000,b,1500,c,3000,d,3500,e,2500.
					-Linear element: a,2000,b,1500,c,3000,d,3500,e.
-v|version        		: display program’s version.
```

# User manual

# Tutorial

# Contact
To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Liliane S. Oliveira (liliane.sntn@gmail.com).
