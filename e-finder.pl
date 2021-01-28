#!/usr/bin/perl

# e-finder.pl v.1.0 (2021-01-25) - A tool to find multigene elements in assembled 
# sequences using profile HMMs

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use Scalar::Util qw(looks_like_number);

# variables

my $output = 'output_dir';
my $version = "1.0";
my $last_update = "2021-01-25";
my $help;
my $transeq_file;
my $fasta_qual;
my $cpu;
my $input_file;
my $input_dir;
my $file_log = "logfile.txt";
my $version_op;
my $hmm_db;
my $usr_score;
my $usr_evalue;
my $conf;
my %config;
my $ignore_cutoff = "yes";
my $flanking = 5000;
my $dist = 5000;
my $minimum_amount = 2;
my $extension_file = "fna";
my $synteny;
my $circular = "no";
my $patric_list;
my $overlap = 0;
my $genetic_code = 0;
my $table = 0;
my $size_filter = 1000;
my $element_distance = 5000;
my $script_command = $0;
my $ignore_genomes;


# Capture the commands given to e-Finder in the command line
foreach (@ARGV) {
    $script_command .= /\s/ ?   " \"" . $_ . "\""
                    :           " "   . $_;
}
 
# Define the genetic code to be used
my %number_codes = (
		    '0' => '0',
                    'standard' => '1',
                    '1' => '1',
                    'vertebrate' => '2',
                    'vertebrate mitochondrial' => '2',
		    'vertebrate mitochondrial code' => '2',
                    '2' => '2',
                    'yeast' => '3',
                    'yeast mitochondrial' => '3',
		    'yeast mitochondrial code' => '3',
                    '3' => '3',
                    'mold' => '4',
                    'protozoan' => '4',
                    'coelenterate mitochondria' => '4',
                    'mycoplasma' => '4',
		    'mycoplasma/spiroplasma' => '4',
                    '4' => '4',
                    'invertebrate' => '5',
                    'invertebrate mitochondrial' => '5',
		    'invertebrate mitochondrial code' => '5',
                    '5' => '5',
                    'ciliate' => '6',
                    'dasycladacean' => '6',
                    'hexamita nuclear' => '6',
                    'hexamita' => '6',
                    '6' => '6',
                    'echinoderm' => '9',
                    'echinoderm nuclear' => '9',
		    'echinoderm mitocondrial' => '9',
		    'flatworm' => '9',
		    'flatworm mitocondrial' => '9',
		    '9' => '9',
                    'euplotic nuclear' => '10',
                    'euplotic' => '10',
                    '10' => '10',
                    'bacterial' => '11',
		    'archaeal' => '11',
	     	    'archaeal code' => '11',
                    'bacterial code' => '11',
                    '11' => '11',
                    'alternative yeast nuclear' => '12',
                    '12' => '12',
                    'ascidian mitochondrial' => '13',
                    'ascidian' => '13',
                    '13' => '13',
                    'alternative flatworm mitochondrial' => '14',
                    'alternative flatworm' => '14',
                    '14' => '14',
                    'chlorophycean mitochondrial' => '16',
                    'chlorophycean' => '16',
                    '16' => '16',
                    'trematode mitochondrial' => '21',
                    'trematode' => '21',
                    '21' => '21',
                    'scenedemus obliquus mitochondrial' => '22',
                    'scenedemus' => '22',
                    '22' => '22',
                    'thraustochytrium' => '23',
                    'thraustochytrium mitochondrial' => '23',
                    '23' => '23');

# Define start codons of each genetic code
my %start_codon = (
		   '0' => 'ATG',
                   '1' => '(ATG|CTG|TTG)',
                   '2' => '(ATA|ATT|ATC|ATG|GTG)',
                   '3' => '(ATG|ATA)',
                   '4' => '(TTA|TTG|CTG|ATT|ATC|ATA|ATG|GTG)',
                   '5' => '(TTG|ATT|ATC|ATA|ATG|GTG)',
                   '6' => 'ATG',
                   '9' => '(ATG|GTG)',
                   '10' => 'ATG',
                   '11' => '(TTG|CTG|ATT|ATC|ATA|ATG|GTG)',
                   '12' => '(CTG|ATG)',
                   '13' => '(TTG|ATA|ATG|GTG)',
                   '14' => 'ATG',
                   '16' => 'ATG',
                   '21' => '(ATG|GTG)',
                   '22' => 'ATG',
                   '23' => '(ATT|ATG|GTG)');

# Defines possible stop codons for each genetic code
my %end_codon = (  '0' => '(TAA|TAG|TGA)',
                   '1' => '(TAA|TAG|TGA)',
                   '2' => '(TAA|TAG|AGA|AGG)',
                   '3' => '(TAA|TAG)',
                   '4' => '(TAA|TAG)',
                   '5' => '(TAA|TAG)',
                   '6' => 'TGA',
                   '9' => '(TAA|TAG)',
                   '10' => '(TAA|TAG)',
                   '11' => '(TAA|TAG|TGA)',
                   '12' => '(TAA|TAG|TGA)',
                   '13' => '(TAA|TAG)',
                   '14' => 'TAG',
                   '16' => '(TAA|TGA)',
                   '21' => '(TAA|TGA)',
                   '22' => '(TCA|TAA|TGA)',
                   '23' => '(TTA|TAA|TAG|TGA)');

my %letters = ('1' => 'a', 
	       '2' => 'b',
               '3' => 'c',
               '4' => 'd',
               '5' => 'e',
               '6' => 'f',
               '7' => 'g',
               '8' => 'h',
               '9' => 'i',
               '10' => 'j',
               '11' => 'k',
               '12' => 'l',
               '13' => 'm',
               '14' => 'n',
               '15' => 'o',
               '16' => 'p',
               '17' => 'q',
               '18' => 'r',
               '19' => 's',
               '20' => 't',
               '21' => 'u',
               '22' => 'v',		
               '23' => 'w',
               '24' => 'x',
               '25' => 'y',
               '26' => 'z'	
	      );

my $help_print = "##### e-Finder - version $version ($last_update) - L. Oliveira & A. Gruber #####

Usage:
e-finder.pl -df <file> -i <file> -s|-e <decimal>  <optional parameters>

Mandatory parameters:
-df|dataset_file  <file>        : Dataset (FASTA file)
-dd|dataset_dir <directory>	: Dataset directory.
-i|input_file                	: Input file (single or multiple profile HMMs).

OPTIONAL PARAMETERS:
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
 \n";

my $optret = GetOptions ("conf=s"               => \$conf,
			 "ce|circular=s"	=> \$circular,
			 "df|dataset_file=s"	=> \$input_file,
			 "dd|dataset_dir=s"	=> \$input_dir,
                         "i|input_file=s"       => \$hmm_db,
                         "s|score=f"            => \$usr_score,
                         "e|e-value=f"          => \$usr_evalue,
                         "o|output=s"       	=> \$output,
			 "ol|overlap"		=> \$overlap,
                         "h|help"               => \$help,
			 "cpu|cpu_threads=i"	=> \$cpu,
			 "ic|ignore_cutoff=s"  	=> \$ignore_cutoff,
			 "ig|ignore_genomes=s"	=> \$ignore_genomes,
			 "v|version"            => \$version_op,
			 "id|intergenic_dist=i"	=> \$dist,
			 "mg|min_gene=i"	=> \$minimum_amount,
			 "fs|flanking_size=i"	=> \$flanking,
			 "ex|extension=s"	=> \$extension_file,
			 "sf|size_filter=i"	=> \$size_filter,
                         "synteny=s"            => \$synteny,
			 "pl|patric_list=s"	=> \$patric_list,
			 "gc|genetic_code=s"	=> \$genetic_code,
			 "ed|element_distance=s"=> \$element_distance);

if($help){
    print $help_print;
    die "\n";    
}

if($version_op){
    die "Version $version.\nLast update: $last_update\n";
}

# If configuration file is not specified, check if the mandatory arguments are defined 
if (not defined $conf) {
    if(!$input_file and !$input_dir){
    	die "ERROR: Missing mandatory argument -df|dataset_file or -dd|dataset_dir.\n$help_print\n";
    }
    
    if(!$input_file and !$input_dir){
	die "ERROR: Please use -df|dataset_file or -dd|dataset_dir, not both.\n";
    }

    if(!$hmm_db){
    	die "ERROR: Missing mandatory argument -i.\n$help_print\n";
    }

    if(($usr_score) and ($usr_evalue)) {
    	die "ERROR: Please use -s or -e, not both.\n$help_print\n";
    }
}

# If configuration file is specified, check if all mandatory arguments are defined
if ($conf) {
    print STDERR "Configuration file specified, command line options will be overridden.\n";

    open(CONFIG, "< $conf") or die("ERROR: Problem opening configuration file $conf: $!\n");

    my $configLine;
    while ($configLine = <CONFIG>) {
        $configLine =~ s/^\s+//;
        $configLine =~ s/\s+\Z//;
        $configLine =~ s/\s+\=/\=/;
        $configLine =~ s/\=\s+/\=/;
        if ($configLine =~ /^\#/ || !($configLine =~ /(.)\=(.)/)) {
            next;
        }
        chomp $configLine;

        if ($configLine =~ m/(.+?)=(.+)/) {
            $config{$1} = $2;
        }
    }
    close(CONFIG);

    # Mandatory arguments
    my $missingArgument = 0;

    if (!$config{"i"} and !$config{"input_file"}) {
        $missingArgument = 1;
        print "Missing mandatory configuration argument or file not found: input file.\n";
    }
    else {
	if($config{"i"}){
            $hmm_db = $config{"i"};
	}
	elsif($config{"input_file"}){
            $hmm_db = $config{"input_file"};
        }
    }
    if (!($config{"df"}) and !($config{"dataset_file"}) and !($config{"dd"}) and !($config{"dataset_dir"})) {
        $missingArgument = 1;
        print "Missing mandatory configuration argument: dataset.\n";
    }
    elsif((($config{"df"}) or ($config{"dataset_file"})) and (($config{"dd"}) or ($config{"dataset_dir"}))){
	$missingArgument = 1;
	print "ERROR: Please use -df|dataset_file or -dd|dataset_dir, not both.\n";
    }
    else {
	if($config{"df"}){
            $input_file = $config{"df"};
	}
	elsif($config{"dataset_file"}){
            $input_file = $config{"dataset_file"};
        }
	elsif($config{"dd"}){
            $input_dir = $config{"dd"};
        }
	elsif($config{"dataset_dir"}){
            $input_dir = $config{"dataset_dir"};
        }
	
    }

        if ($missingArgument) {
        die "\nERROR: Cannot run e-finder.pl, mandatory configuration argument(s) missing (see bellow).\n\n$help_print\n";
    }

    # Optional parameters
    if (((defined($config{"e"}) or defined($config{"e-value"}))) and ((defined($config{"s"}) or (defined($config{"score"}))))){
	die "ERROR: Please use -s|score or -e|e-value, not both.\n$help_print\n";	
    }

    if (defined($config{"ig"})){
        $ignore_genomes = $config{"ig"};
    }
    elsif(defined($config{"ignore_genomes"})){
        $ignore_genomes = $config{"ignore_genomes"};
    }


    if (defined($config{"circular"})){
        $circular = $config{"circular"};
    }
    elsif(defined($config{"ce"})){
        $circular = $config{"ce"};
    }

    if (defined($config{"e"})){
	$usr_evalue = $config{"e"};
    }
    elsif(defined($config{"e-value"})){
        $usr_evalue = $config{"e-value"};
    }
 
    if (defined($config{"s"})){
        $usr_score = $config{"s"};
    }
    elsif(defined($config{"score"})){
        $usr_score = $config{"score"};
    }
    
    if (defined($config{"o"})){
        $output = $config{"o"};
    }
    elsif(defined($config{"output"})){
        $output = $config{"output"};
    }

    if (defined($config{"cpu_threads"})){
        $cpu = $config{"cpu_threads"};
    }
    elsif (defined($config{"cpu"})){
        $cpu = $config{"cpu"};
    }

    if (defined($config{"ic"})){
        $ignore_cutoff = $config{"ic"};
    }
    elsif(defined($config{"ignore_cutoff"})){
        $ignore_cutoff = $config{"ignore_cutoff"};
    }

    if (defined($config{"id"})){
        $dist = $config{"id"};
    }
    elsif(defined($config{"intergenic_dist"})){
        $dist = $config{"intergenic_dist"};
    }
    
    if (defined($config{"fs"})){
        $flanking = $config{"fs"};
    }
    elsif(defined($config{"flanking_size"})){
        $flanking = $config{"flanking_size"};
    }
    
    if (defined($config{"mg"})){
        $minimum_amount = $config{"mg"};
    }
    elsif(defined($config{"min_gene"})){
        $minimum_amount = $config{"min_gene"};
    }
    
    if (defined($config{"ex"})){
        $extension_file = $config{"ex"};
    }
    elsif(defined($config{"extension"})){
        $extension_file = $config{"extension"};
    }
    if(defined($config{"synteny"})){
        $synteny = $config{"synteny"};
    }
    if(defined($config{"pl"})){
        $patric_list = $config{"pl"};
    }
    elsif(defined($config{"patric_list"})){
        $patric_list = $config{"patric_list"};
    }
    if(defined($config{"ol"})){
        $overlap = $config{"ol"};
    }
    elsif(defined($config{"overlap"})){
        $overlap = $config{"overlap"};
    }
    if(defined($config{"gc"})){
        $genetic_code = $config{"gc"};
    }
    elsif(defined($config{"genetic_code"})){
        $genetic_code = $config{"genetic_code"};
    }
    if(defined($config{"sf"})){
        $size_filter = $config{"sf"};
    }
    elsif(defined($config{"size_filter"})){
        $size_filter = $config{"size_filter"};
    }
    if(defined($config{"ed"})){
        $element_distance = $config{"ed"};
    }
    elsif(defined($config{"element_distance"})){
        $element_distance = $config{"element_distance"};
    }
}

if($genetic_code > 23 or $genetic_code < 0){
    die "Invalid genetic code $genetic_code!\n";
}
else{
    $table = $number_codes{$genetic_code};
}

if(!defined $synteny and (lc($circular) eq "yes")){
    print "Parameter -circular is \"yes\" but -synteny is not defined. The value of -circular will be \"no\"\n";
}

if(!$cpu){
    $cpu = `grep CPU -c /proc/cpuinfo`;
    $cpu = $cpu/2;
}

if((lc($circular) eq "yes") and (defined $synteny)){
    if($synteny =~ /^(\+|\-)?(\w+)(\,)(\d+)((\,)(\+|\-)?(\w+)(\,)(\d+))*$/g){}
    else{
	die "Synteny is in a wrong format to circular genome! Example: a,2,b,1.5,c,3,d,4,e,5.\n";
    }    
}

if((lc($circular) eq "no") and (defined $synteny)){
    if($synteny =~ /^(\+|\-)?(\w+)(\,)(\d+)((\,)(\+|\-)?(\w+)(\,)(\d+))*(\,)(\+|\-)?(\w+)$/g){}
    else{
        die "Synteny is in a wrong format to non-circular genome! Example: a,2,b,1.5,c,3,d,4,e\n";
    }
}

my $log = "error.log";
open(my $log_file_handle, ">$log") or die "ERROR: Could not create log file $log : $!\n";

my @hmms_db = split(",", $hmm_db);

my $total_of_models = scalar(@hmms_db);

my $transeq_name;
my %hmms_score = ();

my @aux_print = split(" ",$script_command);

if($element_distance < $dist){
    $element_distance = $dist;
}

my %genomes_not = ();
if(defined $ignore_genomes){
    open(FILE, $ignore_genomes);
    while(<FILE>){
        chomp($_);
        $genomes_not{$_} = 1;
    }
    close(FILE);
}

my $selected_dir;
my $discarded_dir;

# If the input is a FASTA file 
if($input_file){
    my $len = length($input_file);
    my $dir = undef;
    my $pos = rindex($input_file, "/");
    if($pos != -1){
    	$dir = substr $input_file, 0, $pos;
    }
    my $aux_name = substr $input_file, (rindex($input_file, "/") + 1), $len;
    my $file_name = $aux_name;
    my $aux;
    my $resp = 0;
    my $hmmsearch_tabular_file;
    # Verify if the input file exists
    if(!-e $input_file){
	die "Database file $input_file not found!\n";
    }
    else{
	if(!-e $output){
    	    system "mkdir $output";
	}
	# Create the output directory
	my $all_results_dir = $output."/all_results";
	if(!-e $all_results_dir){
	    system "mkdir $all_results_dir";
	}
        my $run_dir = $output."/run_1";
	$run_dir = output_dir_name($run_dir);		
	system "mkdir $run_dir";	
	$file_log = $run_dir."/".$file_log;

	# Print the execution parameters in the logfile
	open(LOG, ">$file_log") or die "ERROR: Could not create log file $file_log : $!\n";
	print LOG "Command line:\n$script_command\n\n";
	print LOG "Parameters: \n\n";
	print LOG "-df dataset file: $input_file\n";
 	print  LOG "-i input file: $hmm_db\n";
	if(defined $extension_file){
	     print LOG "-ex extension file: -ex does not apply to dataset file\n"; 
	}
	if($synteny){
	    print LOG "-synteny: $synteny\n";
	    print LOG "-circular: $circular\n";
	}
	print LOG "-mg minimum genes: $minimum_amount\n";
	print LOG "-fs flanking size: $flanking\n";
	print LOG "-id intergenic_dist: $dist\n";
	print LOG "-ic ignore_cutoff: $ignore_cutoff\n";
	print LOG "-ol overlap: $overlap\n";
	print LOG "-gc genetic_code: $genetic_code\n";
	print LOG "-sf size filter: $size_filter\n";
        print LOG "-ed element_distance: $element_distance\n";
	my %pls = ();

	# Check if a file containing the genome description is provided
	if(defined $patric_list){
	    print LOG "-pl patric list: $patric_list";
	    if(defined $patric_list){
    	    	if(!-e $patric_list){
        	    print "Warning:File $patric_list not found. The final report will be generate without Organism name\n";
    		}
    		else{
		    # Store genome information in memory
        	    open(FILE, $patric_list); 
        	    while(<FILE>){
            	    	chomp($_);
            	    	if($_ =~ /genome_id	genome_name	taxon_id	genome_length/){			
			    print "Warning:File $patric_list in a wrong format for dataset file. The final report will be generate without Organism name\n";
			    print LOG "Warning:File $patric_list in a wrong format for dataset file. The final report will be generate without Organism name\n";
			    $patric_list = undef;
			    last;
			}
			if($_ =~ /contig_id/){}
            		else{
               		     my @aux = split("\t", $_);
                	     $pls{$aux[0]} = $aux[1];
           		}
          	     }
        	     close(FILE);
    		}
	    }
	}
	else{
	    print LOG "-pl patric list: Patric-like list file not provided\n";
	}
	print LOG "\nDataset analysis:\n";

	# Check if the provided file is a FASTA file
	my $ext = verify_file_type($input_file);
	if($ext != 1){
	    close(LOG);
	    system "rm -rf $output";
	    die "File $input_file is not a fasta file! Aborting...\n";	
	}
	if($ext == 1){
	    
	    # Check the content of the provided FASTA file (1 - nucleotide; 2 - amino acid)
            my $type = verifiesFastaFile($input_file);
	    if($type == 1){ # DNA fasta file - verify if transeq file exists
	    	print LOG "Type: Nucleotide Fasta\n";
		my $msg;
		
		# Run transeq program
                ($transeq_name, $msg) = runTranseq($file_name, $input_file, $dir);
                print LOG $msg;
     	    }
    	    else{# Protein FASTA file
	    	print LOG "Type: Protein Fasta\n";
	    	$transeq_name = $input_file;	    
    	    }	

	    # Store the genomes in memory
	    open(FILE, $input_file);
	    my $seq = "";
	    my $name = undef;
	    my %genomes = ();
	    while(<FILE>){
	    	chomp($_);
	    	if($_ =~ /^>/){
		    if(defined $name){
		   	$genomes{$name} = $seq;
		   	$seq = ""; 
		    }
		    my @aux = split(" ", $_);
		    $name = $aux[0];
		    $name =~ s/>//g;	
	    	}
	    	else{
		    $seq .= $_;
	    	}
	    }
	    close(FILE);
	    $genomes{$name} = $seq;
	    my $prefix_dir;

	    # Perform similarity searches for each provided profile HMM
	    foreach my $hmm (@hmms_db){
	    	%hmms_score = ();
	    	my $len = length($hmm);
	    	my $aux_name = substr $hmm, (rindex($hmm, "/") + 1), $len;
	    	$prefix_dir = substr $aux_name, 0, rindex($aux_name, ".");
	    	my $dir = $all_results_dir."/".$prefix_dir;
	    	if(!-e $dir){
	    	     system "mkdir $dir";
	    	}
		my $msg;

		# Check if a cutoff score must be used in similarity searches
	    	($resp, $msg) = verifyScoreEvalue($hmm, $dir);
		print LOG $msg;
	    	my $aux = $dir."/". $file_name;		
		if(!-e $transeq_name){
		    next;
		} 
		
		# Run hmmsearch program
	    	($hmmsearch_tabular_file, $msg) = runHmmsearch($aux, $transeq_name, $resp, $hmm, $dir);
		print LOG $msg;
		if($hmmsearch_tabular_file eq "-1"){
		    next;
		}

		# Save similarity searches results in a file
	    	print "Analysing hmmsearch results...\n";
                my $str_table1 = analyseHmmsearchResults($hmmsearch_tabular_file, $resp);
	    	if(!$str_table1 eq ""){
    		    my $table_1 = "$dir/table1.csv";
    		    open(TBL1, ">$table_1") or die "ERROR: Could not create file $table_1!\n";
    		    print TBL1 "Target Name\tquery_pHMM\tE-value\tScore\n";
    		    print TBL1 $str_table1;
    		    close(TBL1);
	    	}
		my $aux_hmms_dir = $dir."/hmms";
		if(-e $aux_hmms_dir){
		    system "rm -rf $aux_hmms_dir";
		}
	    }

	    # Analysis of the similarity searches results obtained by all profile HMMs
	    print LOG "\nResults: \n\n";
    	    opendir(DIR, "$all_results_dir");
	    my $abs = abs_path($all_results_dir);
    	    my @dirs = readdir(DIR);
    	    closedir(DIR);
    	    my %models = %{selectModels(\@dirs, $abs)};
	    my %content = ();
	    my %contigs = ();
	    foreach my $key (sort keys %models){
	    	my %hash = %{$models{$key}};
	        foreach my $key2 (sort keys %hash){
	    	    push @{$content{$key2}}, $hash{$key2};
		    push @{$contigs{$key2}}, $key;
	    	}
	    }
	    $selected_dir = $run_dir."/selected";	
	    system "mkdir $selected_dir";
	    my $partial_report = $selected_dir."/partial_report.csv";
	    open(PARTIAL, ">$partial_report"); 	 	
	    $discarded_dir = $run_dir."/discarded";
            system "mkdir $discarded_dir";
	    my $final_report = $run_dir."/final_report.csv";
	    open(FINAL, ">$final_report");
	    my $final_fasta = $selected_dir."/elements.fasta";
            open(FASTA, ">$final_fasta");
	    close(FASTA);
            my $str_report = "";
	    my %markers;
	    my %general_markers = ();

	    # Check if the results found for each input sequence meet the criteria chosen in
	    # the user-defined parameters 
	    
	    foreach my $key (sort keys %contigs){
	    	my @coord = ();
	    	my $amount = scalar(@{$contigs{$key}});
	    	
		# Check if the number of genes is equal to or larger than the minimum number 
		# defined by the user 
		
	    	if($amount >= $minimum_amount){
		    my $dir = $selected_dir."/".$key;
		    my $disc = $discarded_dir."/".$key;
		    system "mkdir $dir";
		    my $markers_str = "";
		    my @possible_proteins = ();
		    my $result_name = $hmmsearch_tabular_file;
		    my $len = length($result_name);
		    my $aux_name = substr $result_name, (rindex($result_name, "/") + 1), $len;
                    my $prefix = substr $aux_name, 0, rindex($aux_name, ".");
		    my @contend_marker = @{$content{$key}};
		    my $temp = $dir."/blast_temp_".$key.".fasta";
                    open(FILE, ">$temp");
		    print FILE ">$key\n$genomes{$key}\n";
		    close(FILE); 
		    %markers = ();

		    # Map the genes found in the original sequence to obtain their coordinates
		    for (my $i = 0; $i < scalar(@contend_marker); ++$i){
		    	my @aux_marker = split("_", $contigs{$key}[$i]);
		    	my $marker = $aux_marker[0];
		    	my @aux = @{$contend_marker[$i]};
			my %real_orf = ();
			@aux = @{sortByScoreMax(\@aux)};
			for(my $j = 0; $j < scalar(@aux); ++$j){
                       	    my @a = @{$aux[$j]};
                            my $file = $a[3]."/".$prefix.".txt";
                            my $subject = $a[2];
			    my $fasta = $dir."/".$key."_".$marker.".fasta";
                            generateFastaFileByHmmsearch($file ,$a[1],  $marker, $subject, $fasta);
                            if(!-e $fasta){
                            	next;
                            }
			    my $blast = $dir."/".$key."_".$marker."_blast.txt";
                            runBlast($fasta, $temp, $blast);
			    if(!-e $blast){	
				if(-e $fasta){
                                    system "rm $fasta";
                            	}
				system "rm $blast";
				next;
			    }
                            my @coord_fasta = @{extractCoordinates($blast, $key)};
                            if(scalar(@coord_fasta) > 0){
                            	for(my $index = 0; $index < scalar(@coord_fasta); ++$index){
                                    my @aux = split("\t", $coord_fasta[$index]);
                                    my $real = findORF($aux[0], $aux[1], $aux[2], $genomes{$key});
				    if(defined $real){
                                    	$real .= "\t".$marker;
                                    	$real_orf{$real} = 1;
				    }
                                }
                            }
                            else{
                                next;
                            }
			    if(-e $blast){
			    	system "rm $blast";
			    }
			    if(-e $fasta){
                                system "rm $fasta";
                            }
                            if(!defined $markers{$marker}){
                                $markers{$marker} = 1;
                            }
			    if(!defined $general_markers{$marker}){
                                $general_markers{$marker} = 1;
                            }
                        }
                        foreach my $key (sort keys %real_orf){
                            push @possible_proteins, $key;
                        }
		    } 
		    system "rm $temp";
		    @possible_proteins = @{sortCoordinates(\@possible_proteins)};
                    my @groups = ();
                    my $org_name = undef;
                    if(defined $patric_list){
                    	$org_name = $pls{$key};
                    }
		    my $index = 1;
		    my $prev_element = -1;

		    # Check if the gene order and orientation meet the criteria defined in the 
		    # parameters
		    for(my $i = 0; $i < scalar(@possible_proteins) - 1; ++$i){
                    	my @aux = split("\t", $possible_proteins[$i]);
                        my @aux2 = split("\t", $possible_proteins[$i+1]);
                        my $distance = $aux2[0] - $aux[1];
                        push @groups, $possible_proteins[$i];
                        if($synteny){
			    my @dists = split(",", $synteny);
                            for(my $j = 1; $j < scalar(@dists); $j = $j + 2){
                                if($distance > $dists[$j]){
				    @groups = @{discard_redundant_genes(\@groups)};
                                    my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir, 1, $org_name, $index);
                                    if($resp == 1){
					@groups = ();
                                   	next;
                                    }
				    @groups = ();
                                    $index = $id;
                                    ++$index;
                                    $str_report .= $str;
                                 }
                                 elsif(($distance <= $dist) and ($i == scalar(@possible_proteins) - 2)){
                                     push @groups, $possible_proteins[$i+1];
				     @groups = @{discard_redundant_genes(\@groups)};
                                     my $size_g = scalar(@groups);
                                     my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir, 1, $org_name, $index);
                                     if($resp == 1){
					@groups = ();
                                     	next;
                                     }
				     @groups = ();
                                     $index = $id;
                                     ++$index;
                                     $str_report .= $str;
                                 }
                                 elsif(($distance > $dist) and (scalar(@groups) < $minimum_amount)){
                                     @groups = ();
                                 }
                                 elsif(($distance < 0) and (abs($distance) > $overlap)){
                                     ++$i;
                                 }
                             }
			}
                        else{
			    if(($distance > $dist) and (scalar(@groups) >= $minimum_amount)){
				@groups = @{discard_redundant_genes(\@groups)};
                            	my $size_g = scalar(@groups);
                                my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir, 1, $org_name, $index);
                                if($resp == 1){
				    @groups = ();
                                    next;
                                }
				@groups = ();
                                $index = $id;
                                ++$index;
                                $str_report .= $str;
                            }
                            elsif(($distance <= $dist) and ($i == scalar(@possible_proteins) - 2)){
                            	push @groups, $possible_proteins[$i+1];
				@groups = @{discard_redundant_genes(\@groups)};
                                my $size_g = scalar(@groups);
                                my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir, 1, $org_name, $index);
                                if($resp == 1){
				    @groups = ();
                                    next;
                                }
				@groups = ();
                                $index = $id;
                                ++$index;
                                $str_report .= $str;
                            }
                            elsif(($distance > $dist) and (scalar(@groups) < $minimum_amount)){
                                @groups = ();
                            }
                            elsif(($distance < 0) and (abs($distance) > $overlap)){
                                ++$i;
                            }
                        }
                    }
		    if(scalar(@groups) > 0){
                    	my $i = scalar(@possible_proteins) - 1;
                        push @groups, $possible_proteins[$i];
			@groups = @{discard_redundant_genes(\@groups)};
                        my ($resp, $str, $current) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir, 1, $org_name, $index);
                        if($resp == 1){
                            next;
                        }
                        $str_report .= $str;
                    }		    
	    	}
		else{ # The minimum number of genes was not found in the sequence. The results 
		# will be stored in a specific directory.
		    my @contend_marker = @{$content{$key}};
		    my $markers_str = "";
		    for (my $i = 0; $i < scalar(@contend_marker); ++$i){
                        my @aux_marker = split("_", $contigs{$key}[$i]);
                        my $marker = $aux_marker[0];
                        if($i == scalar(@contend_marker)-2){
                            $markers_str .= $marker." and ";
                        }
                        elsif($i == scalar(@contend_marker)-1){
                            $markers_str .= $marker;
                        }
                        else{
                            $markers_str .= $marker.", ";
                        }                    
		    }
		    if($amount == 0){
                    	print LOG "$dir discarded - it has not markers in any contig.\n";
                    }
                    else{
                    	print LOG "Sequence $key positive to $markers_str discarded - number of markers is lower than $minimum_amount\n";
                        my $disc = $discarded_dir."/".$key;
                        if(!-e $disc){
                            system "mkdir $disc";
                        }
                        my $temp = $disc."/".$key.".fasta";
                        open(FILE, ">$temp");
                        print FILE ">$key\n$genomes{$key}\n";
                        close(FILE);
                        my $result_name = $hmmsearch_tabular_file;
                        my $len = length($result_name);
                        my $aux_name = substr $result_name, (rindex($result_name, "/") + 1), $len;
                        my $prefix = substr $aux_name, 0, rindex($aux_name, ".");
                        my @contend_marker2 = @{$content{$key}};
			for (my $i = 0; $i < scalar(@contend_marker2); ++$i){
                            my @aux_marker = split("_", $contigs{$key}[$i]);
                            my $marker = $aux_marker[0];
                            my @aux = @{$contend_marker[$i]};
			    for(my $j = 0; $j < scalar(@aux); ++$j){
                            	my @a = @{$aux[$j]};
                            	my $file = $a[3]."/".$prefix.".txt";
                            	my $subject = $a[2];
                            	my $fasta = $disc."/".$key."_".$marker.".fasta";
                            	generateFastaFileByHmmsearchDiscarded($file ,$a[1],  $marker, $subject, $fasta);
                            	if(!-e $fasta){
                            	    next;
                            	}
                            	my $blast = $disc."/".$key."_".$marker."_blast.txt";
                            	runBlast($fasta, $temp, $blast);
			        if(!-e $blast){
				    next;
				}
                            	my @coord_fasta = @{extractCoordinates($blast, $key)};
                            	if(scalar(@coord_fasta) > 0){
                            	    for(my $index = 0; $index < scalar(@coord_fasta); ++$index){
                                    	my @aux = split("\t", $coord_fasta[$index]);					
                                    	my $real = findORF($aux[0], $aux[1], $aux[2], $genomes{$key});
					if(defined $real){
                                    	    $real .= "\t".$marker;
                                    	    my @coord;
                                    	    push @coord, $real;
                                    	    generateProteinFasta2(\@coord, $genomes{$key}, $disc, $key);
					}
                                    }	
                            	}
                            	else{
                                    next;
                            	}
				if(-e $blast){
                                    system "rm $blast";
                                }
                                if(-e $fasta){
                                    system "rm $fasta";
                                }
                            }
                   	}
		    }
		}
	    }

	    # Generate the final files containing the positive sequences
	    if(defined $patric_list){
		print FINAL "ContigID_elem#\tOrganism name\tContig size\tElement coord. on contig\tElement size\t";
	    }
            else{
	    	print FINAL "ContigID_elem#\tContig size\tElement coord. on contig\tElement size\t";
	    }
	    my $lenght = keys %general_markers;
	    my $count = 1;	
	    foreach my $key (sort keys %general_markers){
	    	my $protein_file = $selected_dir."/".$key."_proteins.fasta";
		my $aux_file = $selected_dir."/saida";
            	system "ls $selected_dir/*/*$key*_protein.fasta >> $aux_file 2> $aux_file";
            	open(AUX, $aux_file);
            	my $num = 1;
            	while(<AUX>){
            	    chomp($_);
                    if($_ =~ /ls:/){
                    	$num = 0;
                    	last;
                    }
            	}
            	close(AUX);
           	system "rm $aux_file";

            	if($num == 0){		    
                    next;
            	}
	     	my $command = "cat $selected_dir/*/*$key*_protein.fasta > $protein_file";
	     	system "$command";	     
	     	if($count < $lenght){
		    print FINAL "Gene\tContig coord.\tElement coord.\tOrientation\tDistance(bases)\t";
	     	}
	     	else{
		    print FINAL "Gene\tContig coord.\tElement coord.\tOrientation\t";
	     	}


	    my %sequences = ();
            open(FILE, $protein_file);
            my $name = undef;
            my $seq = "";
            while(<FILE>){
                chomp($_);
                if($_ =~ /^>/){
                    if($seq ne ""){                    
			$seq =~ s/\*//g;
                        push @{$sequences{$name}}, $seq;
                    }
                    $name = $_;
                    $seq = "";
                }
                else{
                    $seq .= $_;
                }
            }
	    $seq =~ s/\*//g;
            push @{$sequences{$name}}, $seq;
            close(FILE);
            my $aux_fasta = $selected_dir."/$key\_proteins_temp.fasta";
            open(FILE, ">$aux_fasta");
            foreach my $key (sort {$a cmp $b} keys %sequences){
		my @aux_seq = @{$sequences{$key}};
		for(my $a = 0; $a < scalar(@aux_seq); ++$a){
                    print FILE "$key\n$aux_seq[$a]\n";
		}
            }
            close(FILE);
            system "mv $aux_fasta $protein_file >> /dev/null";

	    }
	    print FINAL "\n$str_report";
	    close(FINAL);
	    close(PARTIAL);
	    opendir(DIR, "$selected_dir");
            my @final_files = readdir(DIR);
            closedir(DIR);
	    if(-e $partial_report){
		system "rm $partial_report";
	    }
            foreach my $fasta_file (@final_files){
            	if($fasta_file eq "." or $fasta_file eq ".."){}
            	else{
                    my $aux = $selected_dir."/".$fasta_file;
                    if(!-d $aux){}
                    elsif(-e $aux){
			my $aux_file = $selected_dir."/saida";
			system "ls $aux/*.fna >> $aux_file 2> $aux_file";
			open(AUX, $aux_file);
			my $num = 1;
			while(<AUX>){
			    chomp($_);
			    if($_ =~ /ls:/){
				$num = 0;
				last;
			    }
			}
			close(AUX);
			system "rm $aux_file";
			if($num > 0){
                    	    my $command = "cat $aux/*.fna >> $final_fasta";
                    	    system $command;
			}
			else{
			    if(-e "$discarded_dir/$fasta_file"){
				system "mv $aux/* $discarded_dir 2> /dev/null";		
				system "rm -rf $aux 2> /dev/null";
			    }
			    else{
			    	system "mv $aux $discarded_dir >> /dev/null";
			    }
			}
                    }
            	}
           }
	    close(LOG);	
    	}
    }
}
# If the input is a directory
elsif($input_dir){
    if(!-e $input_dir){ # Check if the provided directory exists
        die "Database directory $input_file not found!\n";
    }    
    elsif(!-d $input_dir){
        die "$input_dir is not a directory!\n";
    } 
    else{
	my $caracter = substr $output, -1;
        if($caracter eq '/'){
            my $l = length($output);
            $output = substr $output, 0, ($l-1);
        }

	# Create the output directory
	if(!-e $output or (-e $output and !-d $output)){
	    system "mkdir $output";
	}
	my $run_dir = $output."/run_1";
	$run_dir = output_dir_name($run_dir);
	system "mkdir $run_dir";
	$file_log = $run_dir."/".$file_log;
	my $all_results_dir = $output."/all_results";
        if(!-e $all_results_dir){
	    system "mkdir $all_results_dir";
	}

	# Print the execution parameters in logfile
        open(my $fl, ">$file_log") or die "ERROR: Could not create log file $file_log : $!\n";
        print $fl "Command line:\n$script_command\n\n";
	print $fl "Parameters: \n\n";
        print $fl "-dd dataset dir: $input_dir\n";
        print $fl "-i input file: $hmm_db\n";
	print $fl "-ex extension file: $extension_file\n";
        if($synteny){	    
            print $fl "-synteny: $synteny\n";
	    print $fl "-circular: $circular\n";
        }
        print $fl "-mg minimum genes: $minimum_amount\n";
        print $fl "-fs flanking size: $flanking\n";
        print $fl "-id intergenic_dist: $dist\n";
        print $fl "-ic ignore_cutoff: $ignore_cutoff\n-gc genetic code: $genetic_code\n";
	print $fl "-ol overlap: $overlap\n";
	print $fl "-sf size filter: $size_filter\n";
	print $fl "-ed element_distance: $element_distance\n";
	my %pls = ();

	# Check if a file containing the genome description is provided
	if(defined $patric_list){
            print $fl "-pl patric list: $patric_list";
    	    if(!-e $patric_list){
                print "Warning:File $patric_list not found. The final report will be generate without Organism name\n";
    	    }
    	    else{
            	open(FILE, $patric_list);

		# Store genome information in memory
        	while(<FILE>){
            	    chomp($_);
		    if($_ =~ /contig_id/){
		    	print $fl "Warning:File $patric_list in a wrong format for dataset directory.  The final report will be generate without Organism name\n";
			print "Warning:File $patric_list in a wrong format for dataset directory.  The final report will be generate without Organism name\n";
                        $patric_list = undef;
                        last;
		    }
            	    if($_ =~ /genome_id/){}
            	    else{
                    	my @aux = split("\t", $_);
                	$pls{$aux[0]} = $aux[1];
          	    }
       		}
        	close(FILE);
   	    }
        }
        print $fl "\nDataset analysis:\n";
        my $final_dir = $run_dir."/selected";
        system "mkdir $final_dir";
	$discarded_dir = $run_dir."/discarded";
        system "mkdir $discarded_dir";
	close($fl);
	my $partial_head = 0;
	my $final_report = $run_dir."/final_report.csv";
	my $partial_report = $final_dir."/partial_report.csv";
	my $final_fasta = $final_dir."/elements.fasta";
	open(FASTA, ">$final_fasta");
	close(FASTA);
	open(PARTIAL, ">$partial_report");
	if(defined $patric_list){
            print PARTIAL "Organism ID\tOrganism name\tContigID_elem#\t";
        }
        else{
            print PARTIAL "Organism ID\tContigID_elem#\t";
        }
	my $head = 0;
	$caracter = substr $input_dir, -1;
	if($caracter eq '/'){
      	    my $l = length($input_dir);
    	    $input_dir = substr $input_dir, 0, ($l-1);
	}
	opendir(DIR, "$input_dir");
	my @dirs = readdir(DIR);
	closedir(DIR);
	my $hmmsearch_tabular_file;
	my %markers;
	my $str_report = "";
	my %general_markers = ();
	
	# Perform similarity searches for each provided profile HMM
	for my $dir (@dirs){	 
    	    if($dir eq "." or $dir eq ".."){}
            else{
		if(defined $genomes_not{$dir}){
		    my $str = "Genome $dir previously processed. Skipping...\n";
		    open(REP, '>>', $file_log);
		    print REP $str;
		    close(REP);
		    next;
		}	
		my $aux = $input_dir."/".$dir;
        	if(!-d $aux){}
        	else{
		    my $out_dir = $all_results_dir."/".$dir;
		    if(!-e $out_dir){
		    	system "mkdir $out_dir"; 
		    }
		    my $log_temp = $out_dir."/$dir.log";
                    open(LOG, ">$log_temp");
		    print LOG "\nDirectory: $dir\n";
		    my $file = $input_dir."/".$dir."/".$dir.".".$extension_file;
		    if(!-e $file){
			print LOG "\nFile $file not found!\n\n";
			next;
		    }
		    my $type = verifiesFastaFile($file); # Check the FASTA file content
            	    if($type == 1){ # DNA FASTA file
                     	print LOG "Type: Nucleotide Fasta\n";
			my $aux_dir = $input_dir."/".$dir;
			my $aux_name = $input_dir."/".$dir."/".$dir;
			my $msg;
			
			# Run transeq program
                	($transeq_name, $msg) = runTranseq($aux_name, $file, $aux_dir);			
			print LOG $msg;
                    }
            	    else{# Protein FASTA file
               		 print LOG "Type: Protein Fasta\n";
                	$transeq_name = $file;
            	    }
		    if(!-e $transeq_name){
			print LOG "$transeq_name not found. Skipping...\n";
			print STDERR "$transeq_name not found. Skipping...\n";
			close(LOG);
                    	if(-e $log_temp){
                            system "cat  $log_temp >> $file_log";
                            system "rm $log_temp";
                    	}
			next;
		    }
		    open(FILE, $file);
            	    my $seq = "";
		    my %genomes = ();
            	    my $name = undef;
            	    while(<FILE>){
                	chomp($_);
                	if($_ =~ /^>/){
                    	    if(defined $name){
                            	$genomes{$name} = $seq;
                            	$seq = "";
                    	    }
                    	    my @aux = split(" ", $_);
                    	    $name = $aux[0];
                    	    $name =~ s/>//g;
                    	}
                	else{
                      	    $seq .= $_;
                	}
            	    }
            	    close(FILE);
            	    $genomes{$name} = $seq;
		    my $prefix_dir;
		    my $resp;
        	    foreach my $hmm (@hmms_db){
            		%hmms_score = ();
            		my $len = length($hmm);
            		my $aux_name = substr $hmm, (rindex($hmm, "/") + 1), $len;
            		my $aux_prefix = substr $aux_name, 0, rindex($aux_name, ".");
		 	my @aux_marker = split("_", $aux_prefix);
                    	$prefix_dir = $aux_marker[0];
            		my $aux_dir_hmm = $out_dir."/".$prefix_dir; 
			if(!-e $aux_dir_hmm){
            		    system "mkdir $aux_dir_hmm";
			}
			my $msg;
			
			# Check if a cutoff score must be used in similarity searches
            		($resp, $msg) = verifyScoreEvalue($hmm, $aux_dir_hmm);
			print LOG $msg;
			my $aux = $out_dir."/".$prefix_dir."/".$dir; 

			# Run hmmsearch program
            		($hmmsearch_tabular_file,$msg) = runHmmsearch($aux, $transeq_name, $resp, $hmm, $aux_dir_hmm);
			print LOG $msg;
			if($hmmsearch_tabular_file eq "-1"){
			    next;
			}

			# Save similarity searches results in a file
			print "Analysing hmmsearch results...\n";
			my $table_1 = "$aux_dir_hmm/table1.csv";
			my $str_table1 = analyseHmmsearchResults($hmmsearch_tabular_file, $resp);
		 	if(!$str_table1 eq ""){
                	    open(TBL1, ">$table_1") or die "ERROR: Could not create file $table_1!\n";
                	    print TBL1 "Target Name\tquery_pHMM\tE-value\tScore\n";
                	    print TBL1 $str_table1;
                	    close(TBL1);
            		}
			my $aux_hmms_dir = $aux_dir_hmm."/hmms";
                	if(-e $aux_hmms_dir){
                   	    system "rm -rf $aux_hmms_dir";
                	}
		    }
		    opendir(DIR, "$out_dir");
        	    my $abs = abs_path($out_dir);
        	    my @dirs = readdir(DIR);
        	    closedir(DIR);
		    my %models = %{selectModels(\@dirs, $abs)};
		    my %content = ();
        	    my %contigs = ();

		    # Analysis of the similarity searches results found by all profile HMMs

        	    foreach my $key (sort keys %models){
            		my %hash = %{$models{$key}};
            		foreach my $key2 (sort keys %hash){
                	    push @{$content{$key2}}, $hash{$key2};
                	    push @{$contigs{$key2}}, $key;
            		}	
        	    }
		
		# Check if the results found for each input sequence meet the criteria chosen in
	    # the user-defined parameters 

		    foreach my $key (sort keys %contigs){
            		my @coord = ();
            		my $amount = scalar(@{$contigs{$key}});
            		if($amount >= $minimum_amount){
			    my @possible_proteins = ();
                	    my $dir_final = $final_dir."/".$dir;
			    if(!-e $dir_final){
                	    	system "mkdir $dir_final";
			    }
                	    my $result_name = $hmmsearch_tabular_file;
               		    my $len = length($result_name);
                	    my $aux_name = substr $result_name, (rindex($result_name, "/") + 1), $len;
                	    my $prefix = substr $aux_name, 0, rindex($aux_name, ".");
                	    my @contend_marker = @{$content{$key}};
                	    my $temp = $dir_final."/blast_temp_".$key.".fasta";
                	    open(FILE, ">$temp");
               		    print FILE ">$key\n$genomes{$key}\n";	
			    close(FILE);				    
			    my @elements = ();
			    my $seq = $genomes{$key};
			    my %real_orf;
			    %markers = ();
			
			    # Map the genes found in the original sequence to obtain their coordinates
			    for (my $i = 0; $i < scalar(@contend_marker); ++$i){			
                    		my @aux_marker = split("_", $contigs{$key}[$i]);
                   		my $marker = $aux_marker[0];
                    		my @aux = @{$contend_marker[$i]};
				%real_orf = ();
				@aux = @{sortByScoreMax(\@aux)};
				for(my $j = 0; $j < scalar(@aux); ++$j){							
				    my @a = @{$aux[$j]};			       
                    	 	    my $file = $a[3]."/".$prefix.".txt";
                    		    my $subject = $a[2];
                    		    my $fasta = $dir_final."/".$key."_".$marker.".fasta";				
                    		    generateFastaFileByHmmsearch($file ,$a[1],  $marker, $subject, $fasta);
				    if(!-e $fasta){
				    	next;
				    } 
                    		    my $blast = $dir_final."/".$key."_".$marker."_blast.txt";
                    	 	    runBlast($fasta, $temp, $blast);		
				    my @coord_fasta = @{extractCoordinates($blast, $key)};				   
				    if(-e $blast){
                                    	system "rm $blast";
                                    }
                                    if(-e $fasta){
                                    	system "rm $fasta";
                                    }
				    if(scalar(@coord_fasta) > 0){
					for(my $index = 0; $index < scalar(@coord_fasta); ++$index){
					    my @aux = split("\t", $coord_fasta[$index]);
					    my $real = findORF($aux[0], $aux[1], $aux[2],$genomes{$key});					    
					    if(defined $real){
					    	$real .= "\t".$marker;
                                            	$real_orf{$real} = 1;
						if(!defined $markers{$marker}){
                                        	    $markers{$marker} = 1;
                                    		}
						if(!defined $general_markers{$marker}){
                                                    $general_markers{$marker} = 1;
                                                }
					    }
					}
				    }
				    else{
					next;
				    }
				}
				foreach my $key (sort keys %real_orf){
				    push @possible_proteins, $key;
				}
                	    }
			    my @aux_markers = keys %markers;		    
                            my $num = scalar (@aux_markers);

		# If the number of genes is lower than the minimum number defined by the user, the
		# results are stored in a specific directory
		
                            if($num < $minimum_amount){
				if(-e $temp){
                                     system "rm $temp";
                                }
				my $disc = $discarded_dir."/".$dir;
                                if(!-e $disc){
                                    system "mkdir $disc";
                                }
                            	if($num == 0){
                                    print LOG "$dir discarded - it has not markers in any contig.\n";
				    my $temp = $disc."/".$key.".fasta";
                                    open(FILE, ">$temp");
                                    print FILE ">$key\n$genomes{$key}\n";
                                    close(FILE);
                                }
                                else{
                                    my $str_m = "";
                                    if($num > 2){
                                     	my $m;
                                        for($m = 0; $m < $num - 2; ++$m){
                                            $str_m .= $aux_markers[$m].", ";
                                        }
                                        $str_m .= $aux_markers[$m]." and ".$aux_markers[$m+1];
                                    }
                                    elsif($num == 2){
                                        $str_m .= $aux_markers[0]." and ".$aux_markers[1];
                                    }
                                    elsif($num == 1){
                                        $str_m = $aux_markers[0];
                                    }
                                    print LOG "Sequence $key positive to $str_m discarded - number of markers is lower than $minimum_amount\n";
                                    my $temp = $disc."/".$key.".fasta";
                                    open(FILE, ">$temp");
                                    print FILE ">$key\n$genomes{$key}\n";
                                    close(FILE);
				    foreach my $real (sort keys %real_orf){
					my @coord;
                                        push @coord, $real;
                                        generateProteinFasta2(\@coord, $genomes{$key}, $disc, $key);
				    }
				    opendir(DIRHANDLE, $dir_final);
                        	    my @fna = grep {  -f "$dir_final/$_" } readdir DIRHANDLE;
                        	    closedir DIRHANDLE;

				    my $num_files = scalar(@fna);
				    if($num_files == 0){
					system "rm -rf $dir_final";
				    }				     
                                }
			   	next;
                             }			    
			    
			    if($partial_head == 0){
				print PARTIAL "Contig size\tElement coord. on contig\tElement size\t";
				my $len = scalar(@hmms_db);
				for(my $h = 0; $h < $len; ++$h){
            			    if($h < $len - 1){
                			print PARTIAL "Gene\tContig coord.\tElement coord.\tOrientation\tDistance(bases)\t";
            			    }
            			    else{
                			print PARTIAL "Gene\tContig coord.\tElement coord.\tOrientation\t";
            			    }
        			}
				print PARTIAL "\n";
			    }
			    ++$partial_head; 
			    if(-e $temp){
			    	system "rm $temp";       
			    }
			    @possible_proteins = @{sortCoordinates(\@possible_proteins)};
			    my @groups = ();
			    my $disc = $discarded_dir."/".$dir;
			    my $org_name = undef;
			    if(defined $patric_list){
				if(defined $pls{$dir}){
                            	    $org_name = $pls{$dir};
				}
				else{
				    $org_name = "Not defined";
			   	}
			    }
			    my $index = 1;
			    my $pp = scalar(@possible_proteins);
   
			    # Check if the gene order and orientation meet the criteria defined in the 
		    	# parameters
		    	
			    for(my $i = 0; $i < scalar(@possible_proteins) - 1; ++$i){
				my @aux = split("\t", $possible_proteins[$i]);
				my @aux2 = split("\t", $possible_proteins[$i+1]);
				my $distance = $aux2[0] - $aux[1];	
				push @groups, $possible_proteins[$i];		
				if($synteny){
				    my @dists = split(",", $synteny);
				    for(my $j = 1; $j < scalar(@dists); $j = $j + 2){					
				    	if($distance > $dists[$j]){				
					    @groups = @{discard_redundant_genes(\@groups)};	
					    my ($resp, $str, $current_element, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir_final, $dir, $org_name, ($index));
                                            if($resp == 1){
						@groups = ();
                                            	next;
                                            }	
					    @groups = ();
					    $index = $id;
					    ++$index;
                                            $str_report .= $str;
					}
					elsif(($distance <= $dist) and ($i == scalar(@possible_proteins) - 2)){
                                            push @groups, $possible_proteins[$i+1];
					    @groups = @{discard_redundant_genes(\@groups)};
                                            my $size_g = scalar(@groups);
                                            my ($resp, $str, $current_element, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir_final, $dir, $org_name, ($index));
                                            if($resp == 1){
					      	@groups = ();
                                            	next;
                                            }
					    @groups = ();
					    $index = $id;
					    ++$index;
                                            $str_report .= $str;
                                    	}
                                    	elsif(($distance > $dist) and (scalar(@groups) < $minimum_amount)){
                                            @groups = ();
                                    	}
                                    	elsif(($distance < 0) and (abs($distance) > $overlap)){
                                            ++$i;
                                    	}
				    }
				}
				else{
				    if(($distance > $dist) and (scalar(@groups) >= $minimum_amount)){			
					my $size_g = scalar(@groups);
					@groups = @{discard_redundant_genes(\@groups)};
				    	my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir_final, $dir, $org_name, ($index));
					if($resp == 1){
					    @groups = ();
					    next;
					} 
					@groups = ();
					$index = $id;
					++$index;
					$str_report .= $str;					
				    }
				    elsif(($distance <= $dist) and ($i == scalar(@possible_proteins) - 2)){
					push @groups, $possible_proteins[$i+1];
					@groups = @{discard_redundant_genes(\@groups)};
                                        my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir_final, $dir, $org_name, ($index));
                                        if($resp == 1){
					    @groups = ();
                                            next;
                                        }
					@groups = ();
					$index = $id;
					++$index;
                                        $str_report .= $str;
				    }
				    elsif(($distance > $dist) and (scalar(@groups) < $minimum_amount)){
					@groups = ();
				    }	
				    elsif(($distance < 0) and (abs($distance) > $overlap)){
					if($aux[3] ne $aux2[3]){
					    ++$i;		
					}			
				    } 
				    else{
					push @groups, $possible_proteins[$i];
				    }	
				}
			    }
			    my $g = scalar(@groups);
			    if(scalar(@groups) >= $minimum_amount){
				my $i = scalar(@possible_proteins) - 1;
				push @groups, $possible_proteins[$i];
				@groups = @{discard_redundant_genes(\@groups)};
				my ($resp, $str, $id) = analiseGroupOfGenes(\@groups, $disc, $key, $genomes{$key}, $dir_final, $dir, $org_name, $index);
                                if($resp == 1){
                                    next;
                                }
			        $str_report .= $str;
			    }
			    opendir(DIRHANDLE, $dir_final);
                            my @fna = grep {  -f "$dir_final/$_" } readdir DIRHANDLE;
                            closedir DIRHANDLE;
			    my $num_fna = scalar(@fna); 
			    if($num_fna == 0){
				system "rm -rf $dir_final";
			    }
			}
			else{ # The minimum number of genes was not found in the sequence. The results
				  # will be stored in a specific directory
				  
			    my @contend_marker = @{$content{$key}};
                    	    my $markers_str = "";
                    	    for (my $i = 0; $i < scalar(@contend_marker); ++$i){
                        	my @aux_marker = split("_", $contigs{$key}[$i]);
                        	my $marker = $aux_marker[0];
                        	if($i == scalar(@contend_marker)-2){
                            	    $markers_str .= $marker." and ";
                        	}
                        	elsif($i == scalar(@contend_marker)-1){
                            	    $markers_str .= $marker;
                        	}
                        	else{
                           	     $markers_str .= $marker.", ";
                        	}
                    	    }
			    if($amount == 0){
				print LOG "$dir discarded - it has not markers in any contig.\n";
			    }
			    else{
                    	    	print LOG "Sequence $key positive to $markers_str discarded - number of markers is lower than $minimum_amount\n";
			    	my $disc = $discarded_dir."/".$dir;
			    	if(!-e $disc){
				    system "mkdir $disc";
			    	}
			    	my $temp = $disc."/".$key.".fasta";
                            	open(FILE, ">$temp");
                            	print FILE ">$key\n$genomes{$key}\n";
                            	close(FILE);
                            	my $result_name = $hmmsearch_tabular_file;
                           	my $len = length($result_name);
                            	my $aux_name = substr $result_name, (rindex($result_name, "/") + 1), $len;
                            	my $prefix = substr $aux_name, 0, rindex($aux_name, ".");
                            	my @contend_marker2 = @{$content{$key}};
                         	for (my $i = 0; $i < scalar(@contend_marker2); ++$i){
                                    my @aux_marker = split("_", $contigs{$key}[$i]);
                                    my $marker = $aux_marker[0];
                                    my @aux = @{$contend_marker[$i]};
				    for(my $j = 0; $j < scalar(@aux); ++$j){
                                    	my @a = @{$aux[$j]};
                                    	my $file = $a[3]."/".$prefix.".txt";
                                    	my $subject = $a[2];
                                    	my $fasta = $disc."/".$key."_".$marker.".fasta";
                                    	generateFastaFileByHmmsearchDiscarded($file ,$a[1],  $marker, $subject, $fasta);				
				    	if(!-e $fasta){
                                    	    next;
                                    	}
                                    	my $blast = $disc."/".$key."_".$marker."_blast.txt";
                                    	runBlast($fasta, $temp, $blast);
					my @coord_fasta = @{extractCoordinates($blast, $key)};
                                    	if(scalar(@coord_fasta) > 0){
                                            for(my $index = 0; $index < scalar(@coord_fasta); ++$index){
                                            	my @aux = split("\t", $coord_fasta[$index]);
                                            	my $real = findORF($aux[0], $aux[1], $aux[2], $genomes{$key});
						if(defined $real){
                                            	    $real .= "\t".$marker;
						    my @coord;
                                            	    push @coord, $real;
                                            	    generateProteinFasta2(\@coord, $genomes{$key}, $disc, $key);
						    if(-e $blast){
                                            	    	system "rm $blast";
                                        	    }
                                        	    if(-e $fasta){
                                            	    	system "rm $fasta";
                                        	    }
						}
                                            }
                                    	}
                                    	else{
                                            next;
                                    	}
			    	    }
				}
			    }
                	}
        	    }
		    close(LOG);
		    if(-e $log_temp){
		    	system "cat  $log_temp >> $file_log";
        	    	system "rm $log_temp";
		    }
		}
	    } 
	}

	# Generate the final files containing the positive sequences
	
	my $lenght = keys %markers;
        my $count = 1;
        my $str_aux = "";
        open(FINAL, ">$final_report");
	if(defined $patric_list){
	    print FINAL "Organism ID\tOrganism name\tContigID_elem#\t";
	}
	else{
            print FINAL "Organism ID\tContigID_elem#\t";
	}
	print FINAL "Contig size\tElement coord. on contig\tElement size\t";
        my $len = scalar(@hmms_db);
	for(my $h = 0; $h < $len; ++$h){
            if($h < $len - 1){
            	print FINAL "Gene\tContig coord.\tElement coord.\tOrientation\tDistance(bases)\t";
            }
            else{
                print FINAL "Gene\tContig coord.\tElement coord.\tOrientation\t";
            }
        }
        foreach my $key (sort keys %general_markers){
	    my $aux_file = $final_dir."/saida";
            system "ls $final_dir/*/*.fna >> $aux_file 2> $aux_file";
            open(AUX, $aux_file);
            my $num = 1;
            while(<AUX>){
            chomp($_);
                if($_ =~ /ls:/){
                    $num = 0;
                    last;
                }
            }
            close(AUX);
            system "rm $aux_file";
		
	    if($num == 0){
		next;
	    }
            my $protein_file = $final_dir."/".$key."_proteins.fasta";
            my $command = "cat $final_dir/*/*$key\*_protein.fasta > $protein_file 2> /dev/null";
            system "$command";
            $count++;	  
            if(!-z $protein_file){
	    	my %sequences = ();
	    	open(FILE, $protein_file);
	    	my $name = undef;
	    	my $seq = "";
	    	while(<FILE>){
		    chomp($_);
		    if($_ =~ /^>/){
			if($seq ne ""){
			    $seq =~ s/\*//g;
			    push @{$sequences{$name}}, $seq;
		    	}
		    	$name = $_;
		    	$seq = "";
		    }
		    else{
		    	$seq .= $_;
		    }
	    	}
		$seq =~ s/\*//g;
		push @{$sequences{$name}}, $seq;
	    	close(FILE);
	    	my $aux_protein = $final_dir."/".$key."_proteins_temp.fasta";
	    	open(FILE, ">$aux_protein");
	    	foreach my $key (sort {$a cmp $b} keys %sequences){
		    my @aux_seqs = @{$sequences{$key}};
		    for(my $k = 0; $k < scalar(@aux_seqs); ++$k){
		    	print FILE "$key\n$aux_seqs[$k]\n";
		    }
	    	}
	    	close(FILE);
	    	system "mv $aux_protein $protein_file >> /dev/null";
            }
	    else{
	    }
        }
        if($head == 0 and $str_aux ne ""){
            print FINAL $str_aux;
            $head = 1;
        }
        print FINAL "\n$str_report";
	close(FINAL);
	close(PARTIAL);
	system "rm $partial_report";

	opendir(DIR, "$final_dir");
        my @final_files = readdir(DIR);
        closedir(DIR);
	if(-e "$final_dir/*/*blast_temp.txt"){
	    system "rm $final_dir/*/*blast_temp.txt";
	}
	foreach my $fasta_file (@final_files){
	    if($fasta_file eq "." or $fasta_file eq ".."){}
            else{
		my $aux = $final_dir."/".$fasta_file;
                if(!-d $aux){}
		else{		   
		    my $aux_file = $final_dir."/saida";
                    system "ls $aux/*.fna >> $aux_file 2> $aux_file";
                    open(AUX, $aux_file);
                    my $num = 1;
                    	while(<AUX>){
                            chomp($_);
                            if($_ =~ /ls:/){
                                $num = 0;
                                last;
                            }
                        }
                    	close(AUX);
                        system "rm $aux_file";
		    if($num > 0){ 
		    	my $command = "cat $aux/*.fna >> $final_fasta";
		    	system $command;
		    }
		    else{
			if(-e $aux){
			    system "rm -rf $aux";
			}
		    }
		}
	    }
	}
  	if(!-z $final_fasta){
	    my %sequences = ();
            open(FILE, $final_fasta);
            my $name = undef;
            my $seq = "";
            while(<FILE>){
            	chomp($_);
            	if($_ =~ /^>/){
            	    if($seq ne ""){                    
                    	push @{$sequences{$name}}, $seq;
                    }
                    $name = $_;
                    $seq = "";
            	}
            	else{
                    $seq .= $_;
            	}
      	    }
            push @{$sequences{$name}}, $seq;

            close(FILE);
            my $aux_fasta = $final_dir."/elements_temp.fasta";
            open(FILE, ">$aux_fasta");
            foreach my $key (sort {$a cmp $b} keys %sequences){
		my @aux_seqs = @{$sequences{$key}};
		for(my $k = 0; $k < scalar(@aux_seqs); ++$k){
            	    print FILE "$key\n$aux_seqs[$k]\n";
		}
            }
            close(FILE);
            system "mv $aux_fasta $final_fasta >> /dev/null";	
    	}
	else{
	    system "rm $final_fasta >> /dev/null";
  	}
    }
}

print "Done.\n";
close($log);
if(-e $output and -e $log){
    system "mv $log $output >> /dev/null";
}
exit;

####################################################################
###                        Sobroutines                           ###
####################################################################

# This routine checks if a directory with the given output name already 
# exists. If so, a numeric suffix is added to the name.

sub output_dir_name {
    my $output_dir_name = shift;
    my $count = 2;
    my $flag = 0;
    print STDERR "Creating output directory $output_dir_name\n";
    if (-d $output_dir_name) {
	my $len = length($output_dir_name);
	my $dir_name = substr $output_dir_name, 0, rindex($output_dir_name, "/");
	my $subdir_aux = substr $output_dir_name, (rindex($output_dir_name, "/") + 1), $len;
	my @aux = split("_", $subdir_aux);
	my $subdir_name = $aux[0];
        $flag = 1;
	while (-d "$dir_name/$aux[0]\_$count") {
            $count++;
        }
        $output_dir_name = "$dir_name/$aux[0]\_$count"; 
    }
    print STDERR "\nOutput directory already exists, saving results to $output_dir_name instead.\n\n" if $flag;
    return ($output_dir_name);
}

# This routine identifies if the FASTA file contains nucleotide or protein sequences

sub verifiesFastaFile{
    my $file = shift;
    my $count = 0;
    my $num_lines = 0;
    my $type = 1; # 1 - nucleotide; 2 - protein
    my $first = 0;
    open(my $vf, "$file");
    while(<$vf>){
    	chomp($_);
	if($_ =~ /^>/){
	    if($first == 0){
                $first = 1;
            }
            else{
                if($count > 0){
                    $type = 2;
                }
                last;
            }
	}
        else{
	    if ((index($_, 'L') != -1) or (index($_, 'l') != -1)) {
	    	++$count;
	    }
	    if ((index($_, 'V') != -1) or (index($_, 'v') != -1)) {
            	++$count;
            }
	    if ((index($_, 'I') != -1) or (index($_, 'i') != -1)) {
            	++$count;
            }
	    if ((index($_, 'P') != -1) or (index($_, 'p') != -1)) {
            	++$count;
            }
	    if ((index($_, 'F') != -1) or (index($_, 'f') != -1)) {
            	++$count;
            }
	    if ((index($_, 'S') != -1) or (index($_, 's') != -1)) {
            	++$count;
            }
	    if ((index($_, 'Q') != -1) or (index($_, 'q') != -1)) {
            	++$count;
            }
	    if ((index($_, 'D') != -1) or (index($_, 'd') != -1)) {
            	++$count;
            }
	    if ((index($_, 'E') != -1) or (index($_, 'e') != -1)) {
            	++$count;
            }
	    if ((index($_, 'R') != -1) or (index($_, 'r') != -1)) {
            	++$count;
            }
	    if ((index($_, 'H') != -1) or (index($_, 'h') != -1)) {
            	++$count;
            }
	    if ((index($_, 'M') != -1) or (index($_, 'm') != -1)) {
            	++$count;
            }
        }
    }    
    close($vf);
    return $type;
}

# This routine checks if the user-specified tabular input file is a hmmsearch tabular
# output file.

sub verificaArquivoTabular{
    my $file = shift;
    open(my $tab, "$file");
    while(<$tab>){
        chomp($_);
	if($_ =~ /target name/){
	    my @line = split(" ", $_);
	    if(($line[3] eq "accession") and ($line[4] eq "query") and ($line[6] eq "accession") and ($line[7] eq "E-value")){
		close($file);
		return 1;
	    }
	}
    }
    close($tab);
    return 0;
}

# This routine runs transeq program.
sub runTranseq{
    my $name = shift;
    my $file = shift;
    my $dir = shift;
    my $string;
    my $aux;
    if(defined $dir){
    	$string = $dir."/".$name."*_transeq*";
	$aux = `ls $string 2> /dev/null`;
	if($aux eq ""){
	    $string = $name."*_transeq*";
	    $aux = `ls $string 2> /dev/null`;
	}
    }
    else{
	$string = $name."*_transeq*";
	$aux = `ls $string 2> /dev/null`;	
    }

    if($aux eq ""){	
	$aux = undef;
    }
    else{
	my @aux_t = split("\n", $aux);
	$aux = $aux_t[0];
    }
    my $transeq;
    my $auxiliar;
    my $msg;
    $msg = "Step: translation of nucleic acid sequences\n";
    if(defined $aux){ #Verify if transeq file exists
	$msg .= "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	print "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	$transeq = $aux;
    }
    else{ #run the transeq program:
     	$transeq = $name."_transeq.fasta";
	$msg .= "Performing conceptual translation of $file and creating $transeq (Parameters: frame = 6) \n";
        print "Performing conceptual translation of $file and creating $transeq... \n";
        system "transeq -frame 6 $file $transeq";
        print "Done.\n";
	$msg .= "Done\n";
    }
    return ($transeq, $msg);
}

# This routine runs hmmsearch program.

sub runHmmsearch{
    my $name = shift;
    my $transeq = shift;
    my $option = shift;
    my $hmm = shift;
    my $dir = shift;
    my $file1 = $name."_hmmsearch.tab";
    my $full = $name."_hmmsearch.txt";
    my $msg;
    $msg = "Step: similarity search with hmmsearch.\n";
    if(-e $file1){
	my $type = verificaArquivoTabular($file1);
        if($type == 1){
	    $msg .= "Similarity search file found. Skipping similarity search with hmmsearch.\n";
            print "Similarity search file found. Skipping similarity search with hmmsearch.\n";
	}
	else{
	    $msg .= "Step: similarity search with hmmsearch.\n";	    
	    die "Similarity search file found, but it is not in a tabular format.\n";
	}
	
    }
    else{
	$msg .= "Performing similarity search with hmmsearch (Parameters: -E 10 --tblout $file1 -o $full --cpu $cpu $hmm $transeq)\n";
        print "Performing similarity search with hmmsearch ($hmm)... \n";			
        my $error = !system "hmmsearch -E 10 --tblout $file1 -o $full --cpu $cpu $hmm $transeq";
	if($error != 1){
	    $msg .= "ERROR: Could not run hmmsearch : $!\nCommand: hmmsearch --tblout $file1 -o $full $hmm $transeq. Skipping...\n";
	    $file1 = "-1";
	}
	else{
	    print "Done\n";
	    $msg .= "Done\n";
	}
    }
    return ($file1, $msg);
}

# This routine generates FASTA files from the hmmsearch results

sub generateFastaFileByHmmsearch{
    my $file = shift;
    my $hmm_name = shift;
    my $marker = shift;
    my $subject = shift;
    my $output_file = shift;
    my $found = 0;
    my $sub = 0;
    my $final = 0;
    my $score_domain = -100;
    my $get_seq = 0;
    my $found_sub = 0;
    my $final_seq = "";
    open(FILE, $file);
    while(<FILE>){
  	chomp($_);
        if($_ =~ /Query:/){
	    if($_ =~ /$hmm_name/){
            	$found = 1;
            }
	}
	elsif($_ =~ />>/){
            if($found == 1  and $_ =~ /$subject/){
            	$found_sub = 1;
            }
        }
	elsif($found == 1 and $found_sub == 1){
	    if($_ =~ /==/){
		my @aux_domain = split(" ", );
		my $aux_score = $aux_domain[4];
		if($aux_score > $score_domain){
		    $score_domain = $aux_score;
		    $get_seq = 1;
		    $final_seq = "";
		}
		else{
		    $get_seq = 0;
		}
	    }
	    elsif($get_seq == 1 and $_ =~ /$subject/){
		my $seq = $_;
            	$seq =~ s/$subject//;
		if($seq =~ /(\d+)(\s)((\w|\-|\*)+)(\s)(\d+)/){
                    my $aux2 = $3;
                    $aux2 =~ s/-//ig;
		    $final_seq .= $aux2;
		}
	    }
	    elsif($_ =~ /Internal pipeline statistics summary/ and $final_seq ne ""){
		open(OUT, ">$output_file");
                print OUT ">$marker\n$final_seq\n";
                close(OUT);
		$found = 0;
		$found_sub = 0;
	    }
	}
    }

    close(FILE);
}

# This routine generates FASTA files from the hmmsearch results (discarded results).

sub generateFastaFileByHmmsearchDiscarded{
    my $file = shift;
    my $hmm_name = shift;
    my $marker = shift;
    my $subject = shift;
    my $output_file = shift;
    my $found = 0;
    my $sub = 0;
    my $final = 0;
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
        if($_ =~ /Query:/){
            if($_ =~ /$hmm_name/){
                $found = 1;
            }
        }
        elsif($_ =~ />>/){
            if($found == 1 and $final == 0){
                if($_ =~ /$subject/){
                     $sub = 1;
                }
            }
        }
        elsif($sub == 1 and $_ =~ /$subject/ and $final == 0){
            $final = 1;
            my $seq = $_;
            $seq =~ s/$subject//;
            if($seq =~ /(\d+)(\s)((\w|\-|\*)+)(\s)(\d+)/){
                my $aux2 = $3;
                $aux2 =~ s/-//ig;
                my $len = length($aux2);
                open(OUT, ">$output_file");
                print OUT ">$marker\n$aux2\n";
                close(OUT);
            }
        }
    }
    close(FILE);
}

# This routine runs tblastn program

sub runBlast{
    my $query = shift;
    my $subject = shift;    
    my $output_file = shift;
    my $command = "tblastn -query $query -subject $subject -outfmt '6 qseqid qlen qstart qend sseqid sstart send sframe score' -out $output_file >> /dev/null";
    system "$command";     
}

# This routine runs blastn program

sub runBlastn{
    my $query = shift;
    my $subject = shift;
    my $output_file = shift;
    my $command = "blastn -query $query -subject $subject -outfmt '6 qseqid qlen qstart qend sseqid sstart send sframe score' -out $output_file >> /dev/null";
    system "$command";
}

# This routine removes redundant genes found by the similarity searches

sub removeEqualsCoordinates{
    my $aux_vector = shift;
    my @vector = @{$aux_vector};
    my $len = scalar(@vector);
    my %hash = ();
    for(my $i = 0; $i < $len; ++$i){
	my @aux = split("\t", $vector[$i]);
        my $marker = $aux[3];
	push @{$hash{$marker}}, $vector[$i]; 
    }
    my @return = ();
    foreach my $key (sort keys %hash){
	my @aux = @{$hash{$key}};
	my %aux_hash = ();
	for(my $i = 0; $i < scalar(@aux); $i++){
	    my @aux2 = split("\t", $aux[$i]);
	    my $start = $aux2[0];
	    if(!defined $aux_hash{$start}){
		$aux_hash{$start} = 1;
		push @return, $aux[$i];
	    }
	}
    }
    for(my $i = 0; $i < scalar(@return); $i++){
	my @aux = split("\t", $return[$i]);
    }
    return (\@return);
}

# This routine sorts the genes found according to their coordinates

sub sortCoordinates{
    my $aux_vector = shift;
    my @vector = @{$aux_vector};
    my $len = scalar(@vector);
    for(my $i = 0; $i < $len - 1; ++$i){
        my @aux = split("\t", $vector[$i]);
        my $start = $aux[0];
	my $index = $i;
	my $lower = $start;
        for(my $j = $i+1; $j < $len; ++$j){
            my @aux2 = split("\t", $vector[$j]);
            my $start2 = $aux2[0];	
	    if($lower > $start2){
		$lower = $start2;
		$index = $j; 
	    }
        }
	if($index != $i){
	    my $aux3 = $vector[$i];
            $vector[$i] = $vector[$index];
            $vector[$index] = $aux3;
        }
    }
    return \@vector;
}

# This routine sorts the hmmsearch results by their scores

sub sortByScoreMax{
    my $aux_vector = shift;
    my @vector = @{$aux_vector};
    my $len = scalar(@vector);
    for(my $i = 0; $i < $len - 1; ++$i){
        my @aux = @{$vector[$i]};
        my $score = $aux[0];
        my $index = $i;
        for(my $j = $i+1; $j < $len; ++$j){
            my @aux2 = @{$vector[$j]};
            my $score2 = $aux2[0];
            if($score2 > $score){
                $score = $score2;
                $index = $j;
            }
        }
        if($index != $i){
            my $aux3 = $vector[$i];
            $vector[$i] = $vector[$index];
            $vector[$index] = $aux3;
        }
    }
    return \@vector;
}

# This routine checks if the provided profiles HMMs contain the CUTOFF SCORE tag

sub verifyCutoff{
    my $hmm = shift;
    my $value = `grep SCORE $hmm`;
    $value =~ s/CUTOFF SCORE//;
    $value =~ s/\s//g;
    $value =~ s/\://g;
    if((defined $value) and !($value eq "")){
        return $value;
    }
    else{
        return undef;
    }
}

# This routine checks if a score or e-value cutoff value was provided by the user

sub verifyScoreEvalue{
    my $hmm = shift;
    my $dir = shift;
    my $aux_score = verifyCutoff($hmm);
    my $resp = 0;
    my $msg = "";
    if(defined $aux_score){
        $ignore_cutoff =~ s/\s//g;
        if(lc($ignore_cutoff) eq "no"){	    
            my $number = `grep -c NAME $hmm`;
            if($number == 1){
                $usr_score = $aux_score;
                $usr_evalue = undef;
                $msg .= "HMM cutoff score: $usr_score\n";
                return (1, $msg);
            }
            else{
                my $dir = $dir."/hmms";
                system "mkdir $dir";
                $/ = "//\n";
                open(DATA, "<$hmm");
                foreach my $unique_hmm (<DATA>){
                    my @hmm_lines = split("\n", $unique_hmm);
                    my $aux_n;
                    my $achei = 0;
                    foreach my $line (@hmm_lines){
                        if($line =~ m/(NAME\d*)/g){
                            $aux_n = $line;
                            $aux_n =~ s/NAME//;
                            $aux_n =~ s/\s+//g;
                        }
                        elsif($line =~ m/(CUTOFF SCORE\d*)/g){
                            my $value = $line;
                            $value =~ s/CUTOFF SCORE//;
                            $value =~ s/\s//g;
                            $value =~ s/\://g;
                            $hmms_score{lc($aux_n)} = $value;
                            my $file = $dir."/".$aux_n.".hmm";
                            open(HMMFILE, ">$file") or die "ERROR: Couldn't create file $output/$file : $!\n";
                            print HMMFILE $unique_hmm;
                            close (HMMFILE);
                            $achei = 1;
                            next;
                        }
                   }
                   if($achei == 0){
                        my $file = $dir."/".$aux_n.".hmm";
                        open(HMMFILE, ">$file") or die "ERROR: Couldn't create file $output/$file : $!\n";
                        print HMMFILE $unique_hmm;
                        close (HMMFILE);
                   }
              }
            close (DATA);
            $/ = "\n";
            return (2, $msg);
           }
        }
        elsif(lc($ignore_cutoff) eq "yes"){
            if(defined $usr_evalue){
                $msg .= "Cutoff evalue defined by user: $usr_evalue\n";
            }
            elsif(defined $usr_score){
                $msg .= "Cutoff score defined by user: $usr_score\n";
            }
            return (0, $msg);
        }
    }
    else{
        if(($usr_score) and ($usr_evalue)) {
            system "rm -rf $output";
            die "ERROR: Please use -s or -e, not both.\n$help_print\n";
        }
        if(defined $usr_evalue){
            $msg .="Cutoff evalue defined by user: $usr_evalue\n";
        }
        elsif(defined $usr_score){
            $msg .= "Cutoff score defined by user: $usr_score\n";
        }
        return (0, $msg);
    }
}

# This routine checks whether the synteny of the genes meet the user-defined criteria 

sub verifySynteny{
    my $aux_order = shift;
    my $aux_genes = shift;
    my @order = @{$aux_order};
    my @genes = @{$aux_genes};
    my $size = scalar(@order);
    my $amount;
    if($size % 2 == 0){ 
	$amount = $size/2;
    }
    else{
    	my $float = scalar(@order);
    	$amount = int($float/2)+1; 
    }
    my $limit = scalar(@genes);
    if($amount == $limit){
   	my $index = 0;
        for(my $i = 0; $i < ($limit-1); ++$i){
            my @c1 = split("\t", $genes[$i]);
            my @c2 = split("\t", $genes[$i+1]);    
	    my ($g1, $g2, $resp) = verifyGenesDirection($order[$index], $order[$index+2], $c1[2], $c2[2]);	    
	    if($resp == 1){
		return 4;
	    }
	    else{		   
	     	if(($c1[3] eq $g1) and ($c2[3] eq $g2)){
                    my $end1 = $c1[1];
		    my $start2 = $c2[0];		    
                    my $diff = $start2 - $end1;
                    my $real_diff = $order[$index+1];
		    if($diff < 0){
		    	return 2;
		    }
                    elsif($diff < $real_diff){
                    }
                    else{
                    	return 2; 
                    }
            	}
                else{
                    return 3; 
            	}
            	$index += 2;
	    }
        }
        return 0;
    }
    else{
        my $last_indice = -1;
        for(my $i = 0; $i < $limit-1; ++$i){
            my @c = split("\t", $genes[$i]);
            my $index;
            my $found = 0;
            for(my $j = 0; $j < scalar(@order); $j = $j + 2){
		my @aux = split("", $order[$j]);
                my $aux_gene;
		my $orientation = undef;
                if(($aux[0] eq "+") or ($aux[0] eq "-")){
                    my $len = scalar(@aux)-1;
                    $aux_gene = join("", @aux[1..$len]);
		    $orientation = $aux[0];			
                }
                else{
                    $aux_gene = $order[$j];
                }
                if($aux_gene eq $c[3]){
		    if(((defined $orientation) and ($orientation eq $c[2])) or (!defined($orientation))){
                        $index = $j;
                        $found = 1;
		    }
		    elsif((defined $orientation) and ($orientation ne $c[2])){
			return 4;
		    }
                }
            }
            if($found == 0){
                return 3;
            }
            else{
                if($index < $last_indice){
                    return 3;
                }
                else{
                    my @c1 = split("\t", $genes[$i]);
                    my @c2 = split("\t", $genes[$i+1]);
                    my $end1;
                    if($c1[2] eq "+"){
                        $end1 = $c1[1];
                    }
                    else{
                        $end1 = $c1[0];
                    }
                    my $start2;
                    if($c2[2] eq "+"){
                        $start2 = $c2[0];
                    }
                    else{
                        $start2 = $c2[1];
                    }
                    my $genes_distance = $start2 - $end1;
                    my $limit_distance = 0;
                    my $found2 = 0;
                    for(my $k = $index; $k < scalar(@order)-1; $k = $k + 2){
                        $limit_distance += $order[$k+1];
		        my @aux = split("", $order[$k]);
			my $aux_gene;
    			if(($aux[0] eq "+") or ($aux[0] eq "-")){
            		    my $len = scalar(@aux)-1;
            		    $aux_gene = join("", @aux[1..$len]);
        		}
			else{
			    $aux_gene = $order[$k];
			}
                        if($c2[3] eq $aux_gene){
			    my ($g1, $g2, $resp) = verifyGenesDirection($order[$index], $order[$index+2], $c1[2], $c2[2]);
            		    if($resp == 1){
                		return 4;
            		    }
			    else{
                            	if($genes_distance > $limit_distance){
                                    return 2;
                            	}
				$found2 = 1;
				$last_indice = $k;
                            	last;
			    }
                        }
                    }
                }
            }
        }
        return 0;
    }
}

# This routine checks the gene orientation

sub verifyGenesDirection{
    my $g1 = shift;
    my $g2 = shift;
    my $ori1 = shift;
    my $ori2 = shift;
    my $gene1;
    my $gene2;
    my $dir1 = undef;
    my $dir2 = undef;
    my @aux = split("", $g1);
    if(($aux[0] eq "+") or ($aux[0] eq "-")){
	$dir1 = $aux[0];
	my $len = scalar(@aux)-1;
        $gene1 = join("", @aux[1..$len]);
    }
    else{
	$gene1 = $g1;
    }

    @aux = split("", $g2);
    if(($aux[0] eq "+") or ($aux[0] eq "-")){
        $dir2 = $aux[0];
        my $len = scalar(@aux)-1;
        $gene2 = join("", @aux[1..$len]);
    }
    else{
        $gene2 = $g2;
    }

    if(!(defined $dir1) and !(defined $dir2)){
	return ($gene1, $gene2, 0);
    }
    elsif((defined $dir1) and (defined $dir2)){
	if(($dir1 eq $ori1) and ($dir2 eq $ori2)){
	    return ($gene1, $gene2, 0);
	}
	else{
	    return ($gene1, $gene2, 1);
	}
    }
    elsif((defined $dir1) and !(defined $dir2)){
	if($dir1 eq $ori1){
	    return ($gene1, $gene2, 0);
	}
	else{
            return ($gene1, $gene2, 1);
        }
    }
    else{
	if($dir2 eq $ori2){
            return ($gene1, $gene2, 0);
        }
        else{
            return ($gene1, $gene2, 1);
        }
    }
}

# This routine checks whether the synteny of the genes meet the user-defined criteria for 
# a circular topology 

sub verifyCircularSynteny{
    my $aux_order = shift;
    my $aux_genes = shift;
    my @order = @{$aux_order};
    my @genes = @{$aux_genes};
    my @answers = ();    
    for(my $i = 0; $i < scalar(@order); $i = $i + 2){
  	my @syn = ();
    	push @syn, $order[$i];
    	push @syn, $order[$i+1];
   	my $j = $i + 2;
    	if($j >= scalar(@order)){
            $j = 0;
    	}
   	while($j != $i){
            push @syn, $order[$j];
            push @syn, $order[$j+1];
            $j = $j + 2;
            if($j >= scalar(@order)){
            	$j = 0;
            }
   	}
	my $resp = verifySynteny(\@syn, \@genes);
	push @answers, $resp;
    }
    
    my @aux_order = ();
    my $start;
    my $add = 0;
    my $len = scalar(@order);
    if(($len%2) == 0){
        $start = $len - 2;
        $add = 1;
    }
    else{
        $start = $len - 1;
    }
    for(my $i = $start; $i >= 0; $i = $i - 2){
        my @aux_g = split("", $order[$i]);
        if(($aux_g[0] eq "+") or ($aux_g[0] eq "-")){
            my $signal;
            if($aux_g[0] eq "+"){
            	$signal = "-";
            }
            else{
                $signal = "+";
            }
            my $len_g = scalar(@aux_g)-1;
            my $aux_gene = join("", @aux_g[1..$len_g]);
            push @aux_order, $signal.$aux_gene;
            if($i > 0){
                push @aux_order, $order[$i-1];
            }
        }
        else{
            push @aux_order, $order[$i];
            if($i > 0){
                push @aux_order, $order[$i-1];
            }
        }
        my $t = $i - 1;
    }
    if($add == 1){
        my $t = $len - 1;
        push @aux_order, $order[$t];
    }

    for(my $i = 0; $i < scalar(@aux_order); $i = $i + 2){
        my @syn = ();
        push @syn, $aux_order[$i];
        push @syn, $aux_order[$i+1];
        my $j = $i + 2;
        if($j >= scalar(@aux_order)){
            $j = 0;
        }
        while($j != $i){
            push @syn, $aux_order[$j];
            push @syn, $aux_order[$j+1];
            $j = $j + 2;
            if($j >= scalar(@aux_order)){
                $j = 0;
            }
        }
        my $resp = verifySynteny(\@syn, \@genes);
        push @answers, $resp;
    }

    for(my $i = 0; $i < scalar(@answers); ++$i){
	if($answers[$i] == 0){
	    return 0;
	}
    }   
    for(my $i = 0; $i < scalar(@answers); ++$i){
	if($answers[$i] == 2){
            return 2;
        }
    } 
    return 3;  
}

# This routine checks if the gene order and orientation meet the user-defined criteria for
# a circular topology (without considering the intergenic distance)

sub verifyCircularOrder{
    my $aux_order = shift;
    my $aux_genes = shift;
    my @order = @{$aux_order};
    my @genes = @{$aux_genes};
    my @answers = ();
    for(my $i = 0; $i < scalar(@order); $i = $i + 2){
        my @syn = ();
        push @syn, $order[$i];
        push @syn, $order[$i+1];
        my $j = $i + 2;
        if($j >= scalar(@order)){
            $j = 0;
        }
        while($j != $i){
            push @syn, $order[$j];
            push @syn, $order[$j+1];
            $j = $j + 2;
            if($j >= scalar(@order)){
                $j = 0;
            }
        }
        my $resp = verifySyntenyWithoutDistance(\@syn, \@genes);
        push @answers, $resp;
    }

    my @aux_order = ();
    my $start;
    my $add = 0;
    my $len = scalar(@order);
    if(($len%2) == 0){
        $start = $len - 2;
        $add = 1;
    }
    else{
        $start = $len - 1;
    }
    for(my $i = $start; $i >= 0; $i = $i - 2){
        my @aux_g = split("", $order[$i]);
        if(($aux_g[0] eq "+") or ($aux_g[0] eq "-")){
            my $signal;
            if($aux_g[0] eq "+"){
                $signal = "-";
            }
            else{
                $signal = "+";
            }
            my $len_g = scalar(@aux_g)-1;
            my $aux_gene = join("", @aux_g[1..$len_g]);
            push @aux_order, $signal.$aux_gene;
            if($i > 0){
                push @aux_order, $order[$i-1];
            }
        }
        else{
            push @aux_order, $order[$i];
            if($i > 0){
                push @aux_order, $order[$i-1];
            }
        }
        my $t = $i - 1;
    }
    if($add == 1){
        my $t = $len - 1;
        push @aux_order, $order[$t];
    }

    for(my $i = 0; $i < scalar(@aux_order); $i = $i + 2){
        my @syn = ();
        push @syn, $aux_order[$i];
        push @syn, $aux_order[$i+1];
        my $j = $i + 2;
        if($j >= scalar(@aux_order)){
            $j = 0;
        }
        while($j != $i){
            push @syn, $aux_order[$j];
            push @syn, $aux_order[$j+1];
            $j = $j + 2;
            if($j >= scalar(@aux_order)){
                $j = 0;
            }
        }

        my $resp = verifySyntenyWithoutDistance(\@syn, \@genes);

        push @answers, $resp;
    }

    for(my $i = 0; $i < scalar(@answers); ++$i){
        if($answers[$i] == 0){
            return 0;
        }
    }
    return 1;
}

# This routine analyzes hmmsearch results

sub analyseHmmsearchResults{
    my $file = shift;
    my $resp = shift;    
    open(FILE, $file);
    my $str = "";
    while(<FILE>){
    	chomp($_);
        if ($_  =~ m/^#/){
            next;
        }
        else{
            my @aux = split(" ", $_);
            if($resp == 0 || $resp == 1){
            	if(defined $usr_score){
                    if($aux[5] >= $usr_score) {
                    	$str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                    }
              	}
             	elsif(defined $usr_evalue) {
             	    my $calc_str = sprintf("%.10f", $aux[4]);
                    if( $calc_str <= $usr_evalue) {
                    	$str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                    }
             	}
             	else{
                    $str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
             	}
            }
            elsif($resp == 2){
            	my $name = lc($aux[2]);
                $name =~ s/\s+//g;
                if(defined $hmms_score{$name}){
                     my $aux_sc = $hmms_score{$name};
                     if($aux[5] >= $aux_sc) {
                     	$str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                     }
                }
                elsif(defined $usr_score){
                    if($aux[5] >= $usr_score) {
                    	$str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                    }
                }
                elsif(defined $usr_evalue){
                    my $calc_str = sprintf("%.10f", $aux[4]);
                    if( $calc_str <= $usr_evalue) {
		 	$str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                    }
                }
                else{
                    $str .= "$aux[0]\t$aux[2]\t$aux[4]\t$aux[5]\n";
                }
            }
       	}
    }
    close(FILE);
    return $str;
}

# This routine selects the profile HMMs that showed the best results in similarity searches

sub selectModels{
    my $aux = shift;
    my $abs = shift;
    my @dirs = @{$aux};
    my %mod = ();
    foreach my $object (@dirs){
        if($object eq "." || $object eq ".."){}
        else{
            $object = $abs."/".$object;
            if(-d $object){
                my %hash = ();
                my $len = length($object);
                my $aux_name = substr $object, (rindex($object, "/") + 1), $len;
                my $index = rindex($aux_name, ".");
                my $key;
                if($index > -1){
                    $key = substr $aux_name, 0, $index ;
                }
                $key = $aux_name;
                my $table1 = $object."/table1.csv";
                if(-e $table1){
                    open(FILE, $table1);
                    while(<FILE>){
                        chomp($_);
                        if($_ =~ /Target/){}
                        else{
                            my @aux = split("\t", $_);
                            my $index = length($aux[0])-2;
                            my $contig = substr $aux[0],0, $index;
			    my @aux_hash = ();
			    $aux_hash[0] = $aux[3];
                            $aux_hash[1] = $aux[1];
                            $aux_hash[2] = $aux[0];
                            $aux_hash[3] = $object; 
			    push @{$hash{$contig}}, \@aux_hash;
                        }
                    }
                    close(FILE);
                    $mod{$key} = \%hash;
                }
            }
        }
    }
    return \%mod;
}

# This routine returns the coordinates of the similarity search results

sub extractCoordinates{
    my $file = shift;
    my $key = shift;
    my $start_marker;
    my $end_marker;
    my $frame_marker;
    my @vetor = ();
    open(FL, $file);
    while(<FL>){
    	chomp($_);
        my @aux2 = split(" ",$_);
        my $len = $aux2[1];
	my $start;
	my $end;	
        $start = $aux2[2];
        $end = $aux2[3];
        my $sub = $aux2[4];
        if(((($end - $start) + 1) == $len) and ($sub eq $key)){
            $start_marker = $aux2[5];
            $end_marker = $aux2[6];
            if($aux2[7] > 0){
            	$frame_marker = '+';
            }
            else{
            	$frame_marker = '-';
            }
    	    my $string_coord = $start_marker."\t".$end_marker."\t".$frame_marker;
	    push @vetor, $string_coord;
        }
    }
    close(FL);
    return \@vetor;
}

# This routine returns the gene coordinates according to blast similarity searches against
# the genome sequence

sub extractCoordinatesToFastaFile{
    my $file = shift;
    my $key = shift;
    my $number = shift;
    my $start_marker;
    my $end_marker;
    my $frame_marker;    
    my $aux_line = `wc -l $file`;
    my @aux = split(" ", $aux_line);
    my $line = $aux[0];
    open(FILE, $file);
    my $count = 1;
    my $string_coord;
    while(<FILE>){
        chomp($_);
        my @aux2 = split(" ",$_);
        my $len = $aux2[1];
        my $start;
        my $end;
        $start = $aux2[2];
        $end = $aux2[3];
        my $sub = $aux2[4];
        $start_marker = $aux2[5];
        $end_marker = $aux2[6];
        if($aux2[7] > 0){
            $frame_marker = '+';
        }
        else{
            $frame_marker = '-';
        }
 	$string_coord = $start_marker."\t".$end_marker."\t".$frame_marker;
	if($line == 1){
	    close(FILE);
            return $string_coord;
	}
	if($count == $number){      
            close(FILE);
            return $string_coord;
	}
	++$count;
    }
    close(FILE);
    return $string_coord;
}

# This routine checks the intergenic distance, gene order and orientation

sub verifyCoordinates{
    my $aux = shift;
    my @coordinates = @{$aux}; 
    if($synteny){
	if(lc($circular) eq "yes"){
	    my @order = split(",", $synteny);
	    my $resp = verifyCircularSynteny(\@order, \@coordinates); 
	    return $resp;
	}
	else{
            my @order = split(",", $synteny);
            my $resp = verifySynteny(\@order, \@coordinates);
	    if($resp != 0){
		my @aux_order = ();
		my $start;
		my $add = 0;
		my $len = scalar(@order);
		if(($len%2) == 0){
		    $start = $len - 2;	
		    $add = 1;
		}
		else{
		    $start = $len - 1;
		}
		for(my $i = $start; $i >= 0; $i = $i - 2){
		   my @aux_g = split("", $order[$i]);
		   if(($aux_g[0] eq "+") or ($aux_g[0] eq "-")){
                   	my $signal;
                    	if($aux_g[0] eq "+"){
                            $signal = "-";
                    	}
                    	else{
                            $signal = "+";
                    	}
                        my $len_g = scalar(@aux_g)-1;
                        my $aux_gene = join("", @aux_g[1..$len_g]);
                    	push @aux_order, $signal.$aux_gene;
			if($i > 0){
                    	    push @aux_order, $order[$i-1];
			}
                    }
                    else{
                    	push @aux_order, $order[$i];
			if($i > 0){
                    	    push @aux_order, $order[$i-1];
			}
                    }
                    my $t = $i - 1;
		}			
		if($add == 1){
                    my $t = $len - 1;
                    push @aux_order, $order[$t];
                }
	 	my $resp2 = verifySynteny(\@aux_order, \@coordinates);
		if($resp2 == 0){
		    $resp = $resp2;
		}
	    }
	    return $resp;
	}
    }
    else{
    	for(my $i = 0; $i < (scalar(@coordinates) - 1); ++$i){	   
            my @c1 = split("\t", $coordinates[$i]);
            my $end1 = $c1[1];

            my @c2 = split("\t", $coordinates[$i+1]);
            my $start2 = $c2[0];
	    my $int = $start2 - $end1;
            if($int <= $dist){
	    }
            else{
            	return 1;                
            }
      	}
    }
    return 0;
}

# This routine checks the intergenic distance

sub verifyCoordinatesWithoutSynteny{
    my $aux = shift;
    my @coordinates = @{$aux};
    for(my $i = 0; $i < (scalar(@coordinates) - 1); ++$i){
        my @c1 = split("\t", $coordinates[$i]);
        my $end1 = $c1[1];
        my @c2 = split("\t", $coordinates[$i+1]);
        my $start2 = $c2[0];
        my $int = $start2 - $end1;
        if($int <= $dist){
        }
        else{
            return 1;
        }
    }
    return 0;
}

# This routine checks the gene order 

sub verifyOrder{
    my $aux = shift;
    my @coordinates = @{$aux};
    my $count = 0;
    if($synteny){
	my @order = split(",", $synteny);
        if(lc($circular) eq "yes"){
            return verifyCircularOrder(\@order, \@coordinates);
        }
        else{
	    my $resp = verifySyntenyWithoutDistance(\@order, \@coordinates);
            if($resp != 0){
                my @aux_order = ();
                my $start;
                my $add = 0;
                my $len = scalar(@order);
                if(($len%2) == 0){
                    $start = $len - 2;
                    $add = 1;
                }
                else{
                    $start = $len - 1;
                }
                for(my $i = $start; $i >= 0; $i = $i - 2){
                   my @aux_g = split("", $order[$i]);
                   if(($aux_g[0] eq "+") or ($aux_g[0] eq "-")){
                        my $signal;
                        if($aux_g[0] eq "+"){
                            $signal = "-";
                        }
                        else{
                            $signal = "+";
                        }
                        my $len_g = scalar(@aux_g)-1;
                        my $aux_gene = join("", @aux_g[1..$len_g]);
                        push @aux_order, $signal.$aux_gene;
			if($i > 0){
                            push @aux_order, $order[$i-1];
			}
                    }
                    else{
                        push @aux_order, $order[$i];
			if($i > 0){
                            push @aux_order, $order[$i-1];
			}
                    }
		    my $t = $i - 1;
                }
		if($add == 1){
		    my $t = $len - 1;
		    push @aux_order, $order[$t];
		}
                my $resp2 = verifySyntenyWithoutDistance(\@aux_order, \@coordinates);
                if($resp2 == 0){
		    $resp = $resp2;
		}		
	    }
	    return $resp;
        }
    }
    return 0;
}

# This routine checks the synteny of the genes without considering the intergenic distance

sub verifySyntenyWithoutDistance{
    my $aux_order = shift;
    my $aux_genes = shift;
    my @order = @{$aux_order};
    my @genes = @{$aux_genes};
    my $size = scalar(@order);
    my $amount;
    if($size % 2 == 0){ 
        $amount = $size/2;
    }
    else{
        my $float = scalar(@order);
        $amount = int($float/2)+1; 
    }
    my $limit = scalar(@genes);
    if($amount == $limit){
        my $index = 0;
        for(my $i = 0; $i < ($limit-1); ++$i){
            my @c1 = split("\t", $genes[$i]);
            my @c2 = split("\t", $genes[$i+1]);
	    my $g1;
	    my $g2;
	    my @aux = split("", $order[$index]);
	    my $len = scalar(@aux)-1;
	    if(($aux[0] eq "+") or ($aux[0] eq "-")){
		$g1 = join("", @aux[1..$len]);
	    }
	    else{
		$g1 = $order[$index];
	    }
	    @aux = split("", $order[$index+2]);
            $len = scalar(@aux)-1;
            if(($aux[0] eq "+") or ($aux[0] eq "-")){
                $g2 = join("", @aux[1..$len]);
            }
            else{
                $g2 = $order[$index+2];
            }
            if(($c1[3] eq $g1) and ($c2[3] eq $g2)){}
            else{
                return 1; 
            }
            $index += 2;
        }
	return 0;       
    }
    else{
	my $last_index = -1;
        for(my $i = 0; $i < $limit; ++$i){
            my @c = split("\t", $genes[$i]);
            my $index;
            my $found = 0;
            for(my $j = 0; $j < scalar(@order); $j = $j + 2){
		my @aux = split("", $order[$j]);
            	my $len = scalar(@aux)-1;
		my $g;
            	if(($aux[0] eq "+") or ($aux[0] eq "-")){
                    $g = join("", @aux[1..$len]);
            	}
            	else{
                    $g = $order[$j];
            	}
           
                if($g eq $c[3]){
                    $index = $j;
                    $found = 1;
                }
            }
            if($found == 0){
                return 1;
            }
	    elsif($i == 0){
		$last_index = $index;
	    }
            else{
                if($index < $last_index){
                    return 1;
                }
            }
        }
        return 0;
    }

}

# This routine generates a protein sequence FASTA file for discarded results 

sub generateProteinFasta2{
    my $aux = shift;
    my $sequence = shift;
    my $dir = shift;
    my $key = shift;
    my @coordinates = @{$aux};
    my $start_element;
    for(my $i = 0; $i < scalar(@coordinates); ++$i){
        my @aux_coord = split("\t", $coordinates[$i]);
        my $start = $aux_coord[0];
        my $end = $aux_coord[1];
 	my $frame = $aux_coord[2];
        my $marker = $aux_coord[3];
        if(!defined $marker){
            next;
        }
        my $protein_file = $dir."/".$key."_".$marker."_protein.fasta";
	if($i == 0){
	    $start_element = $start;
	}
  	my $temp = $dir."/".$key."_".$marker."_orf.fasta";
	my $len = $end - $start + 1;
	my $t = length($sequence);
        my $orf = substr $sequence, $start, $len;
        open(TMP, ">$temp");
        print TMP ">$key\_$marker\n$orf\n";
        close(TMP);
        my $f;
        if($frame eq "-"){
            $f = -1;
        }
        else{
            $f = +1;
        }
        system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
	system "rm $temp";
	if(-e $protein_file){
	    my $aux_file = $dir."/".$key."_".$marker."_protein_aux.fasta";
            open(AUX, ">$aux_file");
	    open(FIL, $protein_file);
	    while(<FIL>){
	    	chomp($_);
	    	if($_ =~ />/){
		    print AUX $_."\n"; 
	    	}
	    	else{
		    my $seq = $_;
		    $seq =~ s/\*//g;
		    print AUX $seq."\n";
	    	}
	    }
	    close(FIL);
	    close(AUX);
	    system "mv $aux_file $protein_file";
	}
    }
}

# This routine generates a protein sequence FASTA file for positive results 

sub generateProteinFasta{
    my $aux = shift;
    my $sequence = shift;
    my $dir = shift;
    my $key = shift;
    my @coordinates = @{$aux};
    my @ori = split("", $sequence);
    my $ori_len = scalar(@ori);
    my @real_coord = ();
    for(my $i = 0; $i < scalar(@coordinates); ++$i){
        my @aux_coord = split("\t", $coordinates[$i]);
        my $start = $aux_coord[0];
        my $end = $aux_coord[1];
        my $frame = $aux_coord[2];
        my $marker = $aux_coord[3];
	if(!defined $marker){
	    next;
	}
        my $protein_file = $dir."/".$key."_".$marker."_protein.fasta";
	my $stop = 0;
        my $start_pos;
        my $end_pos;
        my $s;
        my $e;
        my $start_print;
        my $end_print;
        if($frame eq "+"){
            $s = $start-1;
            $e = $end-1;
            $start_print = $start;
            $end_print = $end;
        }
        else{
            $s = $end-1;
            $e = $start-1;
            $start_print = $end;
            $end_print = $start;
        }
        my $a;
        my $teste = $ori[$s].$ori[$s+1].$ori[$s+2];
        my $pos = $s - 3;
        while($stop == 0 and $pos > 0){
            my $codon = "";
            $codon = $ori[$pos].$ori[$pos+1].$ori[$pos+2];
            $a = $pos;
            if($frame eq "-"){
                my @aux = split("",$codon);
                my $new = "";
                for(my $j = 2; $j >= 0; --$j){
                    if(lc($aux[$j]) eq 'a'){
                        $new .= 't';
                    }
                    elsif(lc($aux[$j]) eq 't'){
                        $new .= 'a';
                    }
                    elsif(lc($aux[$j]) eq 'c'){
                        $new .= 'g';
                    }
                    elsif(lc($aux[$j]) eq 'g'){
                        $new .= 'c';
                    }
                }
                $codon = $new;
            }
	    if(uc($codon) =~ /$end_codon{$number_codes{$genetic_code}}/){	    
                 $stop = 1;
            }
            else{
                 $pos -= 3;
            }
        }
        if($stop == 0){
            $start_pos = 0;
        }
        else{
            $start_pos = $a;
        }
        $stop = 0;
        $pos = $e + 1;
        $teste = $ori[$e-2].$ori[$e-1].$ori[$e];
	while($stop == 0 and $pos < $ori_len){
            my $codon = "";
            $codon = $ori[$pos];
            if(($pos+1) < $ori_len){
                $codon .= $ori[$pos+1];
            }
            else{
                last;
            }
            if(($pos+2) < $ori_len){
                $codon .= $ori[$pos+2];
            }
            else{
                last;
            }
            $a = $pos + 2;
            if($frame eq "-"){
                my @aux = split("",$codon);
                my $new = "";
                for(my $j = 2; $j >= 0; --$j){
                    if(lc($aux[$j]) eq 'a'){
                        $new .= 't';
                    }
                    elsif(lc($aux[$j]) eq 't'){
                        $new .= 'a';
                    }
                    elsif(lc($aux[$j]) eq 'c'){
                        $new .= 'g';
                    }
                    elsif(lc($aux[$j]) eq 'g'){
                        $new .= 'c';
                    }
                }
                $codon = $new;
            }
	    if(uc($codon) =~ /$end_codon{$number_codes{$genetic_code}}/){
                $stop = 1;
            }
            else{
                $pos += 3;
            }
        }
        if($stop == 0){
            $end_pos = ($ori_len-1);
        }
        else{
            $end_pos = $a;
        }
        my $orf = "";
        for(my $k = $start_pos; $k <= $end_pos; ++$k){
            $orf .= $ori[$k];
        }
	my $add = 0;
        if($start_pos > 0){
	    if($frame eq "+"){
                my @aux_orf = split("", $orf);
                for(my $k = 0; $k < (scalar(@aux_orf)-2); $k = $k + 3){
                    my $codon = "";
                    $codon = $aux_orf[$k].$aux_orf[$k+1].$aux_orf[$k+2];
                    if(uc($codon) =~ /$start_codon{$number_codes{$genetic_code}}/){
                        $add = $k;
                        $orf = "";
                        for(my $j = $k; $j < scalar(@aux_orf); ++$j){
                            $orf .= $aux_orf[$j];
                        }
                        last;
                    }
                }
                $start_pos += $add;
            }
            else{
                my @aux_orf = split("", $orf);
                my $count_nt = 0;
                for(my $k = (scalar(@aux_orf)-1); $k >= 0; $k = $k - 3){
                    my $codon = "";
		    $codon = $aux_orf[$k-2].$aux_orf[$k-1].$aux_orf[$k];
                    my @aux = split("",$codon);
                    my $new = "";
                    for(my $j = 2; $j >= 0; --$j){
                        if(lc($aux[$j]) eq 'a'){
                            $new .= 't';
                        }
                        elsif(lc($aux[$j]) eq 't'){
                            $new .= 'a';
                        }
                        elsif(lc($aux[$j]) eq 'c'){
                            $new .= 'g';
                        }
                        elsif(lc($aux[$j]) eq 'g'){
                            $new .= 'c';
                        }
                    }
                    $codon = $new;
                    if(uc($codon) =~ /$start_codon{$number_codes{$genetic_code}}/){
                        $orf = "";
                        for(my $j = 0; $j <  (scalar(@aux_orf) - $count_nt); ++$j){
                            $orf .= $aux_orf[$j];
                        }
                        last;
                    }
                    $count_nt += 3;
                }
                $end_pos -= $count_nt;
            }
  	}
        my $str = $start_pos."\t".$end_pos."\t".$frame."\t".$marker;
        push @real_coord, $str;
        my $temp = $dir."/".$key."_".$marker."_orf.fasta";
        open(TMP, ">$temp");
        print TMP ">$key\_$marker\n$orf\n";
        close(TMP);
        my $f;
        if($frame eq "-"){
            $f = -1;
        }
        else{
            $f = +1;
        }
        system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
	system "rm $temp";
	my $aux_file = $dir."/".$key."_".$marker."_protein_aux.fasta";
        open(AUX, ">$aux_file");
        open(FIL, $protein_file);
        while(<FIL>){
            chomp($_);
            if($_ =~ />/){
                print AUX $_."\n";
            }
            else{
                my $seq = $_;
                $seq =~ s/\*//g;
                print AUX $seq."\n";
            }
        }
        close(FIL);
        close(AUX);
        system "mv $aux_file $protein_file";    
    }
}

# This routine generates FASTA sequence files for discarded sequences

sub generateFastas{
    my $aux = shift;
    my $sequence = shift;
    my $dir = shift;
    my $key = shift;
    my $index = shift;
    my $discarded = shift;
    my $org_id = shift;
    my $org_name = shift;
    my @coordinates = @{$aux};
    my $start_region = -1;
    my $end_region;
    my $end_partial = -1;
    my $element_coord;
    my %hash = ();
    my @elem_coords = ();
    my @elem_genes = ();
    my $original_index = $index;
    my $count_element = 0;
    my @aux_elements = ();
    for(my $i = 0; $i < scalar(@coordinates); ++$i){
        my @aux_coord = split("\t", $coordinates[$i]);
        my $start = $aux_coord[0];
        my $end = $aux_coord[1];
        my $frame = $aux_coord[2];
        my $marker = $aux_coord[3];
	if(defined $hash{$marker}){
	    my @aux_hash = ();
	    foreach my $m (sort keys %hash){
		push @aux_hash, $hash{$m};		
	    }
	    if(scalar(@aux_hash) > 1){
	    	my @aux1 = split("\t", $hash{$marker});
                my $start_old = $aux1[0];
                my $end_old = $aux1[1];
                @aux_hash = @{sortCoordinates(\@aux_hash)};
                my $resp_aux = verifyOrder(\@aux_hash);
		if($resp_aux == 0){ 
		    $resp_aux = verifyCoordinatesWithoutSynteny(\@aux_hash); 
		    if($resp_aux == 0){
			if(scalar(@aux_hash) == $total_of_models){
			    my @aux_coord = split("\t", $aux_hash[0]);
                            my $start_gene = $aux_coord[0];
                            @aux_coord = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                            my $end_gene = $aux_coord[1];
                            for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                my @aux_h = split("\t", $aux_hash[$g]);
                                my $start_h = $aux_h[0];
                                my $end_h = $aux_h[1];
                                my $frame_h = $aux_h[2];
                                my $marker_h = $aux_h[3];
                                my $protein_file = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                                my $temp = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                                my $len = $end_h - $start_h + 1;
                                my $orf = substr $sequence, $start_h, $len;
                                open(TMP, ">$temp");
                                print TMP ">$key\_$marker_h\_$index\n$orf\n";
                                close(TMP);
                                my $f;
                                if($frame_h eq "-"){
                                    $f = -1;
                                }
                                else{
                                    $f = +1;
                                }
                                system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
                                my $aux_name = ">$key\_$marker_h\_$index";
                                my $aux_file = $dir."/aux_file.fasta";
                                open(AUX, ">$aux_file");
                                open(FAS, $protein_file);
                                while(<FAS>){
                                    chomp($_);
                                    if($_ =~ />/){
                                        print AUX "$aux_name\n";
                                    }
                                    else{
					my $aux_seq = $_;
					$aux_seq =~ s/\*//g;
                                        print AUX "$aux_seq\n";
                                    }
                                }
                                close(AUX);
                                close(FAS);
                                system "mv $aux_file $protein_file";
                            }
			    if(($end_gene - $start_gene + 1) < $size_filter){
                                my $markers_str = generateStringMarkers(\@aux_hash);
                                print LOG "Sequence $key positive to $markers_str discarded - element size is lower than defined size filter ($size_filter)\n";
                                if(!-e $discarded){
                                    system "mkdir $discarded";
                                }
                                my $aux_file = "$dir/$key*\_$index\_protein.fasta";
                                system "mv $aux_file $discarded";
                                system "mv $dir/$key\_$index.fna $discarded";
                            }
                            else{
                                my @aux_start = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                                $end_partial = $aux_start[1];
                                @aux_start = split("\t", $aux_hash[0]);
                                $start_region = $aux_start[0];
				my $suffix;
				if($count_element > 0){
                                    if($count_element == 1){
					my $letter = $letters{1};
					my $prev_index = $index-1;
                                        my $aux_suffix = $prev_index.$letter;						
					for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                	    my @aux_h = split("\t", $aux_hash[$g]);
                                	    my $marker_h = $aux_h[3];
					    my $protein_old = $dir."/".$key."_".$marker_h."_".$prev_index."_protein.fasta";
                                    	    my $protein_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_protein.fasta";
                                    	    system "mv $protein_old $protein_new";
                                    	    my $orf_old = $dir."/".$key."_".$marker_h."_".$prev_index."_orf.fasta";
                                    	    my $orf_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_orf.fasta";
                                    	    system "mv $orf_old $orf_new";

					    my $aux_name = ">$key\_$marker_h\_$aux_suffix";
                                	    my $aux_file = $dir."/aux_file.fasta";
                                	    open(AUX, ">$aux_file");
                                	    open(FAS, $protein_new);
                                	    while(<FAS>){
                                    		chomp($_);
                                    		if($_ =~ />/){
                                        	    print AUX "$aux_name\n";
                                    		}
                                    		else{
                                        	    print AUX "$_\n";
                                    		}
                                	    }
                                	    close(AUX);
                                	    close(FAS);
                                	    system "mv $aux_file $protein_new";
					}
                                    }
				    my $letter = $letters{($count_element+1)};
				    $suffix = $original_index.$letter;
				    for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                        my @aux_h = split("\t", $aux_hash[$g]);
                                    	my $marker_h = $aux_h[3];
				    	my $protein_old = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
				    	my $protein_new = $dir."/".$key."_".$marker_h."_".$suffix."_protein.fasta";
				    	system "mv $protein_old $protein_new";
	                            	my $orf_old = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
				    	my $orf_new = $dir."/".$key."_".$marker_h."_".$suffix."_orf.fasta";
				    	system "mv $orf_old $orf_new";
                                            my $aux_name = ">$key\_$marker_h\_$suffix";
                                            my $aux_file = $dir."/aux_file.fasta";
                                            open(AUX, ">$aux_file");
                                            open(FAS, $protein_new);
                                            while(<FAS>){
                                                chomp($_);
                                                if($_ =~ />/){
                                                    print AUX "$aux_name\n";
                                                }
                                                else{
						    my $aux_seq = $_;
                                        	    $aux_seq =~ s/\*//g;
                                                    print AUX "$aux_seq\n";
                                                }
                                            }
                                            close(AUX);
                                            close(FAS);
                                            system "mv $aux_file $protein_new";
				    }				    
                                }
				else{
				    $suffix = $index;
				}
				my $str_coord = $start_region."\t".$end_partial;
                                push @elem_coords, $str_coord;
                                push @elem_genes, \@aux_hash;
				++$count_element;
                                ++$index;
                            }
			    %hash = ();
			    if(($start_old <= $start) and ($end_old >= $end) and ($start < $end)){ 
                            }
                            elsif(($start_old <= $start) and ($start < $end_old) and ($end >= $end_old)){
                            }
			    else{
				$hash{$marker} = $coordinates[$i];
			    }	
			}
			else{ 
			    if(($start_old <= $start) and ($end_old >= $end) and ($start < $end)){ 
                            }
                            elsif(($start_old <= $start) and ($start < $end_old) and ($end >= $end_old)){
                            }
                            else{
                                push @aux_elements, $coordinates[$i];
                            }
		        }			
		    }
		    else{ 
			$hash{$marker} = $coordinates[$i]; 
			@aux_hash = ();
            		foreach my $m (sort keys %hash){
                	    push @aux_hash, $hash{$m};
            		}
			$resp_aux = verifyCoordinatesWithoutSynteny(\@aux_hash); 
                    	if($resp_aux == 0){
			    if(scalar(@aux_hash) == $total_of_models){
                            	my @aux_coord = split("\t", $aux_hash[0]);
                            	my $start_gene = $aux_coord[0];
                            	@aux_coord = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                            	my $end_gene = $aux_coord[1];
                            	for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                    my @aux_h = split("\t", $aux_hash[$g]);
                                    my $start_h = $aux_h[0];
                                    my $end_h = $aux_h[1];
                                    my $frame_h = $aux_h[2];
                                    my $marker_h = $aux_h[3];
				    
                                    my $protein_file = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                                    my $temp = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                                    my $len = $end_h - $start_h + 1;
                                    my $orf = substr $sequence, $start_h, $len;
                                    open(TMP, ">$temp");
                                    print TMP ">$key\_$marker_h\_$index\n$orf\n";
                                    close(TMP);
                                    my $f;
                                    if($frame_h eq "-"){
                                    	$f = -1;
                                    }
                                    else{
                                    	$f = +1;
                                    }
                                    system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
                                    my $aux_name = ">$key\_$marker_h\_$index";
                                    my $aux_file = $dir."/aux_file.fasta";
                                    open(AUX, ">$aux_file");
                                    open(FAS, $protein_file);
                                    while(<FAS>){
                                    	chomp($_);
                                    	if($_ =~ />/){
                                            print AUX "$aux_name\n";
                                    	}
                                    	else{
					    my $aux_seq = $_;
                                            $aux_seq =~ s/\*//g;
                                            print AUX "$aux_seq\n";
                                    	}
                                    }
                                    close(AUX);
                                    close(FAS);
                                    system "mv $aux_file $protein_file";
                            	}	
				if(($end_gene - $start_gene + 1) < $size_filter){
                                    my $markers_str = generateStringMarkers(\@aux_hash);
                                    print LOG "Sequence $key positive to $markers_str discarded - element size is lower than defined size filter ($size_filter)\n";
                                    if(!-e $discarded){
                                    	system "mkdir $discarded";
                                    }
                                    my $aux_file = "$dir/$key*\_$index\_protein.fasta";
                                    system "mv $aux_file $discarded";
                                    system "mv $dir/$key\_$index.fna $discarded";
                            	}
                            	else{
                                    my @aux_start = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                                    $end_partial = $aux_start[1];
                                    @aux_start = split("\t", $aux_hash[0]);
                                    $start_region = $aux_start[0];
				    my $suffix;
				    if($count_element > 0){
                                    	if($count_element == 1){
                                            my $letter = $letters{1};
                                            my $prev_index = $index-1;
                                            my $aux_suffix = $prev_index.$letter;
					    for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                            	my @aux_h = split("\t", $aux_hash[$g]);
                                            	my $marker_h = $aux_h[3];
                                            	my $protein_old = $dir."/".$key."_".$marker_h."_".$prev_index."_protein.fasta";
                                            	my $protein_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_protein.fasta";
                                            	system "mv $protein_old $protein_new";
                                            	my $orf_old = $dir."/".$key."_".$marker_h."_".$prev_index."_orf.fasta";
                                            	my $orf_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_orf.fasta";
                                            	system "mv $orf_old $orf_new";
                                            	my $aux_name = ">$key\_$marker_h\_$aux_suffix";
                                            	my $aux_file = $dir."/aux_file.fasta";
                                            	open(AUX, ">$aux_file");
                                            	open(FAS, $protein_new);
                                            	while(<FAS>){
                                                    chomp($_);
                                                    if($_ =~ />/){
                                                    	print AUX "$aux_name\n";
                                                    }
                                                    else{
							my $aux_seq = $_;
                                            		$aux_seq =~ s/\*//g;
                                            		print AUX "$aux_seq\n";
                                                    }
                                            	}
                                            	close(AUX);
                                            	close(FAS);
                                            	system "mv $aux_file $protein_new";
					    }
                                    	}
                                    	my $letter = $letters{($count_element+1)};
                                    	$suffix = $original_index.$letter;
					 for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                            my @aux_h = split("\t", $aux_hash[$g]);
                                            my $marker_h = $aux_h[3];
                                    	    my $protein_old = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                                    	    my $protein_new = $dir."/".$key."_".$marker_h."_".$suffix."_protein.fasta";
                                    	    system "mv $protein_old $protein_new";
                                    	    my $orf_old = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                                    	    my $orf_new = $dir."/".$key."_".$marker_h."_".$suffix."_orf.fasta";
                                    	    system "mv $orf_old $orf_new";
                                            my $aux_name = ">$key\_$marker_h\_$suffix";
                                            my $aux_file = $dir."/aux_file.fasta";
                                            open(AUX, ">$aux_file");
                                            open(FAS, $protein_new);
                                            while(<FAS>){
                                                chomp($_);
                                                if($_ =~ />/){
                                                    print AUX "$aux_name\n";
                                                }
                                                else{
						    my $aux_seq = $_;
                                            	    $aux_seq =~ s/\*//g;
                                            	    print AUX "$aux_seq\n";
                                                }
                                            }
                                            close(AUX);
                                            close(FAS);
                                            system "mv $aux_file $protein_new";
					}
                                    }
                                    else{
                                        $suffix = $index;
                                   }
                                   my $str_coord = $start_region."\t".$end_partial;
                                   push @elem_coords, $str_coord;
                                   push @elem_genes, \@aux_hash;
                                   ++$index;
				   ++$count_element;
                                   %hash = ();
                            	}
			    }
		    	}
		    }
		}
		else{
		    $hash{$marker} = $coordinates[$i];
	   	}
	    }
        }
	else{
	    $hash{$marker} = $coordinates[$i];
	}
    }
    my @aux_hash = ();
    foreach my $m (sort keys %hash){
    	push @aux_hash, $hash{$m};
    }
	 
    if(scalar(@aux_hash) >= $minimum_amount){
        @aux_hash = @{sortCoordinates(\@aux_hash)};
        my $resp_aux = verifyOrder(\@aux_hash); 
        if($resp_aux == 0){
	    $resp_aux = verifyCoordinatesWithoutSynteny(\@aux_hash); 
	    if($resp_aux == 0){ 
            	my @aux_coord = split("\t", $aux_hash[0]);
            	my $start_gene = $aux_coord[0];
            	@aux_coord = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
            	my $end_gene = $aux_coord[1];
	    	for(my $g = 0; $g < scalar(@aux_hash); ++$g){
	    	    my @aux_h = split("\t", $aux_hash[$g]);
            	    my $start_h = $aux_h[0];
            	    my $end_h = $aux_h[1];
            	    my $frame_h = $aux_h[2];
            	    my $marker_h = $aux_h[3];
            	    my $protein_file = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
            	    my $temp = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
            	    my $len = $end_h - $start_h + 1;
            	    my $orf = substr $sequence, $start_h, $len;
            	    open(TMP, ">$temp");
            	    print TMP ">$key\_$marker_h\_$index\n$orf\n";
            	    close(TMP);
            	    my $f;
            	    if($frame_h eq "-"){
            	    	$f = -1;
            	    }
            	    else{
            	    	$f = +1;
            	    }
            	    system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
            	    my $aux_name = ">$key\_$marker_h\_$index";
            	    my $aux_file = $dir."/aux_file.fasta";
		    if(-e $protein_file){
            	    	open(AUX, ">$aux_file");
	    	    	open(FAS, $protein_file);
       	    	    	while(<FAS>){
            	    	    chomp($_);
            	    	    if($_ =~ />/){
                    	    	print AUX "$aux_name\n";
            	    	    }
           	    	    else{
				my $aux_seq = $_;
                                $aux_seq =~ s/\*//g;
                                print AUX "$aux_seq\n";
            	    	    }
            	    	}
            	    	close(AUX);
            	    	close(FAS);
            	    	system "mv $aux_file $protein_file >> /dev/null";
		    }
		    else{
			next;
		    }
	    	}
	    	if(($end_gene - $start_gene + 1) < $size_filter){
            	    my $markers_str = generateStringMarkers(\@aux_hash);
            	    print LOG "Sequence $key positive to $markers_str discarded - element size is lower than defined size filter ($size_filter)\n";
		    if(!-e $discarded){
                    	system "mkdir $discarded";
                    }
		    if(-e "$dir/$key\_$index.fna"){
			system "mv $dir/$key\_$index.fna $discarded >> /dev/null";
		    }
           	    my $aux_file = "$dir/$key*\_$index\_protein.fasta";
            	    system "mv $aux_file $discarded >> /dev/null";
            	}
            	else{
            	    my @aux_start = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
            	    $end_partial = $aux_start[1];
            	    @aux_start = split("\t", $aux_hash[0]);
            	    $start_region = $aux_start[0];
		    my $suffix;
		    if($count_element > 0){
                        if($count_element == 1){
                            my $letter = $letters{1};
                            my $prev_index = $index-1;
                            my $aux_suffix = $prev_index.$letter;
			    for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                my @aux_h = split("\t", $aux_hash[$g]);
                                my $marker_h = $aux_h[3];
                                my $protein_old = $dir."/".$key."_".$marker_h."_".$prev_index."_protein.fasta";
                                my $protein_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_protein.fasta";
                                system "mv $protein_old $protein_new";
                                my $orf_old = $dir."/".$key."_".$marker_h."_".$prev_index."_orf.fasta";
                                my $orf_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_orf.fasta";
                                system "mv $orf_old $orf_new";
                                my $aux_name = ">$key\_$marker_h\_$aux_suffix";
                                my $aux_file = $dir."/aux_file.fasta";
                                open(AUX, ">$aux_file");
                                open(FAS, $protein_new);
                                while(<FAS>){
                                    chomp($_);
                                    if($_ =~ />/){
                                    	print AUX "$aux_name\n";
                                    }
                                    else{
					my $aux_seq = $_;
                                        $aux_seq =~ s/\*//g;
                                        print AUX "$aux_seq\n";
                                   }
                               }
                               close(AUX);
                               close(FAS);
                               system "mv $aux_file $protein_new";
			    }
                        }
                        my $letter = $letters{($count_element+1)};
                        $suffix = $original_index.$letter;
			for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                            my @aux_h = split("\t", $aux_hash[$g]);
                            my $marker_h = $aux_h[3];
                            my $protein_old = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                            my $protein_new = $dir."/".$key."_".$marker_h."_".$suffix."_protein.fasta";
                            system "mv $protein_old $protein_new";
                            my $orf_old = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                            my $orf_new = $dir."/".$key."_".$marker_h."_".$suffix."_orf.fasta";
                            system "mv $orf_old $orf_new";
                            my $aux_name = ">$key\_$marker_h\_$suffix";
                            my $aux_file = $dir."/aux_file.fasta";
                            open(AUX, ">$aux_file");
                            open(FAS, $protein_new);
                            while(<FAS>){
                                chomp($_);
                                if($_ =~ />/){
                                    print AUX "$aux_name\n";
                                }
                                else{
				    my $aux_seq = $_;
                                    $aux_seq =~ s/\*//g;
                                    print AUX "$aux_seq\n";
                                }
                           }
                           close(AUX);
                           close(FAS);
                           system "mv $aux_file $protein_new";
			}
                   }
                   else{
                       $suffix = $index;
                   }
                   my $str_coord = $start_region."\t".$end_partial;
		   ++$count_element;
		   ++$index;
                   push @elem_coords, $str_coord;
                   push @elem_genes, \@aux_hash;
	    	}
            }
	    else{
		my $markers_str = generateStringMarkers(\@aux_hash);
		print LOG "Sequence $key positive to $markers_str discarded - distance between markers $markers_str is longer than maximum accepted distance.\n";
		if(!-e $discarded){
                    system "mkdir $discarded";
                }
		if(-e "$dir/$key\_$index.fna"){
		    system "mv $dir/$key\_$index.fna $discarded >> /dev/null";
		}
		my $aux_file = "$dir/$key*\_$index\_protein.fasta";
		if(-e $aux_file){
            	    system "mv $aux_file $discarded >> /dev/null";
		}
	    }
	}
        else{
	    my $markers_str = generateStringMarkers(\@aux_hash);
            print LOG "Sequence $key positive to $markers_str discarded - markers are in a wrong order.\n";
	    if(!-e $discarded){
                system "mkdir $discarded";
            }
            my $temp = $discarded."/".$key.".fasta";
            open(OUT, ">$temp");
            print OUT ">$key\n$sequence\n";
            close(OUT);
            generateProteinFasta2(\@aux_hash, $sequence, $discarded, $key);
	}
    }

    if(scalar(@aux_elements) >= $minimum_amount){
	my %hash_aux = ();
	my %hash2 = ();
	for(my $i = 0; $i < scalar(@aux_elements); ++$i){
	    my @aux_coord = split("\t", $aux_elements[$i]);
            my $start = $aux_coord[0];
            my $end = $aux_coord[1];
            my $frame = $aux_coord[2];
            my $marker = $aux_coord[3];
            if(defined $hash2{$marker}){
	    }
	    else{
		$hash2{$marker} = $aux_elements[$i];
	    }
	}
	my @aux_hash = ();
	foreach my $m (sort keys %hash2){
            push @aux_hash, $hash2{$m};
        }
	if(scalar(@aux_hash) >= $minimum_amount){
            @aux_hash = @{sortCoordinates(\@aux_hash)};
            my $resp_aux = verifyOrder(\@aux_hash);
	    if($resp_aux == 0){
            	$resp_aux = verifyCoordinatesWithoutSynteny(\@aux_hash); 
            	if($resp_aux == 0){
		    my @aux_coord = split("\t", $aux_hash[0]);
                    my $start_gene = $aux_coord[0];
                    @aux_coord = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                    my $end_gene = $aux_coord[1];
		    for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                    	my @aux_h = split("\t", $aux_hash[$g]);
                    	my $start_h = $aux_h[0];
                    	my $end_h = $aux_h[1];
                    	my $frame_h = $aux_h[2];
                    	my $marker_h = $aux_h[3];
                    	my $protein_file = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                    	my $temp = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                    	my $len = $end_h - $start_h + 1;
                    	my $orf = substr $sequence, $start_h, $len;
                    	open(TMP, ">$temp");
                    	print TMP ">$key\_$marker_h\_$index\n$orf\n";
                    	close(TMP);
                    	my $f;
                    	if($frame_h eq "-"){
                            $f = -1;
                    	}
                        else{
                            $f = +1;
                    	}
                    	system "transeq -frame $f -table $table $temp $protein_file 2> /dev/null";
                    	my $aux_name = ">$key\_$marker_h\_$index";
                    	my $aux_file = $dir."/aux_file.fasta";
                    	if(-e $protein_file){
                            open(AUX, ">$aux_file");
                            open(FAS, $protein_file);
                            while(<FAS>){
                            	chomp($_);
                            	if($_ =~ />/){
                                    print AUX "$aux_name\n";
                            	}
                            	else{
				    my $aux_seq = $_;
                                    $aux_seq =~ s/\*//g;
                                    print AUX "$aux_seq\n";
                            	}
                            }
                            close(AUX);
                            close(FAS);
                            system "mv $aux_file $protein_file >> /dev/null";
                    	}
                    	else{
                            next;
                    	}			
		    } 
		    if(($end_gene - $start_gene + 1) < $size_filter){
                    	my $markers_str = generateStringMarkers(\@aux_hash);
                    	print LOG "Sequence $key positive to $markers_str discarded - element size is lower than defined size filter ($size_filter)\n";
                    	if(!-e $discarded){
                            system "mkdir $discarded";
                    	}
                    	if(-e "$dir/$key\_$index.fna"){
                            system "mv $dir/$key\_$index.fna $discarded >> /dev/null";
                    	}
                    	my $aux_file = "$dir/$key*\_$index\_protein.fasta";
                    	system "mv $aux_file $discarded >> /dev/null";
                    }
                    else{
                    	my @aux_start = split("\t", $aux_hash[(scalar(@aux_hash)-1)]);
                    	$end_partial = $aux_start[1];
                    	@aux_start = split("\t", $aux_hash[0]);
                    	$start_region = $aux_start[0];
			my $suffix;
                    	if($count_element > 0){
                            if($count_element == 1){
                            	my $letter = $letters{1};
                            	my $prev_index = $index-1;
                            	my $aux_suffix = $prev_index.$letter;
                            	for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                                    my @aux_h = split("\t", $aux_hash[$g]);
                                    my $marker_h = $aux_h[3];
                                    my $protein_old = $dir."/".$key."_".$marker_h."_".$prev_index."_protein.fasta";
                                    my $protein_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_protein.fasta";
                                    system "mv $protein_old $protein_new";
                                    my $orf_old = $dir."/".$key."_".$marker_h."_".$prev_index."_orf.fasta";
                                    my $orf_new = $dir."/".$key."_".$marker_h."_".$aux_suffix."_orf.fasta";
                                    system "mv $orf_old $orf_new";
                            	}
                            }
                            my $letter = $letters{($count_element+1)};
                            $suffix = $original_index.$letter;
                            for(my $g = 0; $g < scalar(@aux_hash); ++$g){
                            	my @aux_h = split("\t", $aux_hash[$g]);
                            	my $marker_h = $aux_h[3];
                            	my $protein_old = $dir."/".$key."_".$marker_h."_".$index."_protein.fasta";
                            	my $protein_new = $dir."/".$key."_".$marker_h."_".$suffix."_protein.fasta";
                            	system "mv $protein_old $protein_new";
                            	my $orf_old = $dir."/".$key."_".$marker_h."_".$index."_orf.fasta";
                            	my $orf_new = $dir."/".$key."_".$marker_h."_".$suffix."_orf.fasta";
                            	system "mv $orf_old $orf_new";
                            }
                   	}
                   	else{
                           $suffix = $index;
                   	}
                        my $str_coord = $start_region."\t".$end_partial;
                    	push @elem_coords, $str_coord;
                    	push @elem_genes, \@aux_hash;
                    }
		} 
	    } 
	    else{
           	my $markers_str = generateStringMarkers(\@aux_hash);
           	print LOG "Sequence $key positive to $markers_str discarded - markers are in a wrong order.\n";
            	if(!-e $discarded){
                    system "mkdir $discarded";
            	}
            	my $temp = $discarded."/".$key.".fasta";
            	open(OUT, ">$temp");
            	print OUT ">$key\n$sequence\n";
            	close(OUT);
            	generateProteinFasta2(\@aux_hash, $sequence, $discarded, $key);
            }
	}
    }

    my $id = $original_index;
    my @aux = ();
    my $start = -1;
    my $end;
    my $str = "";
    my $use = 0;
    my $len_gen = length($sequence);
    my @resp_vector = ();    
    my $fna;
    my $t = scalar(@elem_coords);
    if(scalar(@elem_coords) == 0){} 
    elsif(scalar(@elem_coords) == 1){
	my @aux1 = split("\t", $elem_coords[0]);
	$start = $aux1[0];
	$end = $aux1[1];
        $start -= $flanking;
        $end += $flanking;
        if($start <= 0){
            $start = 1;
        }
        if($end > $len_gen){
            $end = $len_gen;
        }
 	my $len = $end - $start + 1;
        my $subst = substr $sequence, $start, $len;
        $fna = $dir."/".$key."_".$id.".fna";
        open(FNA, ">$fna");
        print FNA ">$key\_$id\n$subst\n";
        close(FNA);
        my $str_group = "";       
        my $basic_str="";
	my $delete = 0;
	if(scalar(@elem_genes) > 0){
	    if(defined $org_name){
		if($org_id != 1){
		    $basic_str =  "$org_id\t$org_name\t$key\_$id\t$len_gen\t$start-$end\t$len";
		}
		else{
		    $basic_str =  "$key\_$id\t$org_name\t$len_gen\t$start-$end\t$len";
		}
	    }
	    else{
		if(looks_like_number($org_id)){
		    if($org_id != 1){
		    	$basic_str =  "$org_id\t$key\_$id\t$len_gen\t$start-$end\t$len";
		    }
		    else{
			$basic_str =  "$key\_$id\t$len_gen\t$start-$end\t$len";
		    }
		}
		else{
            	    $basic_str =  "$org_id\t$key\_$id\t$len_gen\t$start-$end\t$len";
		}
	    }
            for(my $j = 0; $j < scalar(@elem_genes); ++$j){		
		$use = 0;
                my @aux_g = @{$elem_genes[$j]};		
                my $markers_str = generateStringMarkers(\@aux_g);
		if($synteny){
           	    my $resp =  verifyCoordinates(\@aux_g);
           	    if($resp == 1 || $resp == 4 || $resp == 3){
			if($resp == 1){
                	    print LOG "Sequence $key positive to $markers_str discarded - distance between markers $markers_str is longer than maximum accepted distance.\n";
			}
		   	elsif($resp == 3){
			    print LOG "Sequence $key positive to $markers_str discarded - markers are in a wrong order.\n";	
			}
			elsif($resp == 4){
			    print LOG "Sequence $key positive to $markers_str discarded - markers  $markers_str are not in the right strand.\n";
			}
                    	if(!-e $discarded){
                            system "mkdir $discarded";
                    	}
			if(-e "$dir/$key\_$id.fna"){
                            system "mv $dir/$key\_$id.fna $discarded >> /dev/null";
                        }
			my $num = `ls $dir/$key*\_$id\_protein.fasta | wc -l`;
			if($num > 0){
                   	    system "mv  $dir/$key*\_$id\_protein.fasta $discarded >> /dev/null";
			    system "rm $dir/$key*\_$id\_orf.fasta >> /dev/null";
			}
			next;
                    }
        	}
		my $t = scalar(@aux_g);
                for(my $k = 0; $k < scalar(@aux_g); ++$k){
		    if(!-e $fna){
			last;
		    }
                    my @g1 = split("\t", $aux_g[$k]);
                    my $cont;
                    if($g1[0] == 0){
                    	my $st = $g1[0] + 1;
                        my $en = $g1[1] + 1;
                        $cont = $st."-".$en;
                    }
                    else{
                        $cont = $g1[0]."-".$g1[1];
                    }
                    my $element;
                    if($len == $len_gen){
                    	$element = $cont;
                    }
                    else{
                        $element = extractCoordinatesFromFasta($dir, $key, $g1[3], $id, $id, 1);
			if(!defined $element){
			   last;
			}
                    }
                    if($k == (scalar(@aux_g) - 1)){
                        $str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]";
                    }
                    else{
                        my @g2 = split("\t", $aux_g[$k+1]);
                        my $distance = $g2[0] - $g1[1];
                        if($distance < 0 and abs($distance) > $overlap){
                            $use = 1;
                            last;
                        }
                        $str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]\t$distance";
                    }
                }
                if($use == 0){
                    print PARTIAL $basic_str.$str_group."\n";
                    print LOG "Sequence $key positive to $markers_str valid.\n";
		    ++$count_element;
                    $str .= $basic_str.$str_group."\n";
		    $delete = 1;
                }
                else{
		    if(!-e $discarded){
                    	system "mkdir $discarded";
                    }
		    system "mv $dir/$key*\_$id\_protein.fasta $discarded >> /dev/null";
		    system "rm $dir/$key*\_$id\_orf.fasta >> /dev/null";
                    print LOG "Sequence $key positive to $markers_str discarded - overlap between markers is longer than maximum accepted overlap ($overlap).\n";
                 }
                 $str_group = "";
             }
	}
	if($delete == 0){
	    if(!-e $discarded){
            	system "mkdir $discarded";
            }
            if(-e "$fna"){
                system "mv $dir/$key\_$id.fna $discarded >> /dev/null";
            }
	}
    }
    else{
	my @grouped_genes = ();
	my $basic_str;
	my $str_group = "";
	my $delete = 0;
    	for(my $i = 0; $i < scalar(@elem_coords) - 1; ++$i){
	    my @aux1 = split("\t", $elem_coords[$i]);
	    my @aux2 = split("\t", $elem_coords[$i+1]);    	    
	    push @grouped_genes, $elem_genes[$i];
	    if($start == -1){
	    	$start = $aux1[0];
	    }
	    my $d = $aux2[0] - $aux1[1];
	    my $count_elements = 1;
	    if($d > $element_distance){ 
	    	$end = $aux1[1];
	    	$start -= $flanking;
	    	$end += $flanking;
	    	if($start <= 0){
        	    $start = 1;
    	    	}
    	    	if($end > $len_gen){
        	    $end = $len_gen;
    	    	}
	    	my $len = $end - $start + 1;
   	    	my $subst = substr $sequence, $start, $len;
	    	$fna = $dir."/".$key."_".$id.".fna";
    	    	open(FNA, ">$fna");
    	    	print FNA ">$key\_$id\n$subst\n";
    	    	close(FNA);	    	    
	    	$str_group = "";
		my $replace_basic = 0;
	    	if(scalar(@grouped_genes) > 0){
		    if(scalar(@grouped_genes) == 1){
			if(defined $org_name){
                    	    if($org_id != 1){
                        	$basic_str =  "$org_id\t$org_name\t$key\_$d\t$len_gen\t$start-$end\t$len";
                    	    }
                    	    else{
                        	$basic_str =  "$key\_$id\t$org_name\t$len_gen\t$start-$end\t$len";
                    	    }
                	}
                	else{
                    	    if($org_id != 1){
                        	$basic_str =  "$org_id\t$key\_$id\t$len_gen\t$start-$end\t$len";
                    	    }
                    	    else{
                       		$basic_str =  "$key\_$id\t$len_gen\t$start-$end\t$len";
                    	    }
                	}			
		    }
		    else{
			$replace_basic = 1;
		    }
		    for(my $j = 0; $j < scalar(@grouped_genes); ++$j){
		    	my @aux_g = @{$grouped_genes[$j]};			
		    	my $markers_str = generateStringMarkers(\@aux_g);
			my $aux_index = $j+1;
			my $letter = $letters{$aux_index};
			my $suffix_id;
			if(defined $letter){
      	                    $suffix_id = $id.$letter;
		 	    if(defined $org_name){
                                if($org_id != 1){
                                    $basic_str =  "$org_id\t$org_name\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                                }
                                else{
                                    $basic_str =  "$key\_$suffix_id\t$org_name\t$len_gen\t$start-$end\t$len";
                                }
                            }
                            else{
                                if($org_id != 1){
                                    $basic_str =  "$org_id\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                                }
                                else{
                                    $basic_str =  "$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                                }
                            }
			}
		        else{
			    $suffix_id = $id;
	 		}
			if($replace_basic == 1){
			    $basic_str = "";
                	    my $aux_index = $j+1;
                	    my $letter = $letters{$aux_index};
                	    my $suffix_id = $id.$letter;
                	    if(defined $org_name){
                    	    	if($org_id != 1){
                            	    $basic_str =  "$org_id\t$org_name\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	    	}
                    	    	else{
                        	    $basic_str =  "$key\_$suffix_id\t$org_name\t$len_gen\t$start-$end\t$len";
                    	    	}
                	    }
                	    else{
                    	    	if($org_id != 1){
                            	    $basic_str =  "$org_id\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	    	}
                    	        else{
                        	    $basic_str =  "$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	    	}
                	    }
			}
			if($synteny){
                    	    my $resp =  verifyCoordinates(\@aux_g);
                            if($resp == 1 || $resp == 4 || $resp == 3){
				if($resp == 1){
                            	    print LOG "Sequence $key positive to $markers_str discarded - distance between markers $markers_str is longer than maximum accepted distance.\n";
				}
				elsif($resp == 3){
                            	    print LOG "Sequence $key positive to $markers_str discarded - markers are in a wrong order.\n";
                        	}
				elsif($resp == 4){
                                     print LOG "Sequence $key positive to $markers_str discarded - markers  $markers_str are not in a same strand.\n";
                        	}

				if(!-e $discarded){
                                    system "mkdir $discarded";
                                }
				if(-e "$dir/$key\_$id.fna"){
                            	    system "mv $dir/$key\_$id.fna $discarded >> /dev/null";
                        	}
				my $num = `ls $dir/$key*\_$suffix_id\_protein.fasta | wc -l `;
				if($num > 0){
                            	    system "mv  $dir/$key*\_$suffix_id\_protein.fasta $discarded >> /dev/null";
                                    system "rm $dir/$key*\_$suffix_id\_orf.fasta >> /dev/null";
				}
                            	next;
                    	    }
                        }  
			my $t = scalar(@aux_g);
		    	for(my $k = 0; $k < scalar(@aux_g); ++$k){			
		    	    my @g1 = split("\t", $aux_g[$k]);
                    	    my $cont;
                    	    if($g1[0] == 0){
                    	   	my $st = $g1[0] + 1;
                    	   	my $en = $g1[1] + 1;
                    	   	$cont = $st."-".$en;
                    	    }
                    	    else{
			    	$cont = $g1[0]."-".$g1[1];
                    	    }
			    my $element;			   
                	    if($len == $len_gen){
                    	    	$element = $cont;
                	    }
                	    else{
                    	    	$element = extractCoordinatesFromFasta($dir, $key, $g1[3], $suffix_id, $id, ($j+1));
				if(!defined $element){
				    last;
				}
                	    }
			    if($k == (scalar(@aux_g) - 1)){
                    	    	$str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]";
                	    }
                	    else{
                    	    	my @g2 = split("\t", $aux_g[$k+1]);
                    	    	my $distance = $g2[0] - $g1[1];
                    	    	if($distance < 0 and abs($distance) > $overlap){
                        	    $use = 1;
                        	    last;
                    	    	}
                    	    	$str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]\t$distance";
                	    }
		    	}
		    	if($use == 0){
			    if($str_group ne ""){
                    	        print PARTIAL $basic_str.$str_group."\n";
                    	        print LOG "Sequence $key positive to $markers_str valid.\n";
			        ++$count_element;
                    	        ++$id;
                    	        $str .= $basic_str.$str_group."\n";
			        push @resp_vector, 0;
			        $delete = 1;
			    }
                    	}
                    	else{
                    	    if(!-e $discarded){
                            	system "mkdir $discarded";
                    	    }
			    system "mv $dir/$key*\_$suffix_id\_protein.fasta $discarded >> /dev/null";
			    system "rm $dir/$key*\_$suffix_id\_orf.fasta >> /dev/null";
                    	    print LOG "Sequence $key positive to $markers_str discarded - overlap between markers is longer than maximum accepted overlap ($overlap).\n";
			    push @resp_vector, 1;
                    	}
                    	$str_group = "";
		    }
	    	}
		@grouped_genes = ();
		push @grouped_genes, $elem_genes[$i+1];
	    	$start = -1;
		if($delete == 0){
		    if(!-e $discarded){
                        system "mkdir $discarded";
                    }
                    if(-e "$fna"){
                        system "mv $fna $discarded >> /dev/null";
                    }
		}
	    }
	    else{ 
	    	$end = $aux2[1];
		push @grouped_genes, $elem_genes[$i+1];
		++$count_elements;
	    }
    	}
  	my $t = scalar(@grouped_genes);
	if(scalar(@grouped_genes) > 0){	   
	    $start -= $flanking;
            $end += $flanking;
            if($start <= 0){
            	$start = 1;
            }
            if($end > $len_gen){
            	$end = $len_gen;
            }
            my $len = $end - $start + 1;
            my $subst = substr $sequence, $start, $len;
            my $fna = $dir."/".$key."_".$id.".fna";
            open(FNA, ">$fna");
            print FNA ">$key\_$id\n$subst\n";
            close(FNA);
	    my $replace_basic = 0;
	    if(scalar(@grouped_genes) == 1){
            	if(defined $org_name){
                    if($org_id != 1){
                    	$basic_str =  "$org_id\t$org_name\t$key\_$id\t$len_gen\t$start-$end\t$len";
                    }
                    else{
                        $basic_str =  "$key\_$id\t$org_name\t$len_gen\t$start-$end\t$len";
                    }
                }
                else{
                    if($org_id != 1){
                        $basic_str =  "$org_id\t$key\_$id\t$len_gen\t$start-$end\t$len";
                    }
                    else{
                        $basic_str =  "$key\_$id\t$len_gen\t$start-$end\t$len";
                    }
                }
            }
            else{
                $replace_basic = 1;
            }
            $delete = 0;
            for(my $j = 0; $j < scalar(@grouped_genes); ++$j){
                my @aux_g = @{$grouped_genes[$j]};
                my $markers_str = generateStringMarkers(\@aux_g);
		my $suffix_id = $id;
		if($replace_basic == 1){
		    $basic_str = "";
                    my $aux_index = $j+1;
                    my $letter = $letters{$aux_index};
                    $suffix_id = $id.$letter;
                    if(defined $org_name){
                    	if($org_id != 1){
                            $basic_str =  "$org_id\t$org_name\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	}
                    	else{
                            $basic_str =  "$key\_$suffix_id\t$org_name\t$len_gen\t$start-$end\t$len";
                    	}
                    }
                    else{
                    	if($org_id != 1){
                            $basic_str =  "$org_id\t$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	}
                    	else{
                            $basic_str =  "$key\_$suffix_id\t$len_gen\t$start-$end\t$len";
                    	}
                    }
		}
		if($synteny){
                    my $resp = verifyCoordinates(\@aux_g);
                    if($resp == 1 || $resp == 4 || $resp == 3){
			if($resp == 1){
                            print LOG "Sequence $key positive to $markers_str discarded - distance between markers $markers_str is longer than maximum accepted distance.\n";
			}
			elsif($resp == 3){
                            print LOG "Sequence $key positive to $markers_str discarded - markers are in a wrong order.\n";
                        }
			elsif($resp == 4){
                            print LOG "Sequence $key positive to $markers_str discarded - markers  $markers_str are not in a same strand.\n";
                        }

			if(!-e $discarded){
                            system "mkdir $discarded";
                        }
                        if(-e "$dir/$key\_$id.fna"){
			    system "mv $dir/$key\_$id.fna $discarded >> /dev/null";
			}
                        my $aux_file = "$dir/$key*\_$suffix_id\_protein.fasta";
			my $num = `ls $aux_file | wc -l`;
			if($num > 0){
                            system "mv $aux_file $discarded >> /dev/null";
			    system "rm $dir/$key*\_$suffix_id\_orf.fasta >> /dev/null";
			}
                        next;
                     }
                }
		my $t = scalar(@aux_g);
		$str_group = "";
                for(my $k = 0; $k < $t; ++$k){
                    my @g1 = split("\t", $aux_g[$k]);
                    my $cont;
                    if($g1[0] == 0){
                        my $st = $g1[0] + 1;
                        my $en = $g1[1] + 1;
                        $cont = $st."-".$en;
                    }
                    else{
                        $cont = $g1[0]."-".$g1[1];
                    }
                    my $element;
                    if($len == $len_gen){
                    	$element = $cont;
                    }
                    else{
                        $element = extractCoordinatesFromFasta($dir, $key, $g1[3], $suffix_id, $id, ($j+1));		
			if(!defined $element){
			    last;
			}
                    }
                    if($k == (scalar(@aux_g) - 1)){
                        $str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]";
                    }
                    else{
                        my @g2 = split("\t", $aux_g[$k+1]);
                        my $distance = $g2[0] - $g1[1];
                        if($distance < 0 and abs($distance) > $overlap){
                            $use = 1;
                            last;
                        }
                        $str_group .= "\t$g1[3]\t$cont\t$element\t$g1[2]\t$distance";
                    }
             	}
                if($use == 0){
		    if($str_group ne ""){
                    	print PARTIAL $basic_str.$str_group."\n";
                    	print LOG "Sequence $key positive to $markers_str valid.\n";
		    	++$count_element;
                    	$str .= $basic_str.$str_group."\n";
                    	push @resp_vector, 0;
		    	$delete = 1;
		    }
                }
                else{
		    if(!-e $discarded){
                    	system "mkdir $discarded";
                    }
		    my $num = `ls $dir/$key*\_$suffix_id\_protein.fasta | wc -l`;
		    if($num > 0){
		    	system "mv $dir/$key*\_$suffix_id\_protein.fasta $discarded >> /dev/null";
			system "rm $dir/$key*\_$suffix_id\_orf.fasta >> /dev/null";
		    }
                    print LOG "Sequence $key positive to $markers_str discarded - overlap between markers is longer than maximum accepted overlap ($overlap).\n";
                    push @resp_vector, 1;
                }
		$str_group = "";
      	    }
	    if($delete == 0){
		if(-e $fna){
		    system "mv $fna $discarded >> /dev/null";
		}
	    }
    	}
    }
    my $resp = 1;
    for(my $i = 0; $i < scalar(@resp_vector); ++$i){
	if($resp_vector[$i] == 0){
	    $resp = 0;
	    last;
  	}
    }
    if($count_element == 0){
	--$id;
    }
    return ($str, $resp, $id);
}

# This routine returns the CDS coordinates in the genome

sub extractCoordinatesFromFasta{
    my $dir = shift;
    my $key = shift;
    my $marker = shift;
    my $rand = shift;
    my $rand_fna = shift;
    my $count = shift;
    my $orf = $dir."/".$key."_".$marker."_".$rand."_orf.fasta";
    my $fna = $dir."/".$key."_".$rand_fna.".fna";
    my $blast = $dir."/".$key."_".$marker."_".$rand_fna."_blast_temp.txt";    
    if(!-e $fna or !-e $orf){
	return undef;
    }
    if(!-e $blast){
        runBlastn($orf, $fna, $blast);
    }
    if(-e $blast){      
    	my $string_coord = extractCoordinatesToFastaFile($blast, $key, $count);
    	my @aux = split("\t", $string_coord);
    	system "rm $blast";  
    	if($aux[0] < $aux[1]){
    	    return $aux[0]."-".$aux[1];
    	}
   	else{
	    return $aux[1]."-".$aux[0];
    	}
    }
    return undef;   
}

# This routine finds the ORF coordinates of a similarity search result.
sub findORF{
    my $start = shift;
    my $end = shift;
    my $frame = shift;
    my $sequence = shift;
    my @ori = split("", $sequence);
    my $ori_len = scalar(@ori);
    my $s;
    my $e;
    my $start_final;
    my $end_final;
    if($frame eq "+"){
    	$s = $start-1;
        $e = $end-1;
    }
    else{
        $s = $end-1;
        $e = $start-1;
    }
    my $a;
    my $teste = $ori[$s].$ori[$s+1].$ori[$s+2];
    my $pos = $s - 3; 
    my $stop = 0;
    while($stop == 0 and $pos > 0){
    	my $codon = "";
        $codon = $ori[$pos].$ori[$pos+1].$ori[$pos+2];
        $a = $pos;
        if($frame eq "-"){
            my @aux = split("",$codon);
            my $new = "";
            for(my $j = 2; $j >= 0; --$j){
            	if(lc($aux[$j]) eq 'a'){
                    $new .= 't';
                }
                elsif(lc($aux[$j]) eq 't'){
                    $new .= 'a';
                }
                elsif(lc($aux[$j]) eq 'c'){
                    $new .= 'g';
                }
                elsif(lc($aux[$j]) eq 'g'){
                    $new .= 'c';
                }
            }
            $codon = $new;
        }
        if(uc($codon) =~ /$end_codon{$number_codes{$genetic_code}}/){
            $stop = 1;
        }
        else{
	    if(($pos - 3) > 0){
               $pos -= 3;
	    }
	    else{
		last;
	    }
        }
    }
    if($pos > 0){
	$start_final = $pos;
    } 
    else{
	$start_final = 0;
    }
    $stop = 0;
    $pos = $e + 1;
    $teste = $ori[$e-2].$ori[$e-1].$ori[$e];
    while($stop == 0 and $pos < $ori_len){
    	my $codon = "";
        $codon = $ori[$pos];
        if(($pos+1) < $ori_len){
            $codon .= $ori[$pos+1];
        }
        else{
            last;
        }
        if(($pos+2) < $ori_len){
            $codon .= $ori[$pos+2];
        }
        else{
            last;
        }
        $a = $pos + 2;
        if($frame eq "-"){
            my @aux = split("",$codon);
            my $new = "";
            for(my $j = 2; $j >= 0; --$j){
            	if(lc($aux[$j]) eq 'a'){
                    $new .= 't';
                }
                elsif(lc($aux[$j]) eq 't'){
                    $new .= 'a';
                }
                elsif(lc($aux[$j]) eq 'c'){
                    $new .= 'g';
                }
                elsif(lc($aux[$j]) eq 'g'){
                    $new .= 'c';
                }
            }
            $codon = $new;
        }
	if(uc($codon) =~ /$end_codon{$number_codes{$genetic_code}}/){
            $stop = 1;
        }
        else{
            $pos += 3;
        }
    }
    if($stop == 0){
        $end_final = ($ori_len-1);
    }
    else{
        $end_final = $a;
    }
    my $orf = "";
    for(my $k = $start_final; $k <= $end_final; ++$k){
    	$orf .= $ori[$k];
    }
    my $add = 0;
    if($start_final > 1){
    	if($frame eq "+"){
            my @aux_orf = split("", $orf);
            for(my $k = 0; $k < (scalar(@aux_orf)-2); $k = $k + 3){
            	my $codon = "";
            	$codon = $aux_orf[$k].$aux_orf[$k+1].$aux_orf[$k+2];
            	if(uc($codon) =~ /$start_codon{$number_codes{$genetic_code}}/){		    
            	    $add = $k;
		    last;
            	}
            }
	    if($add == 0){
		return undef;
	    }
	    $start_final += $add;    
      	}
	else{
	    my @aux_orf = split("", $orf);
            my $count_nt = 0;
            for(my $k = (scalar(@aux_orf)-1); $k >= 0; $k = $k - 3){
            	my $codon = "";
                $codon = $aux_orf[$k-2].$aux_orf[$k-1].$aux_orf[$k];
                my @aux = split("",$codon);
                my $new = "";
                for(my $j = 2; $j >= 0; --$j){
                    if(lc($aux[$j]) eq 'a'){
                    	$new .= 't';
                    }
                    elsif(lc($aux[$j]) eq 't'){
                        $new .= 'a';
                    }
                    elsif(lc($aux[$j]) eq 'c'){
                        $new .= 'g';
                    }
                    elsif(lc($aux[$j]) eq 'g'){
                        $new .= 'c';
                    }
                }
                $codon = $new;
                if(uc($codon) =~ /$start_codon{$number_codes{$genetic_code}}/){
                    last;
                }
                $count_nt += 3;
             }
             $end_final -= $count_nt;
	}
    }    
    return "$start_final\t$end_final\t$frame";
}

# This routine analyzes a group of genes and defines if they belong to the same element

sub analiseGroupOfGenes{
    my $aux = shift;
    my $discarded = shift;
    my $key = shift;
    my $sequence = shift;
    my $dir_final = shift;
    my $org_id = shift;
    my $org_name = shift;
    my $index = shift;
    my $previous_elem = shift;
    my @coord = @{$aux};
    my $markers_str = "";
    my $resp = 0;
    my $str_report = "";
    my $current_element;
    my $index_return = $index;
    if(scalar(@coord) < $minimum_amount){
	if(!-e $discarded){
            system "mkdir $discarded";
        }
	my $markers_str = generateStringMarkers(\@coord);
	print LOG "Sequence $key positive to $markers_str discarded - number of markers is lower than $minimum_amount\n";
        my $temp = $discarded."/".$key.".fasta";
        open(OUT, ">$temp");
        print OUT ">$key\n$sequence\n";
        close(OUT);
	
	generateProteinFasta2(\@coord, $sequence, $discarded, $key);
	$resp = 1;	
    }
    else{
	my @aux_coord = @{removeEqualsCoordinates(\@coord)};
  	@aux_coord = @{sortCoordinates(\@aux_coord)};	
	my ($str, $resp, $id) = generateFastas(\@aux_coord, $sequence, $dir_final, $key, $index, $discarded, $org_id, $org_name);
	$str_report .= $str;
	$index_return = $id;
    }     
    return ($resp, $str_report, $index_return);
}

# This routine generates a string containing a list of all genes found in a sequence

sub generateStringMarkers{
    my $aux = shift;
    my @coord = @{$aux};
    my $markers_str = "";
    for(my $i = 0; $i < scalar(@coord); ++$i){
        my @aux_m = split("\t", $coord[$i]);
        if($i == scalar(@coord)-2){
            $markers_str .= $aux_m[3]." and ";

       }
       elsif($i == scalar(@coord)-1){
           $markers_str .= $aux_m[3];
       }
       else{
            $markers_str .= $aux_m[3].", ";
        }
    }
    return $markers_str;
}

# This routine verifies the type of input file

sub verify_file_type{ # 1 - fasta file; 2 - fastq file; 3 - tabular file
    my $file = shift;
    open(FILE, $file);
    my $found = 0;
    my $count = 0;
    while(<FILE>){
        chomp($_);
        if($_ =~ /^>/){
            close(FILE);
            return 1;
        }
        if($_ =~ /^@/){
            $found = 1;
        }
        elsif(($_ =~ /^\+/) and ($found == 1)){
            close(FILE);
            return 2;
        }
        ++$count;
        if($count > 100){
            last;
        }
    }
    close(FILE);
    return 3;
}

# This routine discards redundant genes

sub discard_redundant_genes{
    my $aux = shift;
    my @coordinates = @{$aux};
    my %hash = ();
    for(my $i = 0; $i < scalar(@coordinates); ++$i){
        my @aux_coord = split("\t", $coordinates[$i]);
        my $start = $aux_coord[0];
        my $end = $aux_coord[1];
        my $frame = $aux_coord[2];
        my $marker = $aux_coord[3];    
	push @{$hash{$marker}}, $coordinates[$i];
    }
    my @final = ();
    foreach my $key (sort keys %hash){
	my @aux = @{$hash{$key}};
	@aux = @{sortCoordinates(\@aux)};
        @aux = @{discarded_equal_genes(\@aux)};
	my $size = scalar(@aux);
	if($size > 1){
	    @aux = @{sortCoordinates(\@aux)};
	    @aux = @{discarded_equal_genes(\@aux)};
	    my $current = $aux[0];
	    my @aux_current = split("\t", $current);
	    my $start = $aux_current[0];
	    my $end = $aux_current[1];
	    for(my $i = 1; $i < $size; ++$i){
		my $next = $aux[$i];
		my @aux_next = split("\t",$aux[$i]);
		my $start_n = $aux_next[0];
		my $end_n = $aux_next[1];
		if($start == $start_n and $end == $end_n){
		}
		elsif($start_n > $start and $end_n > $end and $start_n > $end){
		    push @final, $current;
                    @aux_current = split("\t",$aux[$i]);;
                    $current = $aux[$i];
		}
		elsif($start_n >= $start and $end_n <= $end){
		    my $size_p = abs($start - $end);
		    my $size_n = abs($start_n - $end_n);
		    if($size_p >= $size_n){
		    }
		    else{
			@aux_current = split("\t",$aux[$i]);
			$current = $aux[$i];
		    }
		}
	    }
	    push @final, $current;
	}
	else{
	    push (@final, $aux[0]);
	}
    }
    @final = @{sortCoordinates(\@final)};
    return \@final;
}

sub discarded_equal_genes{
    my $aux = shift;
    my @coordinates = @{$aux};
    my @final = ();
    my $c_start;
    my $c_end;
    for(my $i = 0; $i < scalar(@coordinates); ++$i){
        my @aux_coord = split("\t", $coordinates[$i]);
        my $start = $aux_coord[0];
        my $end = $aux_coord[1];
        my $frame = $aux_coord[2];
        my $marker = $aux_coord[3];
	if($i == 0){
	    $c_start = $start;
	    $c_end = $end;
	    push @final, $coordinates[$i];
	}
	elsif($c_start == $start and $c_end == $end){

	}
	else{
	    $c_start = $start;
            $c_end = $start;
            push @final, $coordinates[$i];	
	}
    }
    return \@final;
}
