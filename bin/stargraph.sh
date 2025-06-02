#!/bin/bash
set -euo pipefail

version="v1"

#####################################################################################
############# STEP -1. CREATING THE ENVIRONMENT. NOT ACTUALLY USED NOW ##############
#####################################################################################

### need to create a conda environment for housing all the stargraph dependencies:
## TRIED INSTALLING STARFISH IN THE SAME ENVIRONMENT BUT THE DEPENDENCIES AND NAMING MAKE IT IMPOSSIBLE ESSENTIALLY...
## therefore the starfish portion has to be run outside of stargraph and only input from that pipeline will be accepted

##starfish (egluckthaler::starfish) NOT USED (trying to install it directly with the stargraph conda env)
##pggb/odgi (bioconda::pggb)
##bedtools (bioconda::bedtools)
##mash	(bioconda::mash)
##seqkit (bioconda::seqkit)
##minimap2 (biocojnda::minimap2)
##mummer4 (bioconda::mummer4)
##R (conda-forge::r-base) NOT USED NOW
##gggenomes (conda-forge::r-gggenomes) NOT USED NOW


##starfish conda packages
##samtools (upgrading) samtools=1.6 (should not impact other tools)
##sourmash bioconda::sourmash
##blast bioconda::blast
##mcl bioconda::mcl=14.137 (MOST RECENT VERSION IS MUCH MORE RECENT, THEREFORE SPECIFYING)
##circos bioconda::circos
##hmmer bioconda::hmmer
##metaeuk bioconda::metaeuk (A DIFFERENT VERSION THAN SPECIFIED, SHOULDN'T CHANGE ANYTHING; BUT STARFISH VERSION IS INCOMPATIBLE)
##mafft conda-forge::mafft
##mmseqs2 bioconda::mmseqs2 (MUCH LATER VERSION OF TOOL THAN SPECIFIC IN STARFISH ENV; ONE SPECIFIED BY STARFISH IS INCOMPATIBLE)
##eggnog-mapper bioconda::eggnog-mapper

##created initially as such
#mamba create -n stargraph bioconda::pggb bioconda::bedtools bioconda::mash bioconda::seqkit bioconda::minimap2 bioconda::mummer4 samtools=1.6 bioconda::sourmash bioconda::blast bioconda::mcl=14.137 bioconda::circos bioconda::hmmer bioconda::metaeuk conda-forge::mafft bioconda::mmseqs2 bioconda::eggnog-mapper


##tried installing a starfish joint env but it does not jive
#mamba create -c egluckthaler -n stargraph bioconda::pggb bioconda::bedtools bioconda::mash bioconda::seqkit starfish
##has to install a MUCH earlier version of pggb due to some conflicts


##############################################################
################ STEP 0a: SETTING VARIABLES ##################
##############################################################

#default values, unless denoted when running stargraph
assemblies=""
tyrRs=""
elements=""
threads="1"
identity=""
length="20000"
kmersize="19"
separator="_"
minsize="30000"
prefix="stargraph"

output="stargraph_output"
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-a|--assemblies)
	assemblies="$2"
	shift
	shift
	;;
	-r|--tyrRs)
	tyrRs="$2"
	shift
	shift
	;;
	-e|--elements)
	elements="$2"
	shift
	shift
	;;
	-t|--threads)
	threads="$2"
	shift
	shift
	;;
	-i|--identity)
	identity="$2"
	shift
	shift
	;;
	-l|--length)
	length="$2"
	shift
	shift
	;;
	-k|--kmersize)
	kmersize="$2"
	shift
	shift
	;;
	-s|--separator)
	separator="$2"
	shift
	shift
	;;
	-m|--minsize)
	minsize="$2"
	shift
	shift
	;;
	-p|--prefix)
	prefix="$2"
	shift
	shift
	;;
	-o|--output)
	output="$2"
	shift
	shift
	;;
	-c|--cleanup)
	cleanup="$2"
	shift
	shift
	;;
	-h|--help)
	echo "
	
	stargraph (version: ${version})
 
	stargraph.sh -a assemblies_panSN.txt
	
	Required inputs:
	-a | --assemblies		A txt file with each line containing the path to an assembly using the PanSN-spec-like naming scheme for each contig ([sample][delim][contig/scaffold])
	-r | --tyrRs			Output file from starfish annotate that contains locations for the tyrosine recombinases in all assemblies (geneFinder/*.filt.gff)
	-e | --elements			Output file from starfish insert (preferably manually curated) that contains locations for the Starships (elementFinder/*.elements.ann.feat)


	Recommended inputs:
	-t | --threads			Number of threads for tools that accept this option (default: 1)

	pggb specific inputs:
	-i | --identity			-p option in pggb (Default: Automatically calulated using mash distances)
	-l | --length			-s option in pggb (Default: 20000 ; a conservative value increased from default pggb values)
	-k | --kmersize			-k option in pggb (Default: 19 ; same as pggb)


	Optional parameters:
	-s | --separator		PanSN-spec-like naming separator used (Default: _)
	-m | --minsize			Minimum size of PAVs to be kept (Default: 30000)
	-p | --prefix			Prefix for output (Default: stargraph)
	-o | --output			Name of output folder for all results (Default: stargraph_output)
	-c | --cleanup			Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help			Print this help message

	"
	exit
	;;
	esac
done


#creates error message and exits if these values are not/incorrectly assigned 
[[ $assemblies == "" ]] && echo "ERROR: No file containing assembly paths was given" && exit
#[[ $tyrRs == "" ]] && echo "ERROR: No file containing tyrosine recombinases locations was given" && exit
#[[ $elements == "" ]] && echo "ERROR: No file containing Starship elements locations was given" && exit

[[ $separator == "_" ]] && echo "WARNING: Running with default separator '_' (this is recommended but just a heads up)"
[[ $identity == "" ]] && echo "WARNING: Running with no identity given and therefore will calculate divergence using mash distance (this is recommended but just a heads up)"
[[ $threads == "1" ]] && echo "WARNING: Running with default thread count; depending on the dataset this may take some time"

##############################################################
####################### 0b. SETTING UP #######################
##############################################################

##redefining some variables and checking

##paths to raw data
assembliespath=$( realpath ${assemblies} )

#check if the files given actually exist
[ ! -f "${assembliespath}" ] && echo "ERROR: Cannot find path to assemblies file provided by -a; check path is correct and file exists" && exit

##check if output directory already exists
[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

##create output directory
mkdir ${output}
##get the absolute paths to each of the assemblies and save in a new file
cat ${assembliespath} | while read genome
do
realpath ${genome}
done > ${output}/path_to_assemblies.txt
##move into output directory
cd ${output}

##newly defining some variables for output writing
## leave minimum PAV size as is just setting to kb for some output variables
minsize2=$( echo ${minsize} | awk '{print $1/1000}' )


#####################################################################
################# STEP 1: CREATING THE GENOME GRAPH #################
#####################################################################


##calculate the minimum nucleotide identity between all assemblies provided
##this will made slightly more lentient (adding 5% more divergence) and provided as the -p option in pggb (if not provided manually)
##use mash triangle to generate a all vs all comparison per assembly file
##then use this as a proxy for nucleotide identity distances and add an increase of 5% divergence to be used for the -p option in pggb (can be quite lenient here)
mash triangle -l path_to_assemblies.txt -E > ${prefix}.assemblies.mash_distances.txt 
if [[ ${identity} == "" ]]
then
identity=$( cat ${prefix}.assemblies.mash_distances.txt | cut -f3 | awk '{print 100-$1}' | sort -n | head -n1 | awk -F "." '{print $1-5}' ) 
fi


##concatenate and compress the assemblies using the absolute paths to each assembly
##and removing all soft-masking
cat path_to_assemblies.txt | while read genome
do
if [[ ${genome} =~ ".gz"$ || ${genome} =~ ".bgzip"$ || ${genome} =~ ".gzip"$ ]]
then
zcat ${genome} | awk '{if($1 ~ ">") {print } else {print toupper($0)}}'
else
cat ${genome}  | awk '{if($1 ~ ">") {print } else {print toupper($0)}}'
fi 
done | bgzip > ${prefix}.assemblies.fa.gz

##get the number of assemblies provided to be used as the -n option in pggb
genomecount=$( cat path_to_assemblies.txt | wc -l )
##index the concatenated assemblies file with samtools 
samtools faidx ${prefix}.assemblies.fa.gz

##now build the genome graph
pggb -i ${prefix}.assemblies.fa.gz -o ${prefix}.pggb -t ${threads} -p ${identity} -s ${length} -m -S -n ${genomecount} -k ${kmersize} -Y ${separator}

##move in the output folder so we can try use this graph to extract Starship regions
cd ${prefix}.pggb


################################################################################
################# STEP 2: IDENTIFYING PRESENCE/ABSENCE VARIANT #################
################################################################################


##the simpliest way to do so was to use the odgi presence-absence 'PAV' function 
##it can produce a matrix of whether regions are present or absent in for each strain
##however we need to define the 'regions' for this calculation (i.e the windows as used in coverage)

##first we need to generate a file that associates each path (contigs) with a sample (which in our case is always the first string, split by the separator, in the contig name)
##we wil extract this directly from the odgi file just to be sure
odgi paths -i ${prefix}.assemblies.fa.gz.*.smooth.final.og -L > ${prefix}.paths.txt
cut -f 1 -d '#' ${prefix}.paths.txt > ${prefix}.samples.txt
paste ${prefix}.paths.txt ${prefix}.samples.txt > ${prefix}.path_and_sample.txt

##Now you can just use bins as done for coverage analysis and then give this to the pav analysis
##here I use a 1kb window
bedtools makewindows -g <(cut -f 1,2 ../${prefix}.assemblies.fa.gz.fai) -w 1000 | sed 's/\t0\t/\t1\t/g' > ${prefix}.assemblies.w1kb.bed
odgi pav -i ${prefix}.assemblies.fa.gz.*.smooth.final.og -b ${prefix}.assemblies.w1kb.bed -M -p ${prefix}.path_and_sample.txt > temp

##extract the header, sort the PAV rows then recombine
head -n1 temp > temp2
tail -n+2 temp | sort -k1,1 -k2,2V > temp3
cat temp2 temp3 > ${prefix}.odgi_pav.matrix.tsv
rm temp*

##then we can extract the actual regions (not restricted by the windows), where we can use the absence of any region in any strain (as this is an all-v-all comp)
##this step is only combining windows WITHIN EACH STRAIN-STRAIN COMPARISON
##so it melts the matrix
##first we combine all 1bp overlapping regions with less than 50% covered
##then remove all regions smaller than 10kb (trying to eliminate transposon impact)
##afterwards a second merge allows for 20kb gaps; this distance was based on the distribution of the nearest neighbour for regions which indentified a peak around 20kb so is targetting gaps appear more frequent that what we should expect
##these gaps that are merged generally represent transposon impacted regions that break up previously continuous regions (Can see this in how Starships are often fragmented in assemblies, even some long-read assemblies)
##Also, removing the 10kb regions, prior to this 20kb gap merging, reduces the chances of many little TE fragmemts being combined over large distances 
for i in $( seq 1 1 ${genomecount} )
do
cat ${prefix}.odgi_pav.matrix.tsv | awk -v row="$i" '{if(NR == 1) {header=$(row+4)} else {print $1"\t"$2"\t"$3"\t"$(row+4)"\t"header}}' | awk '{if($4 < 0.5) print}' | bedtools merge -d 1 -c 5 -o distinct -delim ";" | awk '{if($3-$2 > 10000) print}' | sort -k1,1 -k2,2n  | bedtools merge -d 20000 -c 4 -o distinct -delim ";"
done > ${prefix}.PAVs.tsv

##This second step merges all the STRAIN-v-STRAIN PAVs into a single PAV call at each genomic loci
##we merge all the regions into a single call per backbone, keeping track of all the strains within this region and only allowing for a single base gap (previously I allowed for more of a gap but it never made much sense)
##here we can also apply a cut off for the size of the regions that will be of real interest (-m --minsize)

##filter the ROI file and merge accross genomes
cat ${prefix}.PAVs.tsv | awk -v minsize="$minsize" '{if($3-$2 > minsize) print}' | sort -k1,1 -k2,2n | bedtools merge -d 1 -c 4 -o distinct -delim ";"  > ${prefix}.PAVs.${minsize2}kb_min.tsv

##we can now extract the actual region from the genome and create a fasta file of these PAVs
cat ${prefix}.PAVs.${minsize2}kb_min.tsv | awk '{print $1":"$2"-"$3}' | while read region
do
samtools faidx ../${prefix}.assemblies.fa.gz "${region}"
done > ${prefix}.PAVs.${minsize2}kb_min.fa


##how about some stats about these PAVs, per strain
##here we are interested in a few features
##the number of PAVs (however this can be impacted by genome assembly contiguity)
##the total number of bases in the PAVs (less impacted by contiguity so a better measure for comparing these raw results however again these will be filtered for Starship-related regions later)
##we can also look at the average and max size of the events
echo "strain;count;sum_length;mean_length;max_length" | sed 's/;/\t/g' > ${prefix}.PAVs.${minsize2}kb_min.stats.tsv
cat ${prefix}.PAVs.${minsize2}kb_min.tsv | cut -f1 | awk -F "${separator}" '{print $1}' | sort -u | while read sample
do
grep ^$sample ${prefix}.PAVs.${minsize2}kb_min.tsv | awk -v sample="$sample" -v separator="$separator" '{if($1 ~ sample""separator && max < ($3-$2)) {count++ ; sum=sum+($3-$2); max=($3-$2)} else if($1 ~ sample""separator && max > ($3-$2)) {count++ ; sum=sum+($3-$2)}} END{print sample"\t"count"\t"sum"\t"sum/count"\t"max}'
done >> ${prefix}.PAVs.${minsize2}kb_min.stats.tsv



################################################################################
################# STEP 3a: IDENTIFYING STARSHIP-LIKE REGIONS ###################
################################################################################


##this step requires starfish data input
##it will use de-novo annotation of captains (and possibly more starship related genes) in order to predict whether regions are starship-like regions

##need to specify that the tyr tag in the gff file is the the genes to be used for identifying the high-confidence SLRs


##identify lower confidence SLRs using, i.e. regions not containing a captain but containing other Starship-related genes



##we can use mash again to find which SLRs are within 90% similarity
##we can then use this clustering in order to identify 'haplotypes' like with starfish and therefore plot haplotype alignments/insertions below




###################################################################
################# STEP 3b: GENERATING SLR PLOTS ###################
###################################################################







####TO USE LATER ON IF I USE THE ASSEMBLY FILE TO EXTRACT SOMETHING, NEED TO KNOW IF COMPRESSED ESSENTIALLY
#if [[ ${assemblies} =~ ".gz"$ || ${assemblies} =~ ".bgzip"$ || ${assemblies} =~ ".gzip"$ ]]
#then
#XXXX
#else
#XXXX
#fi

