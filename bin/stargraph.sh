#!/bin/bash
set -euo pipefail

version="v1"

#####################################################################################
############# STEP -1. CREATING THE ENVIRONMENT. NOT ACTUALLY USED NOW ##############
#####################################################################################

### need to create a conda environment for housing all the stargraph dependencies:
## TRIED INSTALLING STARFISH IN THE SAME ENVIRONMENT BUT THE DEPENDENCIES AND NAMING MAKE IT IMPOSSIBLE ESSENTIALLY...
## therefore the starfish portion has to be run outside of stargraph and only input from that pipeline will be accepted

##starfish (egluckthaler::starfish) NOT USED (trying to install it directly with the stargraph conda env; see below)
##pggb/odgi (bioconda::pggb)
##bedtools (bioconda::bedtools)
##mash	(bioconda::mash)
##sourmash plugin (conda-forge::sourmash_plugin_branchwater) sourmash is already required with starfish
##seqkit (bioconda::seqkit)
##minimap2 (biocojnda::minimap2)
##mummer4 (bioconda::mummer4)
##R (conda-forge::r-base)
##gggenomes (conda-forge::r-gggenomes)
##ggnewscale (conda-forge::r-ggnewscale)
##IRanges (bioconda::bioconductor-iranges)
##svglite (conda-forge::r-svglite)


##created initially as such
#mamba create -n stargraph bioconda::pggb bioconda::bedtools bioconda::mash bioconda::seqkit bioconda::minimap2 bioconda::mummer4 conda-forge::r-base conda-forge::r-gggenomes conda-forge::r-ggnewscale bioconda::bioconductor-iranges conda-forge::r-svglite samtools=1.6 bioconda::sourmash bioconda::blast bioconda::mcl=14.137 bioconda::circos bioconda::hmmer bioconda::metaeuk conda-forge::mafft bioconda::mmseqs2 bioconda::eggnog-mapper


##tried installing a starfish joint env but it does not jive
#mamba create -c egluckthaler -n stargraph bioconda::pggb bioconda::bedtools bioconda::mash bioconda::seqkit starfish
##has to install a MUCH earlier version of pggb due to some conflicts and also R conflicts

##therefore added the required conda packages manually (using starfish yaml) and then added starfish github manually during conda package creation (therefore not technically installaing the tool and avoiding naming conflicts)
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


##############################################################
################ STEP 0a: SETTING VARIABLES ##################
##############################################################

#default values, unless denoted when running stargraph
assemblies=""
tyrRs=""
elements=""
threads="1"
identifier="tyr"
identity=""
length="20000"
kmersize="19"
separator="_"
minsize="30000"
window="1000"
flank="75000"
poaparam="7919,8069"
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
	-i|--identifier)
	identifier="$2"
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
	-G|--poaparam)
	poaparam="$2"
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
	-w|--window)
	window="$2"
	shift
	shift
	;;
	-f|--flank)
	flank="$2"
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
 
	stargraph.sh -a assemblies_panSN.txt -r starfish_output/geneFinder/*.filt.gff -e starfish_output/elementFinder/*.elements.ann.feat
	
	Required inputs:
	-a | --assemblies		A txt file with each line containing the path to an assembly using the PanSN-spec-like naming scheme for each contig ([sample][delim][contig/scaffold])
	-r | --tyrRs			Output file from starfish annotate that contains locations for the tyrosine recombinases in all assemblies (geneFinder/*.filt.gff)
	-e | --elements			Output file from starfish insert (preferably manually curated) that contains locations for the Starships (elementFinder/*.elements.ann.feat)


	Recommended inputs:
	-t | --threads			Number of threads for tools that accept this option (default: 1)
	-i | --identifier		The identifying tag used for tyrosine recombinases; given as the -i option for starfish annotate (Default: tyr)

	pggb specific inputs:
	-i | --identity			-p option in pggb (Default: Automatically calulated using mash distances)
	-l | --length			-s option in pggb (Default: 20000 ; a conservative value increased from default pggb values)
	-k | --kmersize			-k option in pggb (Default: 19 ; same as pggb)
	-G | --poaparam			-G option in pggb (Default: 7919,8069; a conservative value increased from default pggb values)

	Optional parameters:
	-s | --separator		PanSN-spec-like naming separator used (Default: _)
	-m | --minsize			Minimum size of PAVs to be kept (Default: 30000)
	-w | --window			Size of windows used for PAV detection (Default: 1000)
	-f | --flank			Size of flanking region used when plotting element alignments (Default: 75000)
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
[[ $tyrRs == "" ]] && echo "ERROR: No file containing tyrosine recombinases locations was given" && exit
[[ $elements == "" ]] && echo "ERROR: No file containing Starship elements locations was given" && exit

[[ $separator == "_" ]] && echo "WARNING: Running with default separator '_' (this is recommended but just a heads up)"
[[ $identity == "" ]] && echo "WARNING: Running with no identity given and therefore will calculate divergence using mash distance (this is recommended but just a heads up)"
[[ $threads == "1" ]] && echo "WARNING: Running with default thread count; depending on the dataset this may take some time"

##############################################################
####################### 0b. SETTING UP #######################
##############################################################

##redefining some variables and checking

##paths to raw data
assembliespath=$( realpath ${assemblies} )
tyrRspath=$( realpath ${tyrRs} )
elementspath=$( realpath ${elements} )

#check if the files given actually exist
[ ! -f "${assembliespath}" ] && echo "ERROR: Cannot find path to assemblies file provided by -a; check path is correct and file exists" && exit

##check if output directory already exists
[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

##create output directory
mkdir ${output}
outputpath=$( realpath ${output} )
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
## leave window size size as is just setting to kb for some output variables
window2=$( echo ${window} | awk '{print $1/1000}' )
## leave flank size size as is just setting to kb for some output variables
flank2=$( echo ${flank} | awk '{print $1/1000}' )

#####################################################################
################# STEP 1: CREATING THE GENOME GRAPH #################
#####################################################################


##calculate the minimum nucleotide identity between all assemblies provided
##this will made slightly more lentient (adding 5% more divergence) and provided as the -p option in pggb (if not provided manually)
##use mash triangle to generate a all vs all comparison per assembly file
##then use this as a proxy for nucleotide identity distances and add an increase of 5% divergence to be used for the -p option in pggb (can be quite lenient here)
mash triangle -s 100000 -l path_to_assemblies.txt -E > ${prefix}.assemblies.mash_distances.txt 
if [[ ${identity} == "" ]]
then
identity=$( cat ${prefix}.assemblies.mash_distances.txt | cut -f3 | awk '{print (1-$1)*100}' | sort -n | head -n1 | awk -F "." '{print $1-5}' ) 
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
pggb -i ${prefix}.assemblies.fa.gz -o 1.${prefix}.pggb -t ${threads} -p ${identity} -s ${length} -m -S -n ${genomecount} -k ${kmersize} -Y ${separator} -G ${poaparam}

##move in the output folder so we can try use this graph to extract Starship regions
cd 1.${prefix}.pggb


################################################################################
################# STEP 2: IDENTIFYING PRESENCE/ABSENCE VARIANT #################
################################################################################


##the simpliest way to do so was to use the odgi presence-absence 'PAV' function 
##it can produce a matrix of whether regions are present or absent in for each strain
##however we need to define the 'regions' for this calculation (i.e the windows as used in coverage)

##first we need to generate a file that associates each path (contigs) with a sample (which in our case is always the first string, split by the separator, in the contig name)
##we wil extract this directly from the odgi file just to be sure
odgi paths -i ${prefix}.assemblies.fa.gz.*.smooth.final.og -L > ${prefix}.paths.txt
cat ${prefix}.paths.txt | awk -F "${separator}" '{print $1}' > ${prefix}.samples.txt
paste ${prefix}.paths.txt ${prefix}.samples.txt > ${prefix}.path_and_sample.txt

##Now you can just use bins as done for coverage analysis and then give this to the pav analysis
##here I use the window parameter setting for creating the window size in a bed file that will be analysed for presence/absence
bedtools makewindows -g <(cut -f 1,2 ../${prefix}.assemblies.fa.gz.fai) -w ${window} | sed 's/\t0\t/\t1\t/g' > ${prefix}.assemblies.w${window2}kb.bed
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

cd ..

################################################################################
################# STEP 3a: IDENTIFYING STARSHIP-LIKE REGIONS ###################
################################################################################


mkdir 2.PAVs_to_SLRs
cd 2.PAVs_to_SLRs

##this step requires starfish data input
##it will use de-novo annotation of captains (and possibly more starship related genes) in order to predict whether regions are starship-like regions
##it will use all Starship-related genes (SRGs) given initially to identify all PAVs with SRGs present
##each line is a element and an associated SRG
echo "contig;start;end;SRG_name;SRG_start;SRG_end;SRG_length;SRG_sense" | tr ';' '\t' > ${prefix}.PAVs.all_SRGs.tsv
cat ../1.${prefix}.pggb/${prefix}.PAVs.${minsize2}kb_min.tsv | while read element
do
sample=$( echo "${element}" | awk -F "${separator}" '{print $1}' )
contig=$( echo "${element}" | awk -F "\t" '{print $1}' )
start=$( echo "${element}" | awk -F "\t" '{print $2}' )
end=$( echo "${element}" | awk -F "\t" '{print $3}' )
cat ${tyrRspath} | awk -v contig="$contig" -v start="$start" -v end="$end" -v element="$element" -v sample="$sample" '{if($1 == contig && $4 > start && $5 < end) {print contig"\t"start"\t"end"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}}' | sort -k4n
done | awk -F "Name=" '{print $1"\t"$NF}' | awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7}' >> ${prefix}.PAVs.all_SRGs.tsv

##next the list will be culled down to important SRGs (if there is any distinction to be made)
##need to specify that the 'tyr' tag in the gff file is the the genes to be used for identifying the high-confidence SLRs
##if "tyr" is not given then just all the SRGs given
##will also reduce this list down to a single line per element
echo "SLR;contig;start;end" | tr ';' '\t' > ${prefix}.SLRs.tsv
captains=$( cat ${prefix}.PAVs.all_SRGs.tsv | awk -v identifier="$identifier" '{if($4 ~ "_"identifier) print "present"}' | sort -u )
if [[ $captains == "present" ]]
then
tail -n+2 ${prefix}.PAVs.all_SRGs.tsv | awk -v identifier="$identifier" '{if($4 ~ "_"identifier) print $1"\t"$2"\t"$3}' | sort -u | awk -F "${separator}" '{print $1"\t"$0}' | awk -v count="1" '{print $1"_SLR"count"\t"$2"\t"$3"\t"$4; count++}' >> ${prefix}.SLRs.tsv
else
tail -n+2 ${prefix}.PAVs.all_SRGs.tsv | awk '{print $1"\t"$2"\t"$3}' | sort -u | awk -F "${separator}" '{print $1"\t"$0}' | awk -v count="1" '{print $1"_SLR"count"\t"$2"\t"$3"\t"$4; count++}' >> ${prefix}.SLRs.tsv
fi 

##now get each SLR; get the SRGs associated with it and then classify some important SRGs
##for now the classification is only for those tagged previously as tyr; duf3723 and myb (as done in the starfish wrapper); all others are given NA but this could be manually modified
##these classifications will just be used for plotting
echo "SLR;contig;start;end;SRG_name;SRG_start;SRG_end;SRG_sense;SRG_class" | tr ';' '\t' > ${prefix}.SLRs.plus_SRGs.tsv
tail -n+2 ${prefix}.SLRs.tsv | while read line
do
SLR=$( echo "${line}" | awk '{print $1}' )
contig=$( echo "${line}" | awk '{print $2}' )
start=$( echo "${line}" | awk '{print $3}' )
end=$( echo "${line}" | awk '{print $4}' )
cat ${prefix}.PAVs.all_SRGs.tsv | awk -v SLR="$SLR" -v contig="$contig" -v start="$start" -v end="$end" '{if(contig == $1 && start == $2 && end == $3) {print SLR"\t"contig"\t"start"\t"end"\t"$4"\t"$5"\t"$6"\t"$8}}'
done | awk -v identifier="$identifier" '{if($5 ~ "_"identifier) {print $0"\ttyrR"} else if($5 ~ "_myb") {print $0"\tmyb"} else if($5 ~ "_duf3723") {print $0"\tduf3723"} else {print $0"\tNA"}}' >> ${prefix}.SLRs.plus_SRGs.tsv

##create a fasta file with just the PAVs
tail -n+2 ${prefix}.SLRs.tsv | while read line
do
SLR=$( echo "${line}" | awk '{print $1}'  )
coords=$( echo "${line}" | awk '{print $2":"$3"-"$4}' )
samtools faidx ../${prefix}.assemblies.fa.gz "${coords}" | awk -v SLR="$SLR" '{if($0 ~ ">") {print ">"SLR} else {print}}'
done > ${prefix}.SLRs.fa

##we can use mash again to find which SLRs are clustering together based on overlapping kmers
##sourmash alternative, sketch the signatures (giving a small scalling value (100bp level of detection), small kmer, and sketch for each nucleotide sequence in the file)
sourmash sketch dna -p scaled=100,k=21 ${prefix}.SLRs.fa --singleton -o ${prefix}.SLRs.sig
##now compare it against itself using 'containment' as the metrix (this allows us to easily find smaller elements nested within larger ones)
sourmash compare ${prefix}.SLRs.sig --containment --csv ${prefix}.SLRs.sig.compare.csv --labels-to ${prefix}.SLRs.sig.compare.txt
##convert to pairwise comparisons
cat ${prefix}.SLRs.sig.compare.csv | tr -d '\r'  | awk -F',' 'NR==1{for(i=1;i<=NF;i++)samples[i]=$i;next}{row=NR-1;for(i=row+1;i<=NF;i++)print samples[row],samples[i],$i}' OFS='\t' >  ${prefix}.SLRs.sig.pairwise.tsv
##now use mcl to quickly find the clusters
mcl ${prefix}.SLRs.sig.pairwise.tsv --abc -o ${prefix}.SLRs.sig.pairwise.mcl.txt
##now name the clusters and then append to the summary files
awk -F '\t' '{for (i=1; i <= NF; i++) {print "cluster"NR "\t" $i}}' ${prefix}.SLRs.sig.pairwise.mcl.txt > ${prefix}.SLRs.sig.pairwise.mcl.clusters.txt

echo "SLR;contig;start;end;cluster" | tr ';' '\t' > ${prefix}.SLRs.plus_clusters.tsv
cat ${prefix}.SLRs.tsv | while read line
do
SLR=$( echo "${line}" | awk '{print $1}' )
cat ${prefix}.SLRs.sig.pairwise.mcl.clusters.txt| awk -v SLR="$SLR" -v line="$line" '{if($2==SLR) print line"\t"$1}'
done >> ${prefix}.SLRs.plus_clusters.tsv

cd ..

###################################################################
################# STEP 3b: GENERATING SLR PLOTS ###################
###################################################################

mkdir 3.SLR_plots
cd 3.SLR_plots

##first create a simplified bed file for the annotations
##create header for annotations
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > genes.bed
##now grab the coordinates, sense, name of gene and a label just for tyrs, duf3723 and mybs
cat ${tyrRspath} | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$7"\t"$9}' | awk -F ";Name=" '{print $0"\t"$NF}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' | awk -v identifier="$identifier"  '{if($5 ~ "_"identifier) {print $0"\ttyrR"} else if($5 ~ "_myb") {print $0"\tmyb"} else if($5 ~ "_duf3723") {print $0"\tduf3723"} else {print $0"\tNA"}}' >> genes.bed


##loop through each cluster defined previously to get a list of the SLRs and a list of genomes missing these elements to be used as an insertion site
tail -n+2 ../2.PAVs_to_SLRs/${prefix}.SLRs.plus_clusters.tsv | cut -f5 | sort -u | while read cluster
do

##get the SLRs
cat ../2.PAVs_to_SLRs/${prefix}.SLRs.plus_clusters.tsv | awk -v cluster="$cluster" '{if($5 == cluster) {print $1}}' > ${cluster}.list.txt

##get the genomes missed these SLRs
cat ../2.PAVs_to_SLRs/${prefix}.SLRs.plus_clusters.tsv | awk -v cluster="$cluster" '{if($5 == cluster) print}' | while read line
do
SLR=$( echo "${line}" | awk '{print $1}'  )
position=$( echo "${line}" | awk '{print $2"\t"$3"\t"$4}'  )
cat ../1.${prefix}.pggb/${prefix}.PAVs.${minsize2}kb_min.tsv | grep "${position}" | awk '{print $4}' | tr ';' '\n' | while read absent
do
cat ${cluster}.list.txt | awk -F "${separator}" -v absent="$absent" -v hold="absent" '{if($1 == absent) {hold="present"}} END{if(hold=="absent") print absent}' 
done > ${cluster}.absent.txt
done

###extract the SLRs PLUS the flanking regions around them into a single fasta for alignment (this will be used to identify a good region to visualise the insertion)
##also extract the full contigs in which the SLRs are found (this will be used for the actual alignment)
if [ -f ${cluster}.regions_plus_flank.fa ]
then
rm ${cluster}.regions_plus_flank.fa
fi

if [ -f ${cluster}.contigs.fa ]
then
rm ${cluster}.contigs.fa
fi

echo "contig;start;end;SLR" | tr ';' '\t' > ${cluster}.regions_plus_flank.tsv
cat ${cluster}.list.txt | while read SLR
do
contig=$( cat ../2.PAVs_to_SLRs/${prefix}.SLRs.tsv | awk -F "\t" -v SLR="$SLR" '{if($1 == SLR) print $2}' )
contiglength=$( cat ../${prefix}.assemblies.fa.gz.fai | awk -v contig="$contig" '{if($1 == contig) print $2}' )
cat ../2.PAVs_to_SLRs/${prefix}.SLRs.tsv | awk -F "\t" -v SLR="$SLR" -v flank="$flank" '{if($1 == SLR) print $2"\t"$3-flank"\t"$4+flank}' | awk '{if($2 < 0 ) {print $1"\t1\t"$3} else {print}}' | awk -v SLR="$SLR" -v contiglength="$contiglength" '{if($3 > contiglength ) {print $1"\t"$2"\t"contiglength"\t"SLR} else {print $0"\t"SLR}}'
done >> ${cluster}.regions_plus_flank.tsv
tail -n+2 ${cluster}.regions_plus_flank.tsv  | while read line
do
coords=$( echo "${line}" | awk '{print $1":"$2"-"$3}' )
contig=$( echo "${line}" | awk '{print $1}' )
SLR=$( echo "${line}" | awk '{print $4}' )
samtools faidx ../${prefix}.assemblies.fa.gz "${coords}" | awk -v SLR="$SLR" '{if($0 ~ ">") {print ">"SLR">>"$0} else {print}}' | sed 's/>>>/ /g' >> ${cluster}.regions_plus_flank.fa
samtools faidx ../${prefix}.assemblies.fa.gz "${contig}" >> ${cluster}.contigs.fa
done 

##take just one of the SLRs in the cluster; one with the largest sum of flank lengths on either side or the largest (only one if they are equal)
topSLR=$( cat ${cluster}.list.txt | while read SLR
do
start=$( cat ../2.PAVs_to_SLRs/${prefix}.SLRs.tsv | awk -F "\t" -v SLR="$SLR" '{if($1 == SLR) print $3}' )
end=$( cat ../2.PAVs_to_SLRs/${prefix}.SLRs.tsv | awk -F "\t" -v SLR="$SLR" '{if($1 == SLR) print $4}' )

startmoddiff=$( cat ${cluster}.regions_plus_flank.tsv | awk -F "\t" -v SLR="$SLR" -v start="$start" '{if($4 == SLR) print start-$2}' )
endmoddiff=$( cat ${cluster}.regions_plus_flank.tsv | awk -F "\t" -v SLR="$SLR" -v end="$end" '{if($4 == SLR) print $3-end}' )

echo "${SLR};${startmoddiff};${endmoddiff}" | tr ';' '\t'

done | awk -v flank="$flank" '{if(($2+$3) > sumflank) {sumflank=($2+$3); SLR=$1}} END{print SLR}' )
##save the SLR+flank region to a temporary fasta file
samtools faidx ${cluster}.regions_plus_flank.fa ${topSLR} > ${cluster}.regions_plus_flank.temp.fa

##now align this region to the genomes where this SLR is supposed to be absent
##first extract the genomes into a seperate fasta for alignment (removing any contigs smaller than a single flank size)

cat ${cluster}.absent.txt | while read othergenome
do
zcat ../${prefix}.assemblies.fa.gz | grep ">"${othergenome}"${separator}" | sed 's/>//g' | while read contig
do
length=$( cat ../${prefix}.assemblies.fa.gz.fai | awk -v contig="$contig" '{if($1 == contig) print $2}' )
if [[ ${length} -gt ${flank} ]]
then
samtools faidx ../${prefix}.assemblies.fa.gz ${contig}
fi
done
done > ${cluster}.absent.fa

##create a header for a bed file to be populated
##this bed file will dictate the regions to be aligned (which will be relative to the extracted regions)
echo "contig;start;end;SLR" | tr ';' '\t' > ${cluster}.regions_plus_flank.plotting.bed

##now align the SLR to the empty contigs (here only the flanks should have large aligning regions)
nucmer -t ${threads} --maxmatch --minmatch 100 --delta  ${cluster}.regions_plus_flank.absent.nucmer.delta ${cluster}.regions_plus_flank.temp.fa ${cluster}.absent.fa
paftools.js delta2paf ${cluster}.regions_plus_flank.absent.nucmer.delta > ${cluster}.regions_plus_flank.absent.nucmer.paf

##now use the paf file to find a contigs with good aligning regions
##good aligning can be that at least a single contig has a single alignment larger than half the flank region
##take the best two aligning contigs from the dataset (ideally it'll be large contigs from two different genomes)
##then try to find the edges of the alignments using 20kb seeds (this will be used for the plotting; i.e. only this regions alignment visualised)
if [ -f ${cluster}.absent.contigs.fa ]
then
rm ${cluster}.absent.contigs.fa
fi
cat ${cluster}.regions_plus_flank.absent.nucmer.paf | awk -v flank="${flank}" '{if($11 > (flank/2) ) print}' | cut -f1 | sort -u | while read tempcontig
do
cat ${cluster}.regions_plus_flank.absent.nucmer.paf | awk -v tempcontig="$tempcontig" '{if($1 == tempcontig) sum=sum+$11} END{print tempcontig"\t"sum}'
done | sort -k2nr | head -n2 | awk '{print $1}' | while read insertioncontig
do
strain=$( echo "${insertioncontig}" | awk -F "${separator}" '{print $1}' )
edges=$( cat ${cluster}.regions_plus_flank.absent.nucmer.paf | awk '{if($11 > 20000) print}' | awk -v insertioncontig="$insertioncontig" '{if($1==insertioncontig) print}' | awk -v flank="$flank" 'BEGIN{max=0; min=99999999999999} {if($4 > max) {max=$4}; if($3 < min) {min=$3}} END{print min-(flank/2)"\t"max+(flank/2)}' | awk '{if($1 < 0) {print 0"\t"$2} else {print}}' )
edges2=$( echo "${edges}" | awk '{print $1"-"$2}' )
size=$( echo "${edges}" | awk '{print $2-$1}' )
samtools faidx ${cluster}.absent.fa "${insertioncontig}" >> ${cluster}.absent.contigs.fa
echo "${insertioncontig};${edges};NA" | tr ';' '\t' >> ${cluster}.regions_plus_flank.plotting.bed
done 

##add the coordinates etc for the SLRs to the same bed file (doing it in this order so in the plotting the absent regions will be on top)
##then adding the topSLR first so that it will be plotted downstream and show the insertion site for this element
tail -n+2 ${cluster}.regions_plus_flank.tsv | awk -v topSLR="$topSLR" '{if($4 == topSLR) {print $1"\t"$2"\t"$3"\t"$4}}' >> ${cluster}.regions_plus_flank.plotting.bed
tail -n+2 ${cluster}.regions_plus_flank.tsv | awk -v topSLR="$topSLR" '{if($4 != topSLR) {print $1"\t"$2"\t"$3"\t"$4}}' >> ${cluster}.regions_plus_flank.plotting.bed


##create a simple bed file for the SLRs regions
echo "contig;start;end;SLR" | tr ';' '\t' > ${cluster}.SLRs.bed
cat ../2.PAVs_to_SLRs/${prefix}.SLRs.plus_clusters.tsv | awk -v cluster="${cluster}" '{if($5 == cluster) {print $2"\t"$3"\t"$4"\t"$1}}' >> ${cluster}.SLRs.bed



##add the contigs missing the elements to the complete contigs file
cat ${cluster}.absent.contigs.fa >> ${cluster}.contigs.fa


##generate all vs all alignments for the contigs
##remove self alignment and any alignment smaller than 1kb
nucmer --maxmatch --minmatch 100 --delta  ${cluster}.contigs.nucmer.delta ${cluster}.contigs.fa ${cluster}.contigs.fa
paftools.js delta2paf ${cluster}.contigs.nucmer.delta | awk -F "\t" '{if($1 != $6) {print}}' > ${cluster}.contigs.nucmer.paf


##automate the production of an R script using gggenomes to plot the alignment
##then use gggenome with R script to create the plots
Rscriptpath=$( which gggenomes_skeleton.R )
cat ${Rscriptpath} | sed "s/CLUSTER/${cluster}/g" | sed "s|PATHTOOUTPUT|${outputpath}/3.SLR_plots|g" > ${cluster}.R
Rscript ${cluster}.R

done


cd ..

##########################################################################
################# STEP 4: COMBINING SLRs AND STARSHIPS ###################
##########################################################################


##we now have out SLR dataset  and can combine this with the Starship input from starfish
##in doing so we actually do three simple steps
##Step 1: subtract all the SLR regions that overlap with a starfish region
##Step 2: filter out any remaining SLR region smaller than the initial SLR minimum size
##Step 3: rename the remaining SLRs if split by starships into several smaller SLRs (above the minimum size)

mkdir 4.SLR_starship_combination
cd 4.SLR_starship_combination

##for this we digest the starship and SLR positions into two bed files for the filtering
tail -n+2 ../2.PAVs_to_SLRs/${prefix}.SLRs.plus_clusters.tsv  | cut -f1-4 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > ${prefix}.SLRs.bed
tail -n+2 ${elementspath} | cut -f1,4,6-7 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > ${prefix}.starships.bed


##now do the subtraction from the SLRs and rename if any SLRs were split with multiple sections still remaning
bedtools subtract -a ${prefix}.SLRs.bed -b ${prefix}.starships.bed | sortBed -i - | awk -v minsize="$minsize" '{if($3-$2 >= minsize) print}' | awk '{if(SLR==$4) {print contig"\t"start"\t"end"\t"SLR"_"n ; contig=$1 ; start=$2; end=$3; SLR=$4; line=$0; n++; hold="hold"} else if(SLR != $4 && hold=="hold") {print line"_"n; contig=$1 ; start=$2; end=$3; SLR=$4; line=$0; n=1; hold="nohold"} else {print line; contig=$1 ; start=$2; end=$3; SLR=$4; line=$0; n=1; hold="nohold"}} END{if(hold=="hold") {print line"_"n; contig=$1 ; start=$2; end=$3; SLR=$4; line=$0; n=1; hold="nohold"} else {print line; contig=$1 ; start=$2; end=$3; SLR=$4; line=$0; n=1; hold="nohold"}}' | tail -n+2 > ${prefix}.SLRs.starships_subtracted.bed

##now use this subtracted bedfile to create a starship compatable summary file with captain positions (if present in the new SLR chunk)
echo "SLR;navis-haplotype;contig;start;end;size;captain;captain_start;captain_end;captain_size;captain_sense"  | sed 's/;/\t/g' > ${prefix}.SLRs.starships_subtracted.tyrRs.tsv
cat ${prefix}.SLRs.starships_subtracted.bed | while read line
do
contig=$( echo "${line}" | awk '{print $1}' )
start=$( echo "${line}" | awk '{print $2}' )
end=$( echo "${line}" | awk '{print $3}' )
SLR=$( echo "${line}" | awk '{print $4}' )
SLRshort=$( echo "${SLR}" | awk -F "_" '{print $1}' )
cat ${tyrRspath} | awk -v contig="$contig" -v start="$start" -v end="$end" -v line="$line" '{if($1 == contig && $4 >= start && $5 <= end) {captain="found"; details=line"\t"$4"\t"$5"\t"$6"\t"$7";"$9";\n"details}} END{if(captain == "found") {print details} else {print line"\tNA\tNA\tNA\tNA;;NA"}}' | awk -F ";" '{ print $1"\t"$3}' | sed 's/Name=//g'
done | sort -u  | tail -n+2 | awk '{print $4"\tNA\t"$1"\t"$2"\t"$3"\t"($3-$2)"\t"$9"\t"$5"\t"$6"\t"$7"\t"$8}' >> ${prefix}.SLRs.starships_subtracted.tyrRs.tsv

##aggregate the TyrR coordinates in the SLRs (or just a single tyrR if that is what is present)
##using "_tyr" tag
echo "SLR;navis-haplotype;contig;start;end;size;captain;captain_start;captain_end;captain_size;captain_sense"  | sed 's/;/\t/g' > ${prefix}.SLRs.starships_subtracted.tyrRs_agg.tsv
tail -n+2 ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | cut -f1 | sort -u | while read SLR
do
rest=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" '{if($1 == SLR) {print}}' | cut -f1-6 | sort -u  )
name=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" -v identifier="$identifier"  '{if($1 == SLR && $7 ~ "_"identifier) {sum=sum";"$7}} END{if(sum=="") {print "NA"} else {print sum}}' | sed 's/^;//g' )
start=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" -v identifier="$identifier"  '{if($1 == SLR && $7 ~ "_"identifier) {sum=sum";"$8}} END{if(sum=="") {print "NA"} else {print sum}}' | sed 's/^;//g' )
end=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" -v identifier="$identifier"  '{if($1 == SLR && $7 ~ "_"identifier) {sum=sum";"$9}} END{if(sum=="") {print "NA"} else {print sum}}' | sed 's/^;//g' )
size=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" -v identifier="$identifier"  '{if($1 == SLR && $7 ~ "_"identifier) {sum=sum";"$10}} END{if(sum=="") {print "NA"} else {print sum}}' | sed 's/^;//g' )
sense=$( cat ${prefix}.SLRs.starships_subtracted.tyrRs.tsv | awk -v SLR="$SLR" -v identifier="$identifier"  '{if($1 == SLR && $7 ~ "_"identifier) {sum=sum";"$11}} END{if(sum=="") {print "NA"} else {print sum}}' | sed 's/^;//g' )
echo "${rest},${name},${start},${end},${size},${sense}" | tr ',' '\t'
done >> ${prefix}.SLRs.starships_subtracted.tyrRs_agg.tsv

##now combined this file with an equivalent file for starships
echo "starship_SLR;navis-haplotype;contig;start;end;size;captain;captain_start;captain_end;captain_sense" | tr ';' '\t' > ${prefix}.starships_SLRs.tsv
cat ${elementspath} | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$7-$6"\t"$5}' | while read line
do
tyr=$( echo "${line}" | awk -F "\t" '{print $NF}' )
cat ${tyrRspath} | grep "${tyr}"$ | awk -v line="$line" '{print line"\t"$4"\t"$5"\t"$7}'
done >> ${prefix}.starships_SLRs.tsv
tail -n+2 ${prefix}.SLRs.starships_subtracted.tyrRs_agg.tsv | cut -f1-9,11 >> ${prefix}.starships_SLRs.tsv

##get a fasta for all these sequences
tail -n+2 ${prefix}.starships_SLRs.tsv | awk '{print $3"\t"$4"\t"$5"\t"$1}' > ${prefix}.starships_SLRs.bed
bedtools getfasta -name -bed ${prefix}.starships_SLRs.bed -fi ../${prefix}.assemblies.fa.gz | sed 's/::/ /g' > ${prefix}.starships_SLRs.fa



##redo the clustering approach using SLRs and starships
##convert SLRs to the same navis-haplotype association if they cluster???





######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################


######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################












####TO USE LATER ON IF I USE THE ASSEMBLY FILE TO EXTRACT SOMETHING, NEED TO KNOW IF COMPRESSED ESSENTIALLY
#if [[ ${assemblies} =~ ".gz"$ || ${assemblies} =~ ".bgzip"$ || ${assemblies} =~ ".gzip"$ ]]
#then
#XXXX
#else
#XXXX
#fi

