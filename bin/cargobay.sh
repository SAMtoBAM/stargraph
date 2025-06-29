#!/bin/bash
set -euo pipefail

version="v1"

##############################################################
################ STEP 0a: SETTING VARIABLES ##################
##############################################################

#default values, unless denoted when running stargraph
elements=""
elementsbed=""
assemblies=""
gff3=""
metadata=""
threads="1"
separator="_"
identifier="tyr"
containment="0.5"
flank="50000"
prefix="cargobay"

output="cargobay_output"
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-e|--elements)
	elements="$2"
	shift
	shift
	;;
	-b|--elementsbed)
	elementsbed="$2"
	shift
	shift
	;;
	-a|--assemblies)
	assemblies="$2"
	shift
	shift
	;;
	-g|--gff3)
	gff3="$2"
	shift
	shift
	;;
	-m|--metadata)
	metadata="$2"
	shift
	shift
	;;
	-t|--threads)
	threads="$2"
	shift
	shift
	;;
	-s|--separator)
	separator="$2"
	shift
	shift
	;;
	-i|--identifier)
	identifier="$2"
	shift
	shift
	;;
	-c|--containment)
	containment="$2"
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
 
	cargobay.sh -e elements.fa -b elements.bed -a assemblies.fa -g annotation.gff3 -m metadata.tsv
	
	Required inputs:
	-e | --elements		A multifasta file containing all the elements (Starships and SLRs) to be searched for
	-b | --elementsbed	A bed file containing all the positions of the elements (Starships and SLRs)
	-a | --assemblies	A multifasta file containing all the assemblies used to detect the Starships and SLRs
	-g | --gff3			An annotation file containing all or a subset of genes of interest for plotting (at a minimum the de-novo annotated tyrRs)
	-m | --metadata		A tsv file containing metadata two columns; first column is the sample name and the second an NCBI genus species name (e.g. Aspergillus fumigatus)

	Recommended inputs:
	-t | --threads		Number of threads for tools that accept this option (default: 1)
	-s | --separator	Separator used to split sample and Starship/SLR names (Default: "_")
	-i | --identifier	Identifier for gff3 to find and highlight tyrosine recombinase genes (Default: 'tyr')

	Optional parameters:
	-c | --containment	The minimum proportion containment threshold for identifying candidates for HGT using the sourmash database (Default: 0.5)
	-f | --flank		Number of basepairs up and downstream of the element to be used for plotting (Default: 50000)
	-p | --prefix		Prefix for output (Default: cargobay)
	-o | --output		Name of output folder for all results (Default: cargobay_output)
	-c | --cleanup		Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help			Print this help message

	"
	exit
	;;
	esac
done

[[ $elements == "" ]] && echo "ERROR: No file containing Starship/SLR elements was given (-e)" && exit
[[ $elementsbed == "" ]] && echo "ERROR: No file containing Starship/SLR locations was given (-b)" && exit
[[ $assemblies == "" ]] && echo "ERROR: No file assemblies was given (-a)" && exit
[[ $gff3 == "" ]] && echo "ERROR: No annotations were given (-g)" && exit
[[ $metadata == "" ]] && echo "ERROR: No file metadata was given (-m)" && exit


##############################################################
####################### 0b. SETTING UP #######################
##############################################################

##redefining some variables and checking

##paths to raw data
elementspath=$( realpath ${elements} )
elementsbedpath=$( realpath ${elementsbed} )
assembliespath=$( realpath ${assemblies} )
gff3path=$( realpath ${gff3} )
metadatapath=$( realpath ${metadata} )




#check if the files given actually exist
[ ! -f "${elementspath}" ] && echo "ERROR: Cannot find path to elements file provided by -e; check path is correct and file exists" && exit
[ ! -f "${elementsbedpath}" ] && echo "ERROR: Cannot find path to elements bed file provided by -b; check path is correct and file exists" && exit
[ ! -f "${assembliespath}" ] && echo "ERROR: Cannot find path to assemblies file provided by -a; check path is correct and file exists" && exit
[ ! -f "${gff3path}" ] && echo "ERROR: Cannot find path to gff3 file provided by -g; check path is correct and file exists" && exit
[ ! -f "${metadatapath}" ] && echo "ERROR: Cannot find path to metadata file provided by -m; check path is correct and file exists" && exit


##check if output directory already exists
[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

echo "Welcome to CargoBay!"

[[ $separator == "_" ]] && echo "WARNING: Running with default separator '_' (this is recommended but just a heads up)"
[[ $containment == "0.5" ]] && echo "WARNING: Running with default containment threshold of 50% (just a heads up)"
[[ $threads == "1" ]] && echo "WARNING: Running with default thread count; depending on the dataset this may take some time"
[[ $identifier == "tyr" ]] && echo "WARNING: Running with default identifier for tyrosine recombinases 'tyr' (just a heads up)"


mkdir ${output}
outputpath=$( realpath ${output} )
cd ${output}


######################################################################################
####################### 0c. DOWNLOADING SOURMASH/NCBI DATABASE #######################
######################################################################################

echo "Step 0: Downloading the fungal sourmash database"

####SOURMASH DATABASES

##taken from this github issues post:
#https://github.com/sourmash-bio/sourmash/issues/3649
#The ones linked above are scaled=10k from Jan 2025.

#We also have draft versions from April 2025, scaled 1000, here:

#https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k21.zip
#https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k31.zip
#https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k51.zip

#and the corresponding lineages file: https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/lineages.fungi.csv

mkdir 0.database

##therefore using these links for now

##get the fungal lineages from the sourmash provided db
curl -JLO -s https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/lineages.fungi.csv 
##get the database from april 2025 with 1000 scalling and kmer 21
curl -JLO -s https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2025.04/genbank-20250408-fungi-k21.zip 

mv lineages.fungi.csv  0.database/
mv genbank-20250408-fungi-k21.zip 0.database/

#####################################################################################
####################### STEP 1a. KMER BASED DATABASE SEARCH #########################
#####################################################################################

echo "Step 1a: Sketching signatures for your input sequences then searching for them in the fungi database"

mkdir 1.database_search

##create a 1000 scaled; 21 kmer signature database (needs to be 1000 scaled due to the fungal database being built as so)
sourmash sketch dna -p scaled=1000,k=21 ${elementspath} --singleton -o 1.database_search/${prefix}.elements.sig

##then use gather to find similarities, use the multi version to search for multiple singatures simultaneously and then fast version to speed it up
##didn't use gather here as it was only giving a single match per signature file with all the details of the match
##the prefetch file automatically given per singature presented more matches but without the details of containment etc
#sourmash scripts fastmultigather test.SLRs.1000.sig 0.database/genbank-20250408-fungi-k21.zip -o test_fastgather -k 21 --cores ${threads}

##use multisearch instead which can be multithreaded and accept multiple singatures in a single file
sourmash scripts multisearch --quiet 1.database_search/${prefix}.elements.sig 0.database/genbank-20250408-fungi-k21.zip -o 1.database_search/${prefix}.sourmash_multisearch.csv -k 21 --cores ${threads}


#####################################################################################
###################### STEP 1b. FINDING LIKELY HGT CANDIDATES #######################
#####################################################################################

echo "Step 1b: Filtering out low matches, same species and likely false positives (due to naming errors on NCBI)"


##filter the results to only contain those above the containment threshold (default 0.5) (sometimes the description of the genome has commas...but can only get a csv output...)
cat 1.database_search/${prefix}.sourmash_multisearch.csv | awk -F "," -v containment="$containment" '{if(NR == 1) {print} else if($(NF-6) >= containment) {print $0}}' > 1.database_search/${prefix}.sourmash_multisearch.containment_filt.csv

##filter to only contain those that are from another species as labelled by NCBI
##create a new reduced header for this file and only print out neccesary columns (will also simplify the match name)
echo "element;match_assembly;match_species;containment" | tr ';' '\t' > 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv
tail -n+2 1.database_search/${prefix}.sourmash_multisearch.containment_filt.csv | while read line
do
element=$( echo "${line}" | awk -F "," '{print $1}' )
sample=$( echo "${element}" | awk -F "${separator}" '{print $1}' )
species=$( cat ${metadatapath} | awk -F "\t" -v sample="$sample" '{if(sample == $1) print $2}' )
echo "${line}" | awk -F "," -v species="$species" '{if($3 !~ species) print $1" "$3" "$(NF-6)}' | sed 's/"//g' | awk -F " " '{print $1"\t"$2"\t"$3"_"$4"\t"$NF}'
done >> 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv


##filter out likely false positives (a little tricky but can be very conservative in this filter)
##count all samples with more than 5 elements and then filter if 75% of all those elements are found in a single HGT candidate genome

##first just create a list of accession which and which species based on the metadata they likely resemble
echo "accession;likely_species" | tr ';' '\t' > 1.database_search/incorrectly_named_accession.tsv
tail -n+2 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv | awk -F "\t" '{print $1}' | sort -u | awk -F "$separator" '{print $1}' | uniq -c | awk '{if($1 >=5) print $2}' | while read sample
do
count=$( cat 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv | awk -F "_" -v sample="$sample" '{if($1 == sample) print}' | awk -F "\t" '{print $1}' | sort -u | wc -l )
falsep=$( cat 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv | awk -F "_" -v sample="$sample" '{if($1 ==sample) print}' | awk -F "\t" '{print $1"\t"$2}' | sort -u | awk -F "\t" '{print $2}' | sort | uniq -c | awk -v count="$count" '{if($1 >= (0.75*count)) print $2}' | sort -u )
species=$( cat ${metadatapath} | awk -F "\t" -v sample="$sample" '{if(sample == $1) print $2}' )
if [[ $falsep != "" ]]
then
echo "${falsep}" | awk -v species="$species" '{print $1"\t"species}'
fi
done | sort -u >> 1.database_search/incorrectly_named_accession.tsv

##then go through the previous candidate list and remove any other samples with the same species id from the metadata file and matches with that likely incorrectly labelled accession
echo "element;match_assembly;match_species;containment" | tr ';' '\t' > 1.database_search/${prefix}.sourmash_multisearch.candidates.final.tsv

##locate the lines associated with likely false positives
tail -n+2 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv | while read line2
do
sample=$( echo "${line2}" | awk -F "$separator" '{print $1}' )
species=$( cat ${metadatapath} | awk -F "\t" -v sample="$sample" '{if(sample == $1) print $2}' )
tail -n+2 1.database_search/incorrectly_named_accession.tsv | while read line
do
falsep=$( echo "${line}" | awk -F "\t" '{print $1}' )
altspecies=$( echo "${line}" | awk -F "\t" '{print $2}' )
echo "${line2}" | awk -F "\t" -v falsep="$falsep" -v species="$species" -v altspecies="$altspecies" '{if($2 == falsep && altspecies == species){print}}'
done
done > 1.database_search/toremove.tsv

##subtract the lines flagged as resulting from likely incorrect naming on ncbi
grep -v -F -x -f 1.database_search/toremove.tsv 1.database_search/${prefix}.sourmash_multisearch.candidates.tsv > 1.database_search/${prefix}.sourmash_multisearch.candidates.final.tsv

rm 1.database_search/toremove.tsv




#####################################################################################
########################## STEP 2. ANALYSE HGT CANDIDATES ###########################
#####################################################################################

minsize="10000"
minidentity="80"


##now we have candidates based just on kmer similarity however we now want to go through and look at whether alignments suggest transfer also

echo "Step 2a: Getting sourmash based candidates and generating actual alignments"

mkdir 2.HGT_candidates
mkdir 2.HGT_candidates/alignments

##create a bed file containing all the annotations given
echo "contig;start;end;sense;name;label" | tr ';' '\t' > 2.HGT_candidates/genes.bed
cat ${gff3path} | awk -F "\t" '{if($3 == "mRNA") print $1"\t"$4"\t"$5"\t"$7"\t"$9}' | awk -F ";Name=" '{print $0"\t"$NF}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' | awk -v identifier="$identifier"  '{if($5 ~ "_"identifier) {print $0"\ttyrR"} else if($5 ~ "_myb") {print $0"\tmyb"} else if($5 ~ "_duf3723") {print $0"\tduf3723"} else {print $0"\tNA"}}' >> 2.HGT_candidates/genes.bed

##loop through candidate files for each alignment etc
tail -n+2 1.database_search/${prefix}.sourmash_multisearch.candidates.final.tsv  | cut -f1 | sort -u | while read element
do

mkdir 2.HGT_candidates/alignments/${element}/

##get some metadata
sample=$( echo "${element}" | awk -F "${separator}" '{print $1}' )
species=$( cat ${metadatapath} | awk -F "\t" -v sample="$sample" '{if(sample == $1) print $2}' )


##first extract the contig containing the element to be used for alignment
coords=$( cat ${elementsbedpath} | awk -v element="$element" '{if($1 == element) print $2"\t"$3"\t"$4}'  )
coords2=$( echo "${coords}" | awk '{print $1":"$2"-"$3}' )
contig=$( cat ${elementsbedpath} | awk -v element="$element" '{if($1 == element) print $2}'  )
samtools faidx ${assembliespath} ${contig} > 2.HGT_candidates/alignments/${element}/${element}.contig.fa
samtools faidx ${assembliespath} ${coords2} > 2.HGT_candidates/alignments/${element}/${element}.element.fa


##create bed file for the element position
echo "contig;start;end" | tr ';' '\t' > 2.HGT_candidates/alignments/${element}/${element}.element.bed
echo "${coords}" >> 2.HGT_candidates/alignments/${element}/${element}.element.bed

##loop through the assemblies for each candidate; download them, create a single concantenated fasta file with mild renaming
tail -n+2 1.database_search/${prefix}.sourmash_multisearch.candidates.final.tsv | awk -v element="$element" '{if($1 == element) {print $2}}' | sort -u | while read candidategenome
do
candidategenome2=$( echo $candidategenome | sed 's/_//g' | awk -F "." '{print $1}' )

datasets download genome accession ${candidategenome} --no-progressbar
unzip -qq ncbi_dataset.zip 
rm ncbi_dataset.zip
ls ncbi_dataset/data/ | grep -v json | while read genome
do
cat ncbi_dataset/data/$genome/$genome*.fna | awk -F " " -v genome="$genome" '{if($1 ~ ">") {print ">"genome"_XXX"$1} else {print}}' | sed 's/_XXX>/_/g' 
done > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa
rm -r ncbi_dataset
rm README.md
rm md5sum.txt

##now using nucmer to align the contigs (will be used for plotting)
##and also use these alignments define the contigs/regions of interest to plot
nucmer -t ${threads} --delta 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.delta 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa 2.HGT_candidates/alignments/${element}/${element}.element.fa
##filter alignments for just the best element alignments in the genome
##also filter for only alignment above 80% identity and 10kb
delta-filter -q -i ${minidentity} -l ${minsize} 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.delta > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt
##convert to paf format and remove self matches
##use this paftools format to find the regions of interest for plotting
paftools.js delta2paf 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt | awk -F "\t" '{if($1 != $6) print}' > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt.paf


##now check that the alignment file isn't empty and and continue if so (i.e. if there is any good alignments to work with)
if [ -s "2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt.paf" ]
then

## create a bed file of the regions aligned wanted to be visualised 
## this includes just the starship and gap contig +- flanking regions
## also includes just the aligned contigs with using the max and min regions with >10kb alignment in the paf +- flanking regions
echo "contig;start;end;coords;tag" | tr ';' '\t' > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.aligned_regions.bed
##add starship coordinates (plus the species name and starship name)
##first have to calculate the contig length to make sure the flank extension doesn't exceed the contig size (would plot erroneous extensions)
maxlength=$( echo "${contig}" | seqkit grep --quiet -f - 2.HGT_candidates/alignments/${element}/${element}.contig.fa | grep -v ">" | tr -d '\n' | wc -c )
echo "${coords}" | awk -v flank="$flank" -v element="$element"  '{print $1"\t"$2-flank"\t"$3+flank"\t"element}' | awk '{if($2 < 0) {print $1"\t0\t"$3"\t"$1":"$2"-"$3"\t"$4} else {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}}' | awk -v maxlength="$maxlength" '{if($3 > maxlength) {print $1"\t"$2"\t"maxlength"\t"$1":"$2"-"maxlength"\t"$5} else {print}}'  >> 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.aligned_regions.bed
##add the aligned contigs coordinates (plus the flanking buffer)
cat 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt.paf | awk -F "\t" '{print $6}' | sort -u | while read candidatecontig
do
## calculate length of contig so that the flanks are not over extended
maxlength=$( echo "${candidatecontig}" | seqkit grep --quiet -f - 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa | grep -v ">" | tr -d '\n' | wc -c  )
##find the edges of the alignments
edges=$( cat 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.deltafilt.paf | awk -F "\t" -v candidatecontig="$candidatecontig" '{if($6 == candidatecontig) print}' | awk -F "\t" -v flank="$flank" 'BEGIN{max=0; min=99999999999999} {if($9 > max) {max=$9}; if($8 < min) {min=$8}} END{print min-flank"\t"max+flank}' | awk '{if($1 < 0) {print 0"\t"$2} else {print}}' | awk -F "\t" -v maxlength="$maxlength" '{if($2 > maxlength) {print $1"\t"maxlength} else {print}}' | sed 's/99999999949999/0/g'  )
##get the species to be used as a tag in the plot
species=$( cat 1.database_search/${prefix}.sourmash_multisearch.candidates.final.tsv  | awk -F "\t" -v candidategenome="$candidategenome" '{if($2 == candidategenome) print $3}' | sort -u )
echo "${candidatecontig};${edges};${species}" | tr ';' '\t' | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}'
done >> 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.aligned_regions.bed


##now align the full contigs all against all and use this for the plotting
nucmer -t ${threads} --delta 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.contigs.delta 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa 2.HGT_candidates/alignments/${element}/${element}.contig.fa
##this this alignment for the actual plotting though (pre-filtering so we can see complex regions etc)
paftools.js delta2paf 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.contigs.delta | awk -F "\t" '{if($1 != $6) print}' > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.contigs.delta.paf



##remove the full genome (no need for it now)
rm 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa

##plot the alignments using gggenomes in R
Rscriptpath=$( which gggenomes_skeleton.cargobay.R )
cat ${Rscriptpath} | sed "s/ELEMENT/${element}/g" | sed "s/CANDIDATEGENOME2/${candidategenome2}/g" | sed "s|PATHTOOUTPUT|${outputpath}|g" > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.R
Rscript --no-echo 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.R

else

touch 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.no_alignment.txt
echo "Warning: Although kmer similarities were found; there was no good alignment between ${element} and ${candidategenome2}"
 
fi




done

done


echo "Step 2b: Running BLASTall analyses and plotting alongside alignments"

###NEED TRANSCRIPTOME FOR THIS!!!!!

blastn -dust no -max_target_seqs 1 -max_hsps 1 -target 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa -query ${TRANSCRIPTOME} -outfmt 6 | awk '{print $1"\t"$3"\t"$4}' > 2.HGT_candidates/alignments/${element}/${element}.${candidategenome2}.fa


####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################


################## CREATE BLAST-ALL PLOT DATA ##################

##have to reset a variable with the starship genome to use to find the genes in the 'transcriptome' folder
oggenome=$( echo ${starship} | awk -F "_" '{print $1}' )
##search using all the genes in the blastdb 
## taking only the best hit (-max_target_seqs 1 -max_hsps 1)
##and performing no filtering for repeat regions (-dust no)
blastn -dust no -max_target_seqs 1 -max_hsps 1 -db alignments/${starship}/${starship}.${genome2}.fa -query ../transcriptome/${oggenome}.fa -outfmt 6 | awk '{print $1"\t"$3"\t"$4}' > alignments/${starship}/${starship}.${genome2}.transcriptome_blastn.tsv

conda deactivate
##get list of genes inside the starship to label them as such
start=$( cat ../${dataset}.starships_subtracted.nonredundant_starships_SLRs.UPDATE.tsv | awk -v starship="$starship" '{if($1 == starship) print $4}' )
end=$( cat ../${dataset}.starships_subtracted.nonredundant_starships_SLRs.UPDATE.tsv | awk -v starship="$starship" '{if($1 == starship) print $5}' )
cat ../magory.tyr.mod.consolidated.gff | grep ^"${oggenome}_" | awk -F "\t" -v starshipcontig="$starshipcontig" -v start="$start" -v end="$end" '{if($1==starshipcontig && $4 > start && $5 < end) {print "starship;"$9} else {print "NA;"$9}}' | awk -F ";" '{print $(NF-1)"\t"$1}' | sed 's/Name=//g' > alignments/${starship}/${starship}.${genome2}.transcriptome_starship_association.tsv
##now take this list and grab the blastn data if it exists, if not just mark it down as 0
echo "gene;identity;length;class" | tr ';' '\t' > alignments/${starship}/${starship}.${genome2}.transcriptome_blastn.class.tsv
cat alignments/${starship}/${starship}.${genome2}.transcriptome_starship_association.tsv | while read line
do
gene=$( echo "${line}" | awk '{print $1}' )
class=$( echo "${line}" | awk '{print $2}' )
cat alignments/${starship}/${starship}.${genome2}.transcriptome_blastn.tsv | awk -v gene="$gene" -v class="$class" -v oggenome="$oggenome"  'BEGIN{id="0" ; len="0"} {if(oggenome"_"$1 == gene) {id=$2;len=$3}} END{print gene"\t"id"\t"len"\t"class}'
done >> alignments/${starship}/${starship}.${genome2}.transcriptome_blastn.class.tsv


######################################################


##remove genome now that contig is extracted
rm alignments/${starship}/${starship}.${genome2}.fa
##align all vs all and each starship vs all contigs for each genome (extracted from the all vs all alignment)
conda activate mummer4
nucmer --maxmatch --delta alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.nucmer.delta alignments/${starship}/${starship}.${genome2}.aligned_contigs.fa alignments/${starship}/${starship}.${genome2}.aligned_contigs.fa
conda deactivate
##convert to paf format and remove self matches
paftools.js delta2paf alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.nucmer.delta | awk -F "\t" '{if($1 != $6) print}' > alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.nucmer.paf


## 3. a bed file of the regions aligned wanted to be visualised 
## this includes just the starship and gap contig +- flanking regions
## also includes just the aligned contigs with using the max and min regions with >10kb alignment in the paf +- flanking regions
echo "contig;start;end;coords;species" | tr ';' '\t' > alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.bed
##add starship coordinates (plus the species name and starship name)
cat ../${dataset}.starships_subtracted.nonredundant_starships_SLRs.UPDATE.tsv | awk -v starship="$starship" '{if($1 == starship) print}' | cut -f3-5 | awk '{print $1"\t"$1"\t"$2"\t"$3}'  | awk -v flank="$flank" -v starship="$starship"  '{print $2"\t"$3-flank"\t"$4+flank"\tP."$1":"starship}' | awk '{if($2 < 0) {print $1"\t0\t"$3"\t"$1":"$2"-"$3"\t"$4} else {print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}}' >> alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.bed
##add the aligned contigs coordinates (plus the flanking buffer)
cat alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.nucmer.paf | awk '{if($11 > 7500) print}' | cut -f1 | sort -u | grep -v ${starshipcontig} | while read contig
do
## calculate length of contig so that the flanks are not over extended
maxlength=$( echo "${contig}" | seqkit grep --quiet -f - alignments/${starship}/${starship}.${genome2}.aligned_contigs.fa | grep -v ">" | tr -d '\n' | wc -c  )
edges=$( cat alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.nucmer.paf | awk '{if($11 > 7500) print}' | grep ^"${contig}" | awk -v flank="$flank" 'BEGIN{max=0; min=99999999999999} {if($4 > max) {max=$4}; if($3 < min) {min=$3}} END{print min-flank"\t"max+flank}' | awk '{if($1 < 0) {print 0"\t"$2} else {print}}' | awk -F "\t" -v maxlength="$maxlength" '{if($2 > maxlength) {print $1"\t"maxlength} else {print}}' | sed 's/99999999949999/0/g'  )
species=$( cat starships.pezizomycotina_ncbi.blastn.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv  | awk -F "\t" -v contig="$contig" '{if($2 == contig) print $9}' | sort -u )
echo "${contig};${edges};${species}" | tr ';' '\t' | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}'
done >> alignments/${starship}/${starship}.contig.${genome2}.aligned_contigs.bed


## now generate per combination a specific R file for generating plots using gggenomes
## we will use the alignment_plotting.template.R in scratch/saodonnell/projects/${dataset}/starfish/pezizomycotina_BLAST
## to do so we just need to change a few terms, i.e. the file names for each of the four files used as input
cp /scratch/saodonnell/projects/genomegraphs/${dataset}/alignment_plotting.template.2.R alignments/${starship}/${starship}.${genome2}.R
## first swap in all the positions for the starship name
sed -i "s/STARSHIP/${starship}/g" alignments/${starship}/${starship}.${genome2}.R
##now the genome
sed -i "s/GENOME/${genome2}/" alignments/${starship}/${starship}.${genome2}.R


rm alignments/${starship}/${starship}.${genome2}.aligned_contigs.fa
done

gzip alignments/${starship}/${starship}.contig.fa

done


####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################################







if [[ $cleanup == "yes" ]]
then
##remove sourmash database (as it usually takes up aroud 5Gb)
rm -r 0.database/genbank-20250408-fungi-k21.zip

fi
