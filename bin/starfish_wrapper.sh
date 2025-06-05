#!/bin/bash
set -euo pipefail

version="v1"

##############################################################
################ STEP 0a: SETTING VARIABLES ##################
##############################################################

#default values, unless denoted when running stargraph
assemblies=""
gff3=""
threads="1"
flank="50000"
prefix="starfish"
separator="_"
output="starfish_output"
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
    -g|--gff3)
    gff3="$2"
    shift
    shift
    ;;    
    -t|--threads)
    threads="$2"
    shift
    shift
    ;;
    -f|--flank)
    flank="$2"
    shift
    shift
    ;;
    -s|--separator)
    separator="$2"
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
 
    starfish_wrapper.sh -a assemblies_panSN.txt
    
    Required inputs:
    -a | --assemblies       A txt file with each line containing the path to an assembly using the PanSN-spec naming scheme for each contig

    Recommended inputs:
    -g | --gff3             A gff3 file from the same assemblies with the same panSN-spec contig namings (contig names need to correspond)
    -t | --threads          Number of threads for tools that accept this option (default: 1)

    starfish specific inputs:
    -f | --flank            Size of flanking regions either side of element used during elementViz (Default: 50000)

    Optional parameters:
    -s | --separator        PanSN-spec naming separator used (Default: _)
    -p | --prefix           Prefix for output (Default: starfish)
    -o | --output           Name of output folder for all results (Default: starfish_output)
    -c | --cleanup          Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
    -h | --help         Print this help message

    "
    exit
    ;;
    esac
done


#creates error message and exits if these values are not/incorrectly assigned 
[[ $assemblies == "" ]] && echo "ERROR: No file containing assembly paths was given" && exit

[[ $gff3 == "" ]] && echo "WARNING: Running without any annotation data provided in gff3 format (-g)"
[[ $threads == "1" ]] && echo "WARNING: Running with default thread count; depending on the dataset this may take some time"
[[ $separator == "_" ]] && echo "WARNING: Running with default separator '_' (this is recommended but just a heads up)"


##############################################################
####################### 0b. SETTING UP #######################
##############################################################

##redefining some variables and checking
##get kb version of flank size for file output naming
flank2=$( echo ${flank} | awk '{print $0/1000}' )
##get twice the flank length to use to plotting
flank3=$( echo ${flank} | awk '{print $0*2}' )

##paths to raw data
assembliespath=$( realpath ${assemblies} )
if [[ ${gff3} != "" ]]
then
gff3path=$( realpath ${gff3} )
fi

#check if the files given actually exist
[ ! -f "${assembliespath}" ] && echo "ERROR: Cannot find path to assemblies file provided by -a; check path is correct and file exists" && exit
if [[ ${gff3} != "" ]]
then
[ ! -f "${gff3path}" ] && echo "ERROR: Cannot find path to gff3 file provided by -g; check path is correct and file exists" && exit
fi

##check if output directory already exists
[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

##create output directory
mkdir ${output}
##get the absolute paths to each of the assemblies and save in a new file
cat ${assembliespath} | while read genome
do
realpath ${genome}
done > ${output}/path_to_assemblies.txt

##do the same for the gff3 files is given
if [[ ${gff3} != "" ]]
then
cat ${gff3path} | while read gff3
do
realpath ${gff3}
done > ${output}/path_to_gff3.txt
fi

##move into output directory
cd ${output}

##modify the gff3 files to make sure the sample can be tied to the corresponding genes (as required by starfish's scheme)
##we just add another tag to the end of the file with the name of the gene being accession"_"gene
##and extract only the 'genes' line
if [[ ${gff3} != "" ]]
then
mkdir gff3
cat path_to_gff3.txt | while read file
do
assembly=$( echo $file | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'  )
##########################################SOLVE THIS PART, NOT FINISHED
done
fi
#cat $file | while read line
#do
#gene=$( echo "$line" | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//g'   )
#echo "$line" | awk -v assembly="$assembly" -v gene="$gene" '{if($0 ~ "#") {print $0} else {print assembly"_"$0"Name="assembly"_"gene";"}}' | sed 's/scaffold_/scaffold/g' 
#done > gff3/${assembly}.gff3
#done


##need a list of the assemblies and paths (using the softmasked genomes previously annotated)
cat path_to_assemblies.txt | while read assembly
do
name=$( grep '>' ${assembly} | head -n1 | awk -F "${separator}" '{print $1}' | sed 's/>//' )
echo ${name}";"${assembly} | tr ';' '\t'
done > ome2assembly.txt

##list of annotation files and paths also
#realpath gff3/*.gff3 | while read line
#do
##SOLVE NAMING ISSUE HERE: WANT TO USE THE SAMPLE NAME IN THE GENOME FILE...SHOULD BE THE SAME AS IN THE GFF
##...MAYBE JUST IT DEFINITELY NEEDS TO BE THE SAME THEREFORE GRAB IT FROM THE GFF3 FIRST COLUMN
#name=$( echo $line | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
#echo ${name}";"${line} | tr ';' '\t'
#done > ome2gff.txt

## a consolidated file with all annotations to compare with Starfish
#cat gff3/*.gff3 > ${prefix}.assemblies.gff3

##################################################################
####################### STEP 1: geneFinder #######################
##################################################################


## now run the genefinder step with the tyrosine recombinases
mkdir geneFinder_tyr/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/YRsuperfams.p1-512.hmm -P $CONDA_PREFIX/db/YRsuperfamRefs.faa -i tyr -o geneFinder_tyr/
#starfish annotate -T ${threads} -x ${prefix} -a ome2assembly.txt -g ome2gff.txt -p $CONDA_PREFIX/db/YRsuperfams.p1-512.hmm -P $CONDA_PREFIX/db/YRsuperfamRefs.faa -i tyr -o geneFinder_tyr/

##run the same step for the other starship-related genes
mkdir geneFinder_duf3723/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/duf3723.hmm -P $CONDA_PREFIX/db/duf3723.mycoDB.faa -i duf3723 -o geneFinder_duf3723/
mkdir geneFinder_fre/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/fre.hmm -P $CONDA_PREFIX/db/fre.mycoDB.faa -i fre -o geneFinder_fre/
mkdir geneFinder_nlr/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/nlr.hmm -P $CONDA_PREFIX/db/nlr.mycoDB.faa -i nlr -o geneFinder_nlr/
mkdir geneFinder_plp/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/plp.hmm -P $CONDA_PREFIX/db/plp.mycoDB.faa -i plp -o geneFinder_plp/
mkdir geneFinder_myb/
starfish annotate -T ${threads} -s ${separator} -x ${prefix} -a ome2assembly.txt -p $CONDA_PREFIX/db/myb.hmm -P $CONDA_PREFIX/db/myb.SRG.fa -i myb -o geneFinder_myb/

##generate a combined gff from all the features (to be used by stargraph for Starship-like region classification)
cat geneFinder_*/${prefix}.filt.gff > ${prefix}.filt.SRGs_combined.gff


##get unique regions for each captain (cluster them based on distance)
##first consolidate starfish captain positions with annotations (should spit out a file called XXX_tyr.filt_intersect.consolidated.gff)
#starfish consolidate -o ./ -g ${prefix}.assemblies.gff3 -G geneFinder_tyr/${prefix}.filt_intersect.gff
##now generate a tsv with the path to the gff generated in the previous step
#realpath ${prefix}.filt_intersect.consolidated.gff | awk -v dataset="$dataset" '{print dataset"\t"$0}' > ome2consolidatedGFF.txt
##now run sketch to do the consolidation
#starfish sketch -m 10000 -q geneFinder_tyr/${prefix}.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x ${prefix} -o geneFinder_tyr/
realpath geneFinder_tyr/${prefix}.filt.gff | awk -v prefix="$prefix" '{print prefix"\t"$0}' > ome2GFF.txt
starfish sketch -s ${separator} -m 10000 -q geneFinder_tyr/${prefix}.filt.ids -g ome2GFF.txt -i s -x ${prefix} -o geneFinder_tyr/

##grab just the candidate captain coords
grep -P '\ttyr\t' geneFinder_tyr/${prefix}.bed > geneFinder_tyr/${prefix}.tyr.bed 


##################################################################
###################### STEP 2: elementFinder #####################
##################################################################


##now use the tyrosine elements to find starships
mkdir elementFinder
##generate a blastdb with the assemblies
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/${prefix}.assemblies.fa
makeblastdb -in blastdb/${prefix}.assemblies.fa -out blastdb/${prefix}.assemblies -parse_seqids -dbtype nucl
##look for the actual insertions
starfish insert -s ${separator} -T ${threads} -a ome2assembly.txt -d blastdb/${prefix}.assemblies -b geneFinder_tyr/${prefix}.tyr.bed -i tyr -x ${prefix} -o elementFinder/

##find flanking elements
#starfish flank -a ome2assembly.txt -b elementFinder/${prefix}.insert.bed -x ${prefix} -o elementFinder/

##summarise the info (not using the flank data)
#starfish summarize -a ome2assembly.txt -b elementFinder/${prefix}.insert.bed -x ${prefix} -o elementFinder/ -S elementFinder/${prefix}.insert.stats -g ome2consolidatedGFF.txt -t geneFinder_tyr/${prefix}.filt_intersect.ids 
starfish summarize -s ${separator} -a ome2assembly.txt -b elementFinder/${prefix}.insert.bed -x ${prefix} -o elementFinder/ -S elementFinder/${prefix}.insert.stats -g ome2GFF.txt -t geneFinder_tyr/${prefix}.filt.ids 


####################################################################################
###################### STEP 3a: CLASSIFY CAPTAINS AND ELEMENTS #####################
####################################################################################


## first assign the family using reference captain sequences HMM
hmmsearch --noali --notextw -E 0.001 --max --cpu ${threads} --tblout elementFinder/${prefix}_tyr_vs_YRsuperfams.out $CONDA_PREFIX/db/YRsuperfams.p1-512.hmm geneFinder_tyr/${prefix}.filt.fas
#hmmsearch --noali --notextw -E 0.001 --max --cpu ${threads} --tblout elementFinder/${prefix}_tyr_vs_YRsuperfams.out $CONDA_PREFIX/db/YRsuperfams.p1-512.hmm geneFinder_tyr/${prefix}.filt_intersect.fas
## modified the grep -v '#' part of generic pipeline modified to grep -v ^'#' as it is just trying to get rid of notes and get conflicted with the seperator here
perl -p -e 's/ +/\t/g' elementFinder/${prefix}_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v ^'#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/${prefix}_tyr_vs_YRsuperfams_besthits.txt

## rename the captains by their starships
grep -P '\tcap\t' elementFinder/${prefix}.elements.bed | cut -f4,7 > elementFinder/${prefix}.cap2ship.txt
$CONDA_PREFIX/aux/searchReplace.pl --strict -i elementFinder/${prefix}_tyr_vs_YRsuperfams_besthits.txt -r elementFinder/${prefix}.cap2ship.txt > elementFinder/${prefix}_elements_vs_YRsuperfams_besthits.txt

## group starships into naves using coverage threshold of 25% and identity of 50%
## make sure to use the filt_intersect.fas file for the putative captains as otherwise the name will not correspond downstream and naves information will be lost
mmseqs easy-cluster geneFinder_tyr/${prefix}.filt.fas elementFinder/${prefix}_tyr elementFinder/ --threads 1 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
$CONDA_PREFIX/aux/mmseqs2mclFormat.pl -i elementFinder/${prefix}_tyr_cluster.tsv -g navis -o elementFinder/

## group starships into haplotypes using kmer similarity (95% threshold?)
$CONDA_PREFIX/bin/starfish sim -m element -t nucl -b elementFinder/${prefix}.elements.bed -x ${prefix} -o elementFinder/ -a ome2assembly.txt
#$CONDA_PREFIX/bin/starfish group -m mcl -s elementFinder/${prefix}.element.nucl.sim -i hap -o elementFinder/ -t 0.05
$CONDA_PREFIX/bin/starfish group -s elementFinder/${prefix}.element.nucl.sim -i hap -o elementFinder/ -t 0.05


## replace captain IDs with starship IDs in the naves file
$CONDA_PREFIX/aux/searchReplace.pl -i elementFinder/${prefix}_tyr_cluster.mcl -r elementFinder/${prefix}.cap2ship.txt > elementFinder/${prefix}.element_cluster.mcl
## merge navis with haplotype to create a navis-haplotype label for each Starship
$CONDA_PREFIX/aux/mergeGroupfiles.pl -t elementFinder/${prefix}.element_cluster.mcl -q elementFinder/${prefix}.element.nucl.I1.5.mcl > elementFinder/${prefix}.element.navis-hap.mcl
## convert mcl to gene2og format to simplify downstream parsing:
awk '{ for (i = 2; i <= NF; i++) print $i"\t"$1 }' elementFinder/${prefix}.element.navis-hap.mcl > elementFinder/${prefix}.element.navis-hap.txt
## now add the family and navis-haplotype into to the element.feat file to consolidate metadata:
join -t$'\t' -1 1 -2 2 <(sort -t$'\t' -k1,1 elementFinder/${prefix}.element.navis-hap.txt | grep -P '_e|_s') <(sort -t$'\t' -k2,2 elementFinder/${prefix}.elements.feat) | awk -F'\t' '{print}' > elementFinder/${prefix}.elements.temp.feat
echo -e "#elementID\tfamilyID\tnavisHapID\tcontigID\tcaptainID\telementBegin\telementEnd\telementLength\tstrand\tboundaryType\temptySiteID\temptyContig\temptyBegin\temptyEnd\temptySeq\tupDR\tdownDR\tDRedit\tupTIR\tdownTIR\tTIRedit\tnestedInside\tcontainNested" > elementFinder/${prefix}.elements.ann.feat
join -t$'\t' -1 1 -2 1 <(sort -t$'\t' -k1,1 elementFinder/${prefix}_elements_vs_YRsuperfams_besthits.txt | grep -P '_e|_s' | cut -f1,2) <(sort -t$'\t' -k1,1 elementFinder/${prefix}.elements.temp.feat) | awk -F'\t' '{print}' >> elementFinder/${prefix}.elements.ann.feat



##################################
####not going to do the regionFinder step (NEED ORTHOGROUPS)

##now consider where the starships are in the genome
mkdir regionFinder

## create a file with tyrs that are not found in any elements (this will let us assign them to fragmented haplotypes in the dereplicate analysis, which can be helpful if annotating 'dead' or 'derelict' element copies):
#grep -f <(comm -23 <(cut -f1 geneFinder_tyr/${prefix}.filt_intersect.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/${prefix}.elements.bed | cut -f4| sort)) geneFinder_tyr/${prefix}.tyr.bed > regionFinder/unaffiliated_tyrs.bed
grep -f <(comm -23 <(cut -f1 geneFinder_tyr/${prefix}.filt.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/${prefix}.elements.bed | cut -f4| sort)) geneFinder_tyr/${prefix}.tyr.bed > regionFinder/unaffiliated_tyrs.bed

## filter for confident orthogroups
#$CONDA_PREFIX/aux/filterOG.pl -O orthofinder/orthofinder_results/Results_*/Orthogroups/Orthogroups.txt -a 1 -c 5 -o ./

## now, dereplicate your data to identify independently segregating element insertions
#starfish dereplicate -e elementFinder/${prefix}.element.navis-hap.mcl -t regionFinder/unaffiliated_tyrs.bed -F elementFinder/${prefix}.elements.feat -S elementFinder/${prefix}.elements.named.stats -O Orthogroups.a1.c6.txt -g ome2gff.txt -x ${prefix} -o regionFinder/ --flanking 6 --mismatching 1

## count the number of regions each navis-haplotype is found in (i.e. transposed element)
## modified the grep -v '#' part of generic pipeline modified to grep -v ^'#' as it is just trying to get rid of notes and get conflicted with the seperator here
#grep -v ^'#' regionFinder/${prefix}.fog3.d600000.m1.dereplicated.txt | cut -f2 | sort | uniq -c | perl -pe 's/ +//' | sort -k1,1nr


##################################DON'T REALLY LIKE LOCUS VIZ, SO NOT ADDING IT, PLUS NEED ORTHOGROUP INFO SO CAN'T DO IT IF NOT ANNOTATED


## now visualise the loci
#mkdir locusViz
## add gc content information to be used by the locus visualisation
$CONDA_PREFIX/aux/seq-gc.sh -Nbw 1000 blastdb/${prefix}.assemblies.fa > ${prefix}.assemblies.gcContent_w1000.bed
#rm blastdb/${prefix}.assemblies.fna

#starfish locus-viz -T 2 -m region -a ome2assembly.txt -b elementFinder/${prefix}.elements.bed -x ${prefix} -o locusViz/ -A nucmer -r regionFinder/${prefix}.fog3.d600000.m1.regions.txt -d regionFinder/${prefix}.fog3.d600000.m1.dereplicated.txt -j regionFinder/${prefix}.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder_tyr/${prefix}.filt_intersect.ids --gc ${prefix}.assemblies.gcContent_w1000.bed
#starfish locus-viz -T 2 -m region -a ome2assembly.txt -b elementFinder/${prefix}.elements.bed -x ${prefix} -o locusViz/ -A nucmer -r regionFinder/${prefix}.fog3.d600000.m1.regions.txt -d regionFinder/${prefix}.fog3.d600000.m1.dereplicated.txt -j regionFinder/${prefix}.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2gff.txt --tags geneFinder_tyr/${prefix}.filt.ids --gc ${prefix}.assemblies.gcContent_w1000.bed

##now use the R scripts that were made and the R environment with gggenomes to generate the plots
#ls locusViz/*.R | while read script
#do
#region=$( echo $script | awk -F "/" '{print $NF}' | awk -F "." '{print $2}' )
###modify the Rscript so that the legend is printed at the top and adjusting the width by adding twice the flanking length for each edge
#sed 's/geom_gene_tag/theme(legend.position = "top")+geom_gene_tag/g' $script | sed "s/regionSeqs\$length/regionSeqs\$length+${flank3}/g" > locusViz/${prefix}.${element}.mod.R
#Rscript locusViz/${prefix}.${element}.mod.R
#done

############################################



##########################################################################
###################### STEP 3b: VISUALISE INSERTIONS #####################
##########################################################################


##visualise 
mkdir pairViz
starfish pair-viz -m all -t empty -T ${threads} -A nucmer -a ome2assembly.txt -b elementFinder/${prefix}.elements.bed -S elementFinder/${prefix}.elements.named.stats -o pairViz/


###############################################################################
###################### STEP 3c: VISUALISE NAVIS GROUPINGS #####################
###############################################################################

##now we have all the locations and their starships, but what about visualising the elements and all their regions
##pick the flanking size in kb
mkdir elementViz_${flank2}kbflank

##first go through all the naves and identify which ones have multiple elements
##create a list of those and then loop through that list to create the locus-viz element plots

cat elementFinder/${prefix}.elements.ann.feat | cut -f2,3 | awk -F "-" '{print $1}' | sort -u | while read set
do
set2=$( echo ${set} | awk '{print $1"-"$2}'  )

##now get the haplotypes associated with each element and sort the elements using this (doesn't change anything...?)
elementcount=$( grep "${set}" elementFinder/${prefix}.elements.ann.feat | wc -l )
if [[ ${elementcount} -gt 1 ]]
then
echo "${set}"
fi

done > elementViz_${flank2}kbflank/naves_list.txt


##run the element visualisation by navis grouping (to see if haplotypes are truely different etc)
##to do this the list of elements given is the order of the plot, therefore sort by haplotype and locusViz region
cat elementViz_${flank2}kbflank/naves_list.txt | while read set
do
##new name for the set
set2=$( echo "${set}" | awk '{print $1"-"$2}'  )

##get the list of element names per family-navis clade
#grep "${set}" elementFinder/${prefix}.elements.ann.feat | cut -f1 | while read element
#do
##now get the regions associated with each element and sort by the region (here we lose haplotype information in the sorting...)
#region=$( grep ${element} regionFinder/${prefix}.fog3.d600000.m1.dereplicated.txt | cut -f1 )
#echo ${element}";"${region} | tr ';' '\t' | sort -k2 
#done > elementViz_${flank2}kbflank/${set2}.list

##now get the haplotypes associated with each element and sort the elements using this (doesn't change anything...?)
grep "${set}" elementFinder/${prefix}.elements.ann.feat | cut -f1 | while read element
do
haplotype=$( grep ${element} elementFinder/${prefix}.element.navis-hap.mcl | awk '{print $1}' )
echo ${element}";"${haplotype} | tr ';' '\t' | sort -k2 
done > elementViz_${flank2}kbflank/${set2}.list

##then use locuz-viz with in 'element' mode
#starfish locus-viz -T 2 -m element -a ome2assembly.txt -b elementFinder/${prefix}.elements.bed -x ${set2} -U ${flank} -D ${flank} -l elementViz_${flank2}kbflank/${set2}.list  -o elementViz_${flank}kbflank/ -A nucmer -r regionFinder/${prefix}.fog3.d600000.m1.regions.txt -d regionFinder/${prefix}.fog3.d600000.m1.dereplicated.txt -j regionFinder/${prefix}.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder_tyr/${prefix}.filt_intersect.ids --gc ${prefix}.assemblies.gcContent_w1000.bed

starfish locus-viz -T 2 -m element -a ome2assembly.txt -b elementFinder/${prefix}.elements.bed -x ${set2} -U ${flank} -D ${flank} -l elementViz_${flank2}kbflank/${set2}.list  -o elementViz_${flank2}kbflank/ -A nucmer -g ome2GFF.txt --tags geneFinder_tyr/${prefix}.filt.ids --gc ${prefix}.assemblies.gcContent_w1000.bed

done


####PROBABLY NEEDS TO BE A PART OF STARGRAPH AND IN ANOTHER SCRIPT.

##now use the R scripts that were made and the R environment with gggenomes to generate the plots
ls elementViz_${flank2}kbflank/*.R | while read script
do
element=$( echo $script | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
##modify the Rscript so that the legend is printed at the top and adjusting the width by adding twice the flanking length for each edge
sed 's/geom_gene_tag/theme(legend.position = "top")+geom_gene_tag/g' $script | sed "s/regionSeqs\$length/regionSeqs\$length+${flank3}/g" > elementViz_${flank2}kbflank/${element}.mod.R
Rscript elementViz_${flank2}kbflank/${element}.mod.R
done


echo "################################################################"
echo "Need to now manually validate the candidate starships; have fun"
echo "################################################################"


