<p align="center" >
    <img src="https://github.com/SAMtoBAM/stargraph/blob/main/logo/stargraph_logo.png" width=100%>
</p>

[![Zenodo DOI](https://zenodo.org/badge/992224579.svg)](https://zenodo.org/badge/)
[![Anaconda_version](https://anaconda.org/samtobam/stargraph/badges/version.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_platforms](https://anaconda.org/samtobam/stargraph/badges/platforms.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_downloads](https://anaconda.org/samtobam/stargraph/badges/downloads.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/stargraph/badges/latest_release_date.svg)](https://anaconda.org/samtobam/stargraph)

**_stargraph_** is a tool that detects _Starship_-like regions (SLRs) using a genome-graph based approach and combines this with _Starship_ results from the more conservative tool [starfish](https://github.com/egluckthaler/starfish) <br/>
The combination of both tools provides a comprehensive view of genomic regions impacted by _Starships_ <br/>
**_stargraph_** requires a minimum of two contiguous long-read assemblies in order to identify these _Starship_ and _Starship_-like regions

**_stargraph_** requires the same input as ```starfish``` and some ```starfish``` output <br/>
As ```stargraph``` requies ```starfish``` input; The ```starfish``` pipeline should be run first in order to feed ```stargraph``` with both tyrosine recombinase (TyrR) and _Starship_ positions. See Step 0.

## Pipeline
Using the tool requires 7 steps: <br/>
-1 --> 0 Preprocessing/Set up: <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp;-1: Preprocessing of input assembly data <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 0: Running ```starfish``` <br/>
&nbsp;1 --> 5 ```stargraph```: <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 1: Generating a genome-graph <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 2: Identifying Presence/Absence Variants (PAVs) <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 3: Elevating PAVs to _Starship_-like regions, identifying 'haplotypes' and plotting insertion sites <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 4: Combining SLRs with _Starships_ to generate a non-redundant dataset <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 5: Generating alignments and Network analyses <br/>
&nbsp;6 --> 7 additional stargraph modules: <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 6: ```COMING SOON; allstars``` to classify newly found elements by similarity to a database of described/named elements  <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 7: ```cargobay``` use a database of public fungal assemblies to find evidence of HGT  <br/>

 
## Apptainer usage
```
docker pull ghcr.io/samtobam/stargraph:latest
```

## Conda installation
```
conda install samtobam::stargraph
```

# STEP -1. Preprocessing of input assembly data

To help with building genome-graphs (and starfish compatability) it helps if sample information is stored in the header of each contig for each assembly <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp;	[sample_name][delim][contig/scaffold_name] <br/>
e.g. The fasta header for a contig called 'CP097570.1' from a sample/strain called CEA10 using the underscore as a delimiter: <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp;	>CEA10_CP097570.1  <br/>
The '\_' underscore is highly recommended due to compatability with starfish/other tools and lower frequency in names/use than other separators. <br/>
And just to be sure; the separator (use the underscore please) cannot be in the sample or contig/scaffold name <br/>

This [PanSN-spec-like](https://github.com/pangenome/PanSN-spec) naming modification needs to be done for all assemblies in your dataset <br/>
For assemblies directly downloaded from NCBI, the trailing information (e.g. 'Aspergillus fumigatus CEA10 chromosome 8') after the contig accession can be left as is and will be ignored (due to the space seperation from the contig name) <br/>

Following this you need to create a txt file (e.g.  _assemblies_panSN.txt_) containing one path per line to each of the PanSN-spec-like renamed assemblies <br/>
And voila, the primary input required for ```stargraph``` is ready. <br/>
Feed this _assemblies_panSN.txt_ file to ```stargraph```; input parameter ```-a | --assemblies```.


# STEP 0. Running starfish first (wrapper included)

```stargraph``` requires some ```starfish``` input in order to run in its entirety <br/>
This includes: <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; 1. The _de-novo_ annotations of Tyrosine recombinases used to elevate PAVs to SLRs <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; (usually can use: _starfish_output/geneFinder/\*.filt.gff_ or _starfish_output/${prefix}.filt.SRGs_combined.gff_ ) <br/> 
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; (```stargraph``` input parameter ```-r | --tyrRs```) <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; 2. A final list of curated _Starship_ elements (combined with SLRs to generate the non-redundant dataset) <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; (usually can use: _'starfish_output/elementFinder/\*.elements.ann.feat'_) <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; (```stargraph``` input parameter ```-e | --elements```) <br/>

Therefore ```starfish``` needs to be run first (installed in stargraph environment) <br/>
You can follow the starfish tutorials running; provided on the [github/wiki](https://github.com/egluckthaler/starfish) <br/>
Or <br/>
To simplify running ```starfish``` and ensure compatability etc: you can use the wrapper ```starfish_wrapper.sh``` provided by ```stargraph``` <br/>
The wrapper runs the primary steps required with most default parameters <br/>
In this case the input used will be the same list of paths to the PanSN-spec-like renamed assemblies as used for ```stargraph``` e.g. _assemblies_panSN.txt_ <br/>
_Note: You can use all available assemblies for starfish but only long-contiguous assemblies for stargraph_

	starfish_wrapper.sh -a assemblies_panSN.txt

The final wrapper output includes a set of putative _Starships_.

Additional steps in the wrapper include the _de-novo_ detection of DUF3723 and MYB/SANT genes associated with _Starships_ to be used in _Starship_ and SLR visualisation <br/>
These annotations are combined in the output file _starfish_output/${prefix}.filt.SRGs_combined.gff_ (best input for stargraph --tyrRs ) <br/>
Due to the putative assocaition of these MYB/SANT elements near the opposite _Starship_ edge to the captain; they are particularly helpful in visualising the ends of some elements

# STEP 1-5 Running stargraph

**_stargraphs_** initial module will find _Starship_-like regions and combine them with _Starships_ <br/>
To do this combine the list of assembly paths with ```starfish``` output: <br/>
the tyrosine recombinase annotations (_starfish_output/geneFinder/\*.filt.gff_ or _starfish_output/${prefix}.filt.SRGs_combined.gff_  | --tyrRs) <br/>
the _Starship_ annotations (_starfish_output/elementFinder/\*.elements.ann.feat_ | --elements)

 
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
	-w | --window			Size of windows used for PAV detection (Default: 1000)
	-m | --minsize			Minimum size of PAVs to be kept (Default: 30000)
	-x | --maxsize			Maximum size of SLRs to be kept; filter only applied after starship merging (Default: 2000000)
	-k | --kmerthreshold	The minimum 'max-containment/jaccard-similarity' value to be used for clustering of elements for visualiation (Default: 0.3)
	-f | --flank			Size of flanking region used when plotting element alignments (Default: 75000)
	-p | --prefix			Prefix for output (Default: stargraph)
	-o | --output			Name of output folder for all results (Default: stargraph_output)
	-c | --cleanup			Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help			Print this help message



The final output **_stargraph_output/${prefix}.starships_SLRs.tsv_** contains the final results of ```stargraph```; a nonredundant list of _Starships_ and _Starship_-like elements <br/>
Additionally: <br/>
**1.*.pggb** : all genome-graph output from both pggb and odgi <br/>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; information about the identified regions of Presence/Absence variation <br/>
**2.PAVs_to_SLRs** : information on the elevation of PAVs to SLRs using the provided tyrosine recombinase locations <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; contains information on the elevation of PAVs not elevated to SLR status but containing DUF3723 or MYB genes <br/>
**3.SLR_plots** : plots showing the alignment of SLRs clustered based on k-mer max-containment, including one insertion site per cluster <br/>
**4.SLR_starship_combination** : information on generating the non-redundant dataset combining the provided _Starships_ with the newly identified SLRs <br/>
**5.SLR_starship_network_alignments** : plots showing the alignment of _Starships_ and SLRs clustered together based on k-mer max-containment <br/>
&nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; plots of networks using both Jaccard similarity and containment <br/>


# COMING SOON STEP 6 Running allstars 

Use the allstars module in order to classify your elements using a manually curated database of named elements (_db/named_starships_database.curated.fa_)

	allstars.sh -e stargraph_output/${prefix}.starships_SLRs.fa -l $CONDA_PREFIX/db/named_starships_database.curated.fa

# STEP 7 Running cargobay

Use a database of all public fungal assemblies on NCBI (thank you sourmash team!) in order to look for your elements in other species 

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
