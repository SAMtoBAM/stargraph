# stargraph
A genome-graph based starship detection plugin to be combined with starfish

<p align="center" >
    <img src="https://github.com/SAMtoBAM/stargraph/blob/main/logo/stargraph.png" width=100%>
</p>

[![Zenodo DOI](###.svg)](###)
[![Anaconda_version](https://anaconda.org/samtobam/stargraph/badges/version.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_platforms](https://anaconda.org/samtobam/stargraph/badges/platforms.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_downloads](https://anaconda.org/samtobam/stargraph/badges/downloads.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/stargraph/badges/latest_release_date.svg)](https://anaconda.org/samtobam/stargraph)

**_stargraph_** is a tool that detects _Starship_-like regions (SLRs) using a genome-graph based and combines this with results from the more conservative tool ```starfish``` <br/>
The combination of both tools provides a more comprehensive view of genomic regions impacted by _Starships_

**_stargraph_** requires only the same input as ```starfish``` and some ```starfish``` output and therefore can be used easily in conjunction <br/>
As ```stargraph``` requies ```starfish``` input; The ```starfish``` pipeline should be run first in order to feed ```stargraph``` with both tyrosine recombinase (TyrR) and _Starship_ positions.

# Easy installation

	conda install samtobam::stargraph


# Preprocessing of input assembly data

To help with building genome-graphs each contig for each assembly needs sample/haplotype information stored in the header <br/>
To do this we use the PanSN specifications (explained in detail [here](https://github.com/pangenome/PanSN-spec)): <br/>
&nbsp; &nbsp; &nbsp; &nbsp;	[sample_name][delim][haplotype_id][delim][contig/scaffold_name] <br/>
e.g. For a sample/strain called CEA10 of haplotype 1 and a contig called 'CP097570.1' using the "#" seperator (recommended as it is rarely used anywhere else in sample names etc): <br/>
&nbsp; &nbsp; &nbsp; &nbsp;	>CEA10#1#CP097570.1  <br/>
_Note If you do not have haplotypes just put '1' is this space for all samples_ <br/>

This PanSN naming modification needs to be done for all assemblies in your dataset <br/>
For assemblies directly downloaded from NCBI, the trailing information (e.g. 'Aspergillus fumigatus CEA10 chromosome 8') after the contig accession can be left as is and will be ignored (due to the space seperation from the contig name) <br/>

Following this you can concatenate all assemblies into a single fasta file then compress with bgzip <br/>
e.g. ```cat assemblies_panSN/*.fa | bgzip > assemblies_panSN.concatenated.fa.gz``` <br/>
And voila, the primary input required for ```stargraph``` is ready. <br/>
Feed this _assemblies_panSN.concatenated.fa.gz_ file to ```stargraph``` (-a).


# Running stargraph first (wrapper included)

```stargraph``` requires some ```starfish``` input in order to run in its entirety <br/>
This includes:
	1. The _de-novo_ annotations of Tyrosine recombinases (used to elevate Presence/Absence variants (PAVs) to _Starship_-like regions (SLRs))
 	2. A final list of curated _Starship_ elements (combined with SLRs to generate the non-redundant dataset)

Therefore ```starfish``` needs to be run first <br/>
You can follow the starfish tutorial provided on the [github/wiki](https://github.com/egluckthaler/starfish) <br/>
Or you can use the wrapper provided by stargraph to run the primary steps required with default parameters <br/>
In this case the input used will be a list of paths to the PanSN renamed assemblies

	starfish_wrapper.sh -a assembly_files.txt


# Now running stargraph

 
	stargraph.sh -a assemblies_panSN.concatenated.fa.gz
	
	Required inputs:
	-a | --assemblies	XXXX

	Recommended inputs:
	-X | --XXXX		XXXX
	-t | --threads		Number of threads for tools that accept this option (default: 1)
	
	Optional parameters:
	-X | --XXXX  XXXX
	-p | --prefix		Prefix for output (default: name of assembly file (-a) before the fasta suffix)
	-o | --output		Name of output folder for all results (default: fusemblr_output)
	-c | --cleanup		Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
	-h | --help		Print this help message

XXXXX
