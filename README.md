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


# Preprocessing of input data

To help with building genome-graphs each contig for each assembly needs sample/haplotype information stored in the header <br/>
To do this we use the PanSN specifications [sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]; explained in detail [here](https://github.com/pangenome/PanSN-spec) <br/>
This means the header for each contig will look something like this _>CEA10#1#CP097570.1_ for a sample called CEA10 of haplotype 1 and a contig called 'CP097570.1' using the "#" seperator (recommended as it is rarely used anywhere else in sample names etc) <br/>
If you do not have haplotypes just put '1' is this space for all samples <br/>
This needs to be done for all assemblies wanting to be analysed <br/>
For assemblies directly downloaded from NCBI, the trailing information (e.g. 'Aspergillus fumigatus CEA10 chromosome 8') after the contig accession can be left as is and will be ignored (due to the space seperation from the contig name) <br/>

Following this you can concatenate all the assemblies used into a single fasta file then compress with bgzip <br/>
e.g. ```cat assemblies_panSN/*.fa | bgzip > assemblies_panSN.concatenated.fa.gz``` <br/>
And voila, the only input required for _stargraph_ is ready. Feed this _assemblies_panSN.concatenated.fa.gz_ file to stargraph (-a).


# How to run

 
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
