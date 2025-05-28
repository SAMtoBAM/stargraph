# stargraph
A genome-graph based starship detection plugin to be combined with starfish

<p align="center" >
    <img src="https://github.com/SAMtoBAM/stargraph/blob/main/logo/stargraph.png" width=100%>
</p>

#[![Zenodo DOI](https://zenodo.org/badge/963417762.svg)](https://doi.org/10.5281/zenodo.15190276)
[![Anaconda_version](https://anaconda.org/samtobam/fusemblr/badges/version.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_platforms](https://anaconda.org/samtobam/fusemblr/badges/platforms.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda_downloads](https://anaconda.org/samtobam/fusemblr/badges/downloads.svg)](https://anaconda.org/samtobam/stargraph)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/fusemblr/badges/latest_release_date.svg)](https://anaconda.org/samtobam/stargraph)

**_stargraph_** is a tool that detects _Starship_-like regions (SLRs) using a genome-graph based and combines this with results from the more conservative tool ```starfish``` <br/>
The combination of both tools provides a more comprehensive view of genomic regions impacted by _Starships_

**_stargraph_** requires only the same input as ```starfish``` and some ```starfish``` output and therefore can be used easily in conjunction <br/>
As ```stargraph``` requies ```starfish``` input; The ```starfish``` pipeline should be run first in order to feed ```stargraph``` with both tyrosine recombinase (TyrR) and _Starship_ positions.

# Easy installation

	conda install samtobam::stargraph



# How to run

 
	stargraph.sh XXXXXXXX
	
	Required inputs:
	-X | --XXXX	XXXX

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
