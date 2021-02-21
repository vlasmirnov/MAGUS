# MAGUS
Multiple Sequence Alignment using Graph Clustering

- - - -

## Purpose and Functionality
MAGUS is a tool for piecewise large-scale multiple sequence alignment.  
The dataset is divided into subsets, which are independently aligned with a base method (currently MAFFT -linsi). These subalignments are merged together with the Graph Clustering Merger (GCM). GCM builds the final alignment by clustering an alignment graph, which is constructed from a set of backbone alignments. This process allows MAGUS to effectively boost MAFFT -linsi to over a million sequences.

The basic procedure is outlined below. Steps 4-7 are GCM.
1. The input is a set of unaligned sequences. Alternatively, the user can provide a set of multiple sequence alignments and skip the next two steps.
2. The dataset is decomposed into subsets.
3. The subsets are aligned with MAFFT -linsi. 
4. A set of backbone alignments are generated with MAFFT -linsi (or provided by the user).
5. The backbones are compiled into an alignment graph.
6. The graph is clustered with MCL.
7. The clusters are resolved into a final alignment.

- - - -

## Dependencies
MAGUS requires
* Python 3
* MAFFT (linux version is included)
* MCL (linux version is included)
* FastTree and Clustal Omega are needed if using these guide trees (linux versions included) 

If you would like to use some other version of MAFFT and/or MCL (for instance, if you're using Mac),
you will need to edit the MAFFT/MCL paths in configuration.py  
(I'll pull these out into a separate config file to make it simpler).

- - - -

## Getting Started
Please navigate your terminal to the "example" directory to get started with some sample data.  
A few basic ways of running MAGUS are shown below.  
Run "magus.py -h" to view the full list of arguments. 

**Align a set of unaligned sequences from scratch**  
*python3 ../magus.py -d outputs -i unaligned_sequences.txt -o magus_result.txt*  

*-o* specifies the output alignment path  
*-d* (optional) specifies the working directory for GCM's intermediate files, like the graph, clusters, log, etc.  

**Merge a prepared set of alignments**  
*python3 ../magus.py -d outputs -s subalignments -o magus_result.txt*  

*-s* specifies the directory with subalignment files. Alternatively, you can pass a list of file paths.   

- - - -

## Controlling the pipeline

**Specify subset decomposition behavior**  
*python3 ../magus.py -d outputs -i unaligned_sequences.txt -t fasttree --maxnumsubsets 100 --maxsubsetsize 50 -o magus_result.txt*  

*-t* specifies the guide tree method to use, and is the main way to set the decomposition strategy.  
Available options are fasttree (default), parttree, clustal (recommended for very large datasets), and random.  
*--maxnumsubsets* sets the desired number of subsets to decompose into (default 25).  
*--maxsubsetsize* sets the threshold to stop decomposing subsets below this number (default 50).  
Decomposition proceeds until maxnumsubsets is reached OR all subsets are below maxsubsetsize.

**Specify beckbones for alignment graph**  
*python3 ../magus.py -d outputs -i unaligned_sequences.txt -r 10 -m 200 -o magus_result.txt*  
*python3 ../magus.py -d outputs -s subalignments -b backbones -o magus_result.txt*  

*-r* and *-m* specify the number of MAFFT backbones and their maximum size, respectively. Default to 10 and 200.  
Alternatively, the user can provide his own backbones; *-b* can be used to provide a directory or a list of files.

**Specify graph trace method**  
*python3 ../magus.py -d outputs -i unaligned_sequences.txt --graphtracemethod mwtgreedy -o magus_result.txt*  

*--graphtracemethod* is the flag that governs the graph trace method. Options are minclusters (default and recommended), fm, mwtgreedy (recommended for very large graphs), rg, or mwtsearch.

**Unconstrained alignment**  
*python3 ../magus.py -d outputs -i unaligned_sequences.txt -c false -o magus_result.txt*  

By default, MAGUS constrains the merged alignment to induce all subalignments. This constraint can be disabled with *-c false*.  
This drastically slows MAGUS and is strongly not recommended above 200 sequences. 

- - - -

## Things to Keep in Mind

* MAGUS will not overwrite existing backbone, graph and cluster files.  
Please delete them/specify a different working directory to perform a clean run.
* Related issue: if MAGUS is stopped while running MAFFT, MAFFT's output backbone files will be empty.  
This will cause errors if MAGUS reruns and finds these empty files.
* A large number of subalignments (>100) will start to significantly slow down the ordering phase, especially for very heterogenous data.  
I would generally disadvise using more than 100 subalignments, unless the data is expected to be well-behaved.  

- - - -

## Related Publications

* Original MAGUS paper: ___Smirnov, V. and Warnow, T., 2020. MAGUS: Multiple Sequence Alignment using Graph Clustering. Bioinformatics. https://doi.org/10.1093/bioinformatics/btaa992___
* GCM-MWT paper:
* MAGUS on ultra-large datasets: 
