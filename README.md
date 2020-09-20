# MAGUS
Multiple Sequence Alignment  using Graph Clustering

- - - -

## Dependencies
MAGUS requires
* Python 3
* MAFFT (linux version is included)
* MCL (linux version is included)

If you would like to use some other version of MAFFT and/or MCL (for instance, if you're using Mac),
you will need to edit the MAFFT/MCL paths in configuration.py  
(I'll pull these out into a separate config file to make it simpler).

- - - -

## Getting Started
Please navigate to the "example" directory to get started with some sample data.  
A few basic ways of running MAGUS are shown below.  

**-o** specifies the output merged alignment path  
**-s** specifies the directory with subalignment files  
**-d** specifies the working directory for GCM's intermediate files, like the graph, clusters, log, etc.  

**python3 ../magus.py -d outputs -s subalignments -o magus_result.txt**

**-d** can be omitted; the working directory will be the directory of the ouput file

**python3 ../magus.py -s subalignments -o magus_result.txt**

**-b** specifies the directory with user-provided backbone alignment files. If omitted, they will be generated with MAFFT

**python3 ../magus.py -d outputs -s subalignments -b backbones -o magus_result.txt**

Instead of passing a directory with **-b**, you can pass a list of backbone files

**python3 ../magus.py -d outputs -s subalignments -b backbones/backbone_1.txt backbones/backbone_2.txt -o magus_result.txt**

**-r** and **-m** specify the number of MAFFT backbones and their maximum size, respectively. Default to 10 and 200.  
**-f** specifies the MCL inflation factor. Defaults to 4.0

**python3 ../magus.py -d outputs -s subalignments -r 10 -m 200 -f 2.5 -o magus_result.txt**

Instead of passing a directory with **-s**, you can pass a list of subalignment files

**python3 ../magus.py -d outputs -s subalignments/s_1.txt subalignments/s_2.txt subalignments/s_3.txt subalignments/s_4.txt subalignments/s_5.txt subalignments/s_6.txt subalignments/s_7.txt subalignments/s_8.txt subalignments/s_9.txt subalignments/s_10.txt subalignments/s_11.txt subalignments/s_12.txt subalignments/s_13.txt subalignments/s_14.txt -o magus_result.txt**

- - - -

## Things to Keep in Mind

* GCM will not overwrite existing backbone, graph and cluster files.  
Please delete them/specify a different working directory to perform a clean run.
* Related issue: if GCM is stopped while running MAFFT, MAFFT's output backbone files will be empty.  
This will cause errors if GCM reruns and finds these empty files.
* A large number of subalignments (>100) will start to significantly slow down the ordering phase, especially for very heterogenous data.  
I would generally disadvise using more than 50 subalignments, unless the data is expected to be well-behaved.  
