## MELT runs with simulated TE insertion data

Workflow for estimating MELT accuracy and sensitivity for your desired TE or TE type (i.e. DNA transposons) in a given genome.
This is a version of the ESAT pipeline used by [Lammers et al. 2019](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0143-2), modified for more automation.

Use batch_loop_ESAT.py to generate job submission scripts for your files. View script for required arguments.
Note: Directories and job queues are hard-coded variables.

Required software:
* SimSeq (https://github.com/jstjohn/SimSeq)
* MELT
* SAMtools
* Picard toolkit
* BWA
* BEDtools

Required files are:
* TE consensus FASTA
* A representative scaffold or chunk of the desired genome.
* Index file of the representative scaffold or chunk.
* Other required input files for a MELT run (i.e. a genetrack.bed file (can be an empty file); MELT ZIP files)
