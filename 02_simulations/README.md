## MELT runs with simulated TE insertion data

Workflow for estimating MELT accuracy and sensitivity for your desired transposable element in a given genome.
This is a version of the ESAT pipeline used by Lammers et al. 2019 (https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0143-2), modified for more automation.

Use batch_loop_ESAT.py and batch_ESAT_MELT_runs.py to generate job submission scripts for your files.

Required software:
* MELT
* SAMtools
* Picard toolkit
* SimSeq
* BWA
* BEDtools

Required files are:
1. TE consensus FASTA
2. A representative scaffold or chunk of the desired genome.
3. Index file of the representative scaffold or chunk.
4. Other required input files for a MELT run (i.e. a genetrack.bed file (can be an empty file); MELT ZIP files)

Note: Directories are hard-coded variables in some scripts.

