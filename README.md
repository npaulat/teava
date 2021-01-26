# TEAVA - TE-Associated Variant Analysis
A pipeline for a local single nucleotide variant (SNV) and indel enrichment analysis near polymorphic transposable element (TE) sites.

This workflow is based on 1) TE annotation via RepeatMasker, 2) BWA-mapped Illumina paired-end reads, 3) polymorphic TE site identification via MELT, and 4) indel and SNV calls via FreeBayes. It uses a combination of Shell scripts, Python, and R.

## General workflow
1. TE annotation of reference genome with RepeatMasker and selection of desired TEs
2. BWA read mapping and FreeBayes variant calling
3. Polymorphic TE identification with MELT
4. Polymorphic TE quality filtering and reformatting
5. SNV/Indel quality filtering and enrichment analysis

## Necessary programs
* Python v2.7/v3.6.0 (https://www.python.org/)
* R v3.6.2 (https://www.r-project.org/)
* RepeatMasker v4.0.9 (http://repeatmasker.org/)
* Bedtools2 v2.26.0 (https://bedtools.readthedocs.io/en/latest/)
* BWA v0.7.17 (http://bio-bwa.sourceforge.net/)
* Picard toolkit v2.5.0 (http://broadinstitute.github.io/picard/)
* SAMtools v1.9 (http://www.htslib.org/)
* BCFtools v1.9 (http://www.htslib.org/)
* FreeBayes v1.2.0 (https://github.com/freebayes/freebayes) - uses Python v2.7 until v1.3.2
* GATK v3.8.0 (https://gatk.broadinstitute.org/hc/en-us)
* VCFtools v0.1.16 (https://vcftools.github.io/index.html)
* MELT v2.2.0 (https://melt.igs.umaryland.edu/) - requires Java and Bowtie2

NOTE: The program versions listed are those used in the original pipeline, but newer versions should be compatible.
