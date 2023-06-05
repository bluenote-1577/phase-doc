Important concepts
=================

This page introduces a basic overview of how floria works and the important terms we will work with. 

What is floria?
--------------

In microbial communities, strain-level genomic heteogeneity is often present. For example, not all E. coli in a human gut may have the exact same genome, but they will differ slightly. If we have two distinct E. coli populations, then there are two **strains** presents in our community. When we sequence this community, reads may come from one of the two different E. coli strains. 


The problem  
^^^^^^^^^^^^

If we want to reconstruct the genomes in a metagenomic sample, we perform *metagenomic assembly* with software like MEGAHIT or metaflye. However, when heterogeneity is present, assemblies may become a mosaic of the different strains and do not represent an individual strain. 

Recently, long, accurate PacBio HiFi reads give a promising solution to this problem, allowing for strain-level de novo assemblies. However, HiFi PacBio reads are more expensive than Nanopore or Illumina reads and may give lower throughputs. 

One can also map the reads and call variants to understand variation within the sample. However, variant calling only gives a set of variants without any linkage information. That is, we can not guarantee which alleles reside physically on the same chromosome. This is the classical **haplotype phasing problem**. 

The solution: phasing
^^^^^^^^^^^^^^^^^^^^^^

floria solves this problem by **phasing** the reads in a metagenomic community. In our case, this means that it groups (i.e. **clusters**) reads such that every group corresponds from exactly **one strain**. This information can then be fed downstream for a variety of purposes. For example

#. We can assemble the groups that floria outputs, giving strain-level contigs instead of a mosaic contig
#. We can forgo assembly and look at the variation within the group of reads. This tells us the sequence of SNPs that physically reside on one strain, for example. 

What floria can and can not do
-----------------------------

What is a "strain"?
------------------

Overview of floria's algorithm
------------------------------

Inputs to floria
------------

Outputs of floria
----------------
