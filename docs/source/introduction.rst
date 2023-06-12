Introduction
=================

This page introduces a basic overview of how floria works and the problem it solves.

What is floria?
--------------

In microbial communities, strain-level genomic heteogeneity is often present. For example, not all E. coli in a human gut may have the exact same genome, but they will differ slightly.

If we have two distinct E. coli populations, then there are two **strains** presents in our community. *Strain* is an ambiguous term, but we will let a strain be a genome of a clonal population, all with "similar" genomes. When we sequence this community, reads may come from one of the two different E. coli strains. 

The problem  
^^^^^^^^^^^^

If we want to reconstruct the genomes in a metagenomic sample, we perform *metagenomic assembly* with software like MEGAHIT or metaflye. However, when heterogeneity is present, assemblies may become a mosaic of the different strains and do not represent an individual strain. 

One can also map the reads and call variants to understand variation within the sample. However, variant calling only gives a set of variants without any linkage information. That is, we can not guarantee which alleles reside physically on the same chromosome. This is the classical **haplotype phasing problem**. 

The solution: phasing
^^^^^^^^^^^^^^^^^^^^^^

floria solves this problem by **phasing** the reads in a metagenomic community. In our case, this means that it groups (i.e. **clusters**) reads such that every group corresponds to exactly **one strain**. This information can then be fed downstream for a variety of purposes. For example

#. We can assemble the groups that floria outputs, giving strain-level contigs instead of a mosaic contig
#. We can forgo assembly and look at the variation across the groups of reads. This tells us the sequence of SNPs that physically reside on one strain, for example. 

Overview of floria's algorithm
------------------------------

.. image:: img/method_diagram.png
  :width: 600
  :alt: floria algorithm diagram.

The above diagram outlines floria's main steps. floria requires at minimum

#.  A BAM file with reads mapped a set of contigs that we will try to phase and
#.  A VCF file with SNP information (other variants will be ignored).
#.  The fasta files for the contigs in the BAM file.

Given the above information, floria first uses SNP information in the VCF to annotate the reads in the BAM. floria then groups the reads into overlapping blocks along the contig, and phases each of the blocks. We group together these blocks using graph algorithms with network flows to form our final **haplosets**, which we define to be groups of reads coming from the indistinguishable haplotypes. 

The final output from floria are a collection of haplosets. We also provide SNP information along the haplosets (i.e. what are the likely SNPs for the strain represented by the haploset?), as well as visualization and QC information. 

What floria can and can not do
-----------------------------

What floria **can** do:

#. Efficient single-sample shotgun metagenomic phasing of short or long-reads on entire communities

    *   Short-reads should at least be paired-end, 2x100bp at minimum (longer is more ideal)  

    *   Long-reads can either be HiFi accurate reads or even noisier reads (still works even when identity is < 90%).

#.  Phase up to 5 strains by default and more strains if specified
#.  QC metrics and phasing statistics to manually visualize and confirm putative phasings.

What floria **can not** do:

#.  Phase consistently at very low coverage. Less than 5x becomes quite difficult.
#.  Too many strains at once. floria was designed to phase not more than single digit number of strains (5 by default).
#.  Call SNPs or map reads. To generate your VCF + BAM, you need to map reads and call variants by yourself. 
#.  floria does not give base-level assemblies. However, you can run an assembler on each of the output haplosets quite easily. 
