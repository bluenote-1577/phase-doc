floria - metagenomic read-based strain phasing
===================================

floria is a software package for phasing metagenomic shotgun sequencing samples. 

Given a bam and a vcf file, floria clusters reads into strain-level clusters. floria offers the following:

#. Minimal parameter tuning and **automatic determination of strain number**

#. Works with both **short** or **long** reads

#. Can phase entire metagenomic communities,

   *  up to tens of gigabases of reads in less than a day with reasonable RAM (< 20 GB)
   *  phase a single organism at high depth in minutes 

#. A multitude of downstream information for interpreting clusters, including

   *  haplotagging reads for viewing in `IGV <https://software.broadinstitute.org/software/igv/>`_.
   *  summary statistics, such as putative number of strains for each contig, HAPQ -- an analogue of MAPQ for haplotyping quality, etc

To get started, check out the table of contents below. 

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   quick_start
   introduction
   usage
   tutorials
   how-to-guides
