Usage: input and output information
=================

.. code-block:: sh

    # required options
    floria -b sorted-and-indexed.bam -v vcf-with-header.vcf -r indexed-contigs.fa

    # with 20 threads, phasing only contig_1 to 'output_folder'
    floria -b map.bam -v snp.vcf -r refs.fa -G contig_1 -t 20 -o output_folder

    # output reads in fastq format, filter contigs with < 1000 SNPs, do not use supplementary alignments, require MAPQ 60.
    floria -b map.bam -v snp.vcf -r refs.fa --output-reads --snp-count-filter 1000 -X -m 60

Mandatory inputs
-------

*   ``-b sorted-and-indexed.bam``: **Sorted** and **indexed** bam file (i.e. samtools sort, samtools index). Can be either paired-end short reads or long-reads. We assume one primary alignment per read, unique read names, and use MAPQs. Mainly tested with minimap2, but standard read aligners should be fine (i.e. one primary aligment per read, supplementary flags, reasonable 0-60 MAPQs). 

*   ``-v vcf-with-header.vcf``: VCF with header info and SNPs. Only SNPs are considered in the vcf. Works with gzipped, bgzipped vcfs as well. 

*   ``-r indexed-contigs.fa``: **Indexed** fasta file (i.e. samtools faidx). Any contig in the bam file must appear in this fasta file. This is used for re-aligning reads to improve allele calls. 

For optional parameter descriptions, do ``floria -h``. For tuning floria's parameters to your problem of interest, see :doc:`how-to-guides/htg1`. 

.. _usage-outputs:

Outputs
---------

You can specify the output folder name with the ``-o out-dir`` option. The structure of the output directory is shown below.

.. code-block:: sh

    out-dir
    |   contig_ploidy_info.tsv
    |   cmd.log
    │
    └───contig1_in_bam
    │   │   contig1.haplosets
    |   |   contig1.vartigs
    |   |   reads_without_snps.tsv
    │   │   vartig_info.txt
    |   |
    │   └───long_reads (only when --output-reads is enabled)
    │   |   │   0_part.fastq (reads output only when the partition is non-empty)
    │   |   |   ...
    │   |   │   snpless.fastq
    |   |
    │   └───short_reads (only when --output-reads is enabled)
    │       │   0_part_paired1.fastq
    │       |   0_part_paired2.fastq
    |       |   ...
    |       |   snpless_paired1.fastq
    |       |   snpless_paired2.fastq
    |   
    └───contig2_in_bam
        │   ...
        │   ...


Each contig in the bam file is phased independently of the other contigs and has its own set of outputs. Below we describe the outputs for each contig. 

For a deeper discussion of the output information and potential pitfalls, see :doc:`how-to-guides/htg2`. 

Haplosets
--------

For the contig with name ``contig1``, the ``contig1.haplosets`` file is output and describes the strain-level haplosets (clusters of reads) for this contig. The file looks as follows:

For a contig named contig1, the output file ``contig1.haplosets`` provides detailed information about the strain-level haplosets (clusters of reads) for this contig. Here is an example of what the file might look like:

.. code-block:: sh

    >HAP0.out-dir/contig1   CONTIG:contig1  SNPRANGE:1-15   BASERANGE:771-3416  COV:2.429   ERR:0.088   HAPQ:30  REL_ERR:1.591
    read_name1  first_snp_covered   last_snp_covered
    read_name2  first_snp_covered   last_snp_covered
    ...
    >HAP1.out-dir/contig1   CONTIG:contig1  SNPRANGE:16-21  BASERANGE:5055-6500 COV:25.012  ERR:0.050   HAPQ:15 REL_ERR:1.11
    ...

In this file, lines beginning with ``>`` provide statistics about the haploset, and the subsequent lines list the reads contained within that haploset.

Here's a breakdown of the data provided:

#. ``>HAP0.out-dir/contig1``: The > symbol marks the start of a new group of reads. HAP(XXX) is the haploset identifier, where XXX is an integer. 'out-dir' refers to the name of the output directory, and 'contig1' is the contig name.
#. ``SNPRANGE``: This refers to the range of SNPs considered for this haploset during the algorithm. This range is inclusive. For example, 1-6 implies this haploset covers the 1st SNP to the 6th SNP (starting from index 1).
#. ``BASERANGE::: Similar to SNPRANGE but uses base-level locations (1-indexed) rather than SNP numbers.
#. ``COV``: This is an estimation of the coverage for this haploset, calculated by counting the number of times a SNP is covered by the reads in this haploset. This estimate might not be as accurate as coverage determined by considering base-level alignments.
#. ``ERR``: This represents the SNP error rate within this haploset. For instance, if 9 out of 10 reads carry the reference allele and one read has the alternate, the ERR would be 1/10.
#. ``HAPQ``: This is a confidence score ranging from 0-60 that indicates how likely it is that this haploset is not a duplicate or spurious haploset. This is similar to MAPQ from read mapping. Note that HAPQ is not an estimate of phasing quality, just as MAPQ differs from a Smith-Waterman score.
#. ``REL_ERR``: This denotes the relative error of this haploset compared to all haplosets within this contig. For example, 1.35 means the error is 35% higher than the average error (ERR).
#. ``read_name1 first_snp_covered last_snp_covered``: This provides the name of the read and the range of SNPs covered by that read (inclusive).

reads_without_snps.tsv
^^^^^^^^^^^^^^^^^^^^^

Reads in blocks that have little-to-no variation are written to this file. This happens due to the following reasons:

#. There is little variation present in your contig.
#. A strain and your reference genome have long, similar regions.
#. Your reference contig has a large inserted, unique sequence relative to another strain, so no reads from the other strain get mapped here. 
#. Variant calling was not successful for this region, so it looks like your reads don't cover any SNPs


Vartigs
------

The file named ``contig1.vartigs`` provides the vartigs, which are analogous to base-level contigs but display only the SNPs rather than all bases.

.. code-block:: sh

    >HAP0.out-dir/contig1   SNPRANGE:1-6    CONTIG:contig1   BASERANGE:772-5000    COV:49.371  ERR:0.075   HAPQ:47   REL_ERR:1.35
    ?11111
    >HAP1.out-dir/contig1   SNPRANGE:7-11  CONTIG:contig2   BASERANGE:5055-6500    COV:25.012  ERR:0.050   HAPQ:15   REL_ERR:1.11
    01111


The line starting with ``>`` follows the same format as for the haplosets. The subsequent line provides the consensus alleles for this haploset. Here, ``0`` always represents the reference allele, ``1`` stands for the first alternate allele, ``2`` for the second alternate allele, and so on. ``?`` implies that no read covers this allele. Thus, for the strain represented by HAP0.out-dir/contig1, the alleles are alternate for all SNPs in the range [2,6], except for the first SNP, which is not covered by a read.

vartig_info.txt
--------------

You can find additional vartig information in the ``out-dir/contig/vartig_info.txt`` file. For each vartig HAPX, where X = 0,1,2,..., floria generates additional vartig information.

.. code-block:: 

    >HAP0.out-dir/contig1       SNPRANGE:1-6
    1:770   ?       NA      
    2:1022  1       1:1     
    3:2007  1       0:1|1:2 
    4:2034  1       1:3  
    ...
    >HAP1.out-dir/contig1       SNPRANGE:4-7
    4:2034  0       0:25
    5:2037  0       0:25
    6:2217  0       0:28
    7:2220  0       0:29

The lines after the header are of the form ``snp_number:base    consensus_allele    NA_or_allele_and_support``. The first two columns are straightforward. The third column indicates how strongly each allele is supported. For example, SNP 2 has only 1 read supporting the 1 allele. SNP 3 has 1 read supporting the 0 allele ``(0:1)`` and it has 2 reads supporting the 1 allele ``(1:2)``, hence why the conensus is 1 for SNP 3. 

Contig ploidy information
-----------------------

The ``out-dir/contig_ploidy_info.tsv`` file is extremely useful for characterizing the strain heterogeneity of your community at a glance. 

.. code-block:: sh

    contig  average_local_ploidy    average_global_ploidy   approximate_coverage_ignoring_indels    total_vartig_bases_covered    average_local_ploidy_min1hapq   average_global_ploidy_min1hapq
    contig1   1.706   0.971   17.739  194971  1.680   0.741
    contig2   2.509   2.351   69.065  3438158 2.437   2.231
    ...


The following are the interpretations of each column:

#. ``contig``: The contig's name.
#. ``average_local_ploidy``: This refers to the estimated ploidy of the blocks (see algorithm details in :doc:introduction) that pass floria's filtering thresholds. This value is always greater than 1.
#. ``average_global_ploidy``: This represents the average SNP multiplicity across the contig, which is the estimated ploidy. SNP multiplicity is the frequency of a SNP's coverage by haplosets. This can be less than 1 because blocks with 0 ploidy, i.e., blocks lacking any SNPs or reads that pass filters, are included in this metric.
#. ``approximate_coverage_ignoring_indels``: This is the average coverage of the SNPs. Reads with many indels can slightly decrease this metric as they may not properly cover SNPs.
#. ``total_vartig_bases_covered``: This is the total number of bases covered by vartigs. For instance, if a contig has 4 strains, this number will be approximately four times the contig length. However, it may be lower if certain parts of the contig are not covered by some strains.
#. ..._min1hapq``: This represents the same statistics but disregards vartigs with 0 HAPQ.

Interpreting the Ploidy Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Based on our experience, the ``average_global_ploidy`` metric is particularly useful as it offers a good indication of the number of strains present. As short reads tend to capture fewer strains, a rough guideline is that an ``average_global_ploidy`` of 2.5 likely indicates the presence of three strains for short reads.

The ``total_vartig_bases_covered metric is also critical``. Occasionally, you may observe a contig with high ploidies but a small ``total_vartig_bases_covered value``. This could suggest mismappings or the appearance of false strains due to repetitive elements. If the ``total_vartig_bases_covered`` value is large, it's likely that multiple strains are present.

For instance, consider contig1, which is a genome larger than 2,000,000 bases. Its global ploidy is 1, suggesting the presence of only one strain. Furthermore, the number of bases covered is considerably lower than the genome size. This indicates that the variants and mappings could be spurious, or there might be only slight heterogeneity. On the other hand, contig2 is more likely to be a multi-strain contig.


.. _read-outputs:

Read output
----------

If you specify the ``--output-reads option``: 

#. Reads will be output in the ``out-dir/long_reads`` or ``out-dir/short_reads folders``. 
#. Paired-end reads are placed in the short_reads folder, and single-end reads are in the long-reads folders (even for short single-end reads). 
#. A special file ``out-dir/long_reads/snpless_reads.fastq`` represent reads that pass filtering thresholds but reside in blocks on the reference that have little to no variation (i.e. from ``reads_without_snps.tsv``. If you attempt to phase a genome with almost no variation (e.g. only false positive, scattered SNPs) then most of your reads will be here. 

For example, the file ``1_part.fastq`` contains all of the reads in fastq format for the 1th haplotype, also labeled as HAP1 in other files.

Importantly, **long-reads are trimmed against the haplosets**. For example, consider a read of length 1500 that is mapped from [200,1700] on the contig. If this read is part of HAP1, and HAP1's BASERANGE is 500-1000, we only output the portions of the read that are mapped between [500,1000] on the genome. This approach helps to prevent poor assemblies near the edges of the haplosets.


