Usage: input and output information
=================

.. code-block:: sh

    # required options
    floria -b sorted-and-indexed.bam -v vcf-with-header.vcf -r indexed-contigs.fa

    # with 20 threads and only phasing contig_1 to 'output_folder'
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
    |   (additional debugging folders/files)
    │
    └───contig1_in_bam
    │   │   contig1.haplosets
    |   |   contig1.vartigs
    │   │
    │   └───vartig_info
    │   |   │   0_hap.txt
    │   |   │   ...
    |   |
    │   └───long_reads (only when --output-reads is enabled)
    │   |   │   0_part.fastq
    │   |   |   ...
    |   |
    │   └───short_reads (only when --output-reads is enabled)
    │       │   0_part_paired1.fastq
    │       |   0_part_paired2.fastq
    |   
    └───contig2_in_bam
        │   ...
        │   ...


Each contig in the bam file is phased independently of the other contigs and has its own set of outputs. Below we describe the outputs for each contig. 

For a deeper discussion of the output information and potential pitfalls, see :doc:`how-to-guides/htg2`. 

Contig ploidy information
-----------------------


The ``out-dir/contig_ploidy_info.tsv`` file is extremely useful for characterizing the strain heterogeneity of your community at a glance. 

.. code-block:: sh

    contig  average_local_ploidy    average_global_ploidy   approximate_coverage_ignoring_indels    total_vartig_bases_covered    average_local_ploidy_min1hapq   average_global_ploidy_min1hapq
    contig1   1.706   0.971   17.739  194971  1.680   0.741
    contig2   2.509   2.351   69.065  3438158 2.437   2.231
    ...
#. ``contig``: Contig name
#. ``average_local_ploidy``: The local ploidy is the estimated ploidy of the blocks (see algorithm details in :doc:`introduction`) given that that the block passes floria's filtering thresholds. This is always greater than 1. 
#. ``average_global_ploidy``: The global ploidy is the estimated ploidy is the average SNP multiplicity over the contig. The SNP multiplicity is how many times a SNP is covered by haplosets. This can be < 1 because blocks which have 0 ploidy, i.e. do not have any SNPs or reads passing filters, are included in this metric. 
#. ``approximate_coverage_ignoring_indels``: The average coverage of the SNPs given that the SNP is covered by at least one read. Reads with many indels bias this metric down slightly because they may not cover SNPs properly. 
#. ``total_vartig_bases_covered``: How many bases are covered by vartigs. For example, if a contig has 4 strains, this will be about 4 times the contig length. However, it will be lower then certain parts of the contig are not covered by some of the strains. 
#. ``..._min1hapq``: The same statistics, but ignoring vartigs with 0 HAPQ.

Interpreting the ploidy information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In our experience, the most useful metric is the ``average_global_ploidy`` metric, which gives a good sense of how many strains there are. Short-reads tend to capture less strains, so a rough rule of thumb is that ``average_global_ploidy`` = 2.5 means there are probably 3 strains present for short reads.

The ``total_vartig_bases_covered`` metric is also important. Occasionally, you will see a contig with high ploidies, but ``total_vartig_bases_covered`` small. This may indicate mismappings or false strains appearing due to repetitive elements. If ``total_vartig_bases_covered`` is large, you probably have multiple strains present. 

For example, contig1 is a genome of size > 2,000,000 bases. Its global ploidy is 1, which seems to indicate that there's only 1 strain present. Furthermore, the number of bases covered is much less than the genome. Therefore, it is likely that the variants and mappings are spurious, or there is only a little bit of heterogeneity. contig2 is much more likely to be a multi-strain contig. 

Haplosets
--------

For the contig with name ``contig1``, the ``contig1.haplosets`` file is output and describes the strain-level haplosets (clusters of reads) for this contig. The file looks as follows:

.. code-block:: sh

    >HAP0_out-dir/contig1   SNPRANGE:1-6    BASERANGE:772-5000    COV:49.371  ERR:0.075   HAPQ:47   REL_ERR:1.35
    read_name1  first_snp_covered   last_snp_covered
    read_name2  first_snp_covered   last_snp_covered
    ...
    >HAP1_out-dir/contig1   SNPRANGE:7-11    BASERANGE:5055-6500    COV:25.012  ERR:0.050   HAPQ:15   REL_ERR:1.11
    ...

The lines with ``>`` give statistics about the haploset, and the lines below are reads within the haploset. 

#.  ``>HAP0_out-dir/contig1``: The ``>`` symbol delimits groups of reads. HAP(XXX) indicates the haploset identifier, where XXX is an integer. out-dir is the name of the output directory and contig1 is the contig name.

#.  ``SNPRANGE``: Which SNPs are considered for this haploset during the algorithm. Inclusive range. For example, 1-6 indicates this haploset covered the 1st SNP and the 6th SNP (starting from index 1). 

#. ``BASERANGE``: Same as SNPRANGE, but using the base level locations (1-indexed) instead of SNP number. 

#. ``COV``: An estimate of the coverage for this haploset. This estimate comes from counting how many times a SNP was covered by the reads on this haploset. This may not be as accurate as a coverage one would get by considering base-level alignments.

#. ``ERR``: The SNP error rate within this haploset. For example, if 9/10 reads had the reference allele but the other read had the alternate, the ERR would be 1/10.

#. ``HAPQ``: A number from 0-60 indicating the confidence (higher is better) that this haploset is **not a duplicated, spurious haploset**. Analogous to MAPQ from read mapping. HAPQ is **not** an estimate of phasing goodness, just like how MAPQ is different than a Smith-Waterman score.

#. ``REL_ERR``: The relative error of this haploset compared to all haplosets within this contig. 1.35 means 35% higher than the average error (ERR), for example. 

#. ``read_name1  first_snp_covered   last_snp_covered``: The name of the read and the first/last SNP covered by the read (inclusive). 

Vartigs
------

The ``contig1.vartigs`` file gives the **vartigs**, which are analogous to base-level contigs but only display the SNPs instead of the bases. 

.. code-block:: sh

    >HAP0_out-dir/contig1  SNPRANGE:1-6    BASERANGE:772-5000    COV:49.371  ERR:0.075   HAPQ:47   REL_ERR:1.35
    ?11111
    >HAP1_out-dir/contig1   SNPRANGE:7-11    BASERANGE:5055-6500    COV:25.012  ERR:0.050   HAPQ:15   REL_ERR:1.11
    01111

The line with ``>`` is the same as for the haplosets. 

The line below indicates the consensus alleles on this haploset. ``0`` always indicates the reference allele, and ``1`` indicates the first alternate allele, ``2`` the second alternate allele, etc. ``?`` indicates this allele is not covered by any read. So for the strain represented by ``HAP0_out-dir/contig1``, the alleles are alternate for all SNPs between [2,6] except for the first SNP.

Additional vartig info is available in the ``out-dir/contig/vartig_info/`` folder. For the vartig HAPX, X = 0,1,2,.. floria outputs a vartig information file `X_hap.txt` in the following format:

.. code-block:: 

    >HAP0_out-dir/contig1       SNPRANGE:1-6
    1:770   ?       NA      
    2:1022  1       1:1     
    3:2007  1       0:1|1:2 
    4:2034  1       1:3  

The lines after the header are of the form ``snp_number:base    consensus_allele    NA_or_allele_and_support``. The first two columns are straightforward. The third column indicates how strongly each allele is supported. For example, SNP 2 has only 1 read supporting the 1 allele. SNP 3 has 1 read supporting the 0 allele ``(0:1)`` and it has 2 reads supporting the 1 allele ``(1:2)``, hence why the conensus is 1 for SNP 3. 

Read output
----------

If the option ``--output-reads`` is specified, then reads will be output in the ``out-dir/long_reads`` or ``out-dir/short_reads`` folders. We put paired-end reads in the short_reads folder, and single-end reads in the long-read folders (even if you use short single-end reads). 

For example, the read ``0_part.fastq`` are all of the reads in fastq format for the 0th haplotype, or HAP0 as labelled in the other files. 

Importantly, the **long-reads are trimmed against the haplosets**. This means the following: consider a read of length 1500 that is mapped from [200,1700] on the contig. If this read is in HAP0, and HAP0's BASERANGE is 500-1000, we only output the read's bases that are mapped between [500,1000] on the genome. This is so that we don't get poor assemblies near the edges of the haplosets. 
