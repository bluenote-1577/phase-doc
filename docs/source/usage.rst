Usage
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

Outputs
---------

You can specify the output folder name with the ``-o results_folder_name`` option. The structure of the output directory is shown below.

.. code-block:: sh

    results
    |   contig_ploidy_info.txt
    |   cmd.log
    |   (additional debugging folders/files)
    │
    └───contig1_in_bam
    │   │   contig1_haplosets.txt
    |   |   contig1_vartigs.txt
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

Haplosets
^^^^^^^^

For the contig with name ``contig1``, the ``contig1_haplosets.txt`` file is output and describes the strain-level haplosets (clusters of reads) for this contig. The file looks as follows:

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

#. ``REL_ERR``: The relative error of this haploset compared to all haplosets within this contig. 1.35 means 35% higher than the average error, for example. 

#. ``read_name1  first_snp_covered   last_snp_covered``: The name of the read and the first/last SNP covered by the read (inclusive). 
