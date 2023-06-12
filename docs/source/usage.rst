Usage
=================

``floria -b sorted-and-indexed.bam -v vcf-with-header.vcf -r indexed-contigs.fa``

Mandatory inputs
-------

*   ``-b sorted-and-indexed.bam``: **Sorted** and **indexed** bam file (i.e. samtools sort, samtools index). Can be either paired-end short reads or long-reads. We assume one primary alignment per read, unique read names, and use MAPQs. Mainly tested with minimap2, but standard read aligners should be fine. 

*   ``-v vcf-with-header.vcf``: VCF with header info and SNPs. Only SNPs are considered in the vcf. Works with gzipped, bgzipped vcfs as well. 

*   ``-r indexed-contigs.fa``: **Indexed** fasta file (i.e. samtools faidx). Any contig in the bam file must appear in this fasta file. This is used for re-aligning reads to improve allele calls. 


For parameter descriptions, do ``floria -h``. Below are some useful parameters are controlling inputs/outputs for floria.

Outputs
---------

.. code-block:: sh
    results
    |   contig_ploidy_info.txt
    |   cmd.log
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
    │   |   │   0_part_paired1.fastq
    │   |   |   0_part_paired2.fastq
    |   |
    |   └───(debugging information)
    └───contig2_in_bam
        │   ...
        │   ...
```


