Guide: how to generate a good .vcf and .bam file for phasing
=======================================

This guide goes over how to generate good variant calls and mappings for floria. 

We have also implemented a pipeline called Floria-PL at https://github.com/jsgounot/Floria_analysis_workflow (a.k.a the "Production" pipeline). This pipeline goes from reads-to-phasings/assemblies and generates a VCF and BAM file as an intermediate step, including an optional kraken mapping approach. This may be useful if you want an end-to-end pipeline, but requires more configuration. 

Alignment: Use-case 1 - metagenomic contig phasing
-------------------------------------------------

If you're phasing a metagenomic assembly, simply map your reads to the resulting assembly using your favourite aligner to generate a BAM file. We used minimap2 in all of our experiments. 

For example:

.. code-block:: sh

  minimap2 -a assembly.fa my_nanopore_reads.fq -t 20 | samtools sort -o nanopore_on_assembly.bam
  samtools index nanopore_on_assembly.bam
  samtools faidx assembly.fa 

.. note::

  You must have an **indexed and sorted** bam file and an **indexed** fasta file

Metagenome-assembled genomes?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Floria works on contigs and does not use binning information. However, you can choose to only map reads to metagenome-assembled genomes (MAGs) if you want. In this case, concatenate all MAGs into one fasta file and map all reads to that fasta file. 


Alignment: Use-case 2 - reference-based phasing
--------------------------------------

In the following cases, you may consider mapping reads to reference genomes instead of use-case 1 (assembling + mapping). 

* You already have a set of reference genomes that represent your metagenome well (e.g. by using a taxonomic profiler)
* Your reads come from a single species but contains multiple strains 

In both cases simply concatenate all genomes (if there are multiple) into one reference file and align:

.. code-block:: sh

  cat my_genomes/* > concatenated_genomes.fa
  minimap2 -a concatenated_genomes.fa nanopore_reads.fq -t 20 | samtools sort -o result.bam
  samtools index result.bam
  samtools faidx concatenated_genomes.fa
  ...

.. note::

  You must have an **indexed and sorted** bam file and an **indexed** fasta file

Making sure your reference genomes are comprehensive (metagenomics)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is important that your fasta file represents your metagenome well. 

If you're missing reference genomes, **cross-mappings** may occur. Cross-mapping is when reads from an unrelated species maps to the wrong species' reference genome (e.g. Klebsiella africana reads may map to a Klebsiella pneumoniae genome). 

To handle cross-mappings, consider running a taxonomic profiler first to find what species are present, and then using those reference genomes. One option is to  use sylph (https://github.com/bluenote-1577/sylph) as a way of finding reference genomes are in your community. 

Variant calling: using a good low-frequency SNP caller
---------------------------------------

After obtaining a BAM file, you must call SNPs (single nucleotide polymorphisms) before using floria. 

Unlike traditional diploid variant calling, SNPs may occur in a community at low frequency (instead of being 50 or 100 percent frequency). It is important to use a ploidy-agnostic variant caller. 

Short reads
^^^^^^^^^^^

For short reads consider:

* freebayes using the --pooled-continuous option. 
* LoFreq also works well.  

Long reads
^^^^^^^^^^

We have used longshot (https://github.com/pjedge/longshot) as a SNP caller with success, even though longshot assumes diploidy. Other potential variant callers that look promising, but we have not tested in depth:

* Clair-TO (https://github.com/HKU-BAL/ClairS-TO): a neural-network variable-frequency variant caller for long reads.
* Clair3 (https://github.com/HKU-BAL/Clair3): a **diploid** neural-network variant caller. 
* NanoCaller (https://github.com/WGLab/NanoCaller): a **diploid** neural-network variant caller. 

