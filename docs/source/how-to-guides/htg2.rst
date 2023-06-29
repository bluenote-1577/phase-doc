Guide: understanding output phasings and their pitfalls
=======================================

This guide gives a deeper explanation for how to interpret the phasings output by floria and common pitfalls. 

It would be beneficial to do :doc:`../tutorials/tut1` before reading this guide. 

How many strains are there in my contig?
---------------------------------------

This is the question that the ``contig_ploidy_info.tsv`` file is trying to answer. A reference for this file can be found at :ref:`contig_ploidy_info`.

Based on our experience, the ``average_global_ploidy`` metric and the ``whole_contig_multiplicity`` metrics offers a good indication of the number of strains present. 

Essentially, ``average_global_ploidy`` tells you how many times a SNP is covered by vartigs on average. ``whole_contig_multiplicity`` tells you how many vartigs cover your contig. It is important you use both pieces of information. 

**For long reads**: The ``average_global_ploidy`` gives a pretty good indication of how many strains there are, and the ``whole_contig_multiplicity`` tends to be very similar to ``average_global_ploidy`` **IF** strains are actually present. 

However, it can be the case that ``average_global_ploidy`` >> ``whole_contig_multiplicity``, which indicates that only part of the genome is being phased, in which case it is likely that only part of the genome has variation -- perhaps from a very distant, conserved sequence or a mobile element. 

**For short reads**: The ``average_global_ploidy`` gives a pretty good indication of how many strains there are, but it may underestimate the number of strains slightly due to short-reads being able to resolve less strains. 

The ``whole_contig_multiplicity`` tends to be much less than ``average_global_ploidy`` -- this lies in short-read vartigs covering less of the entire contig. We recommend making sure that ``whole_contig_multiplicity`` is high enough, say > 0.1 at least. 

Spurious variants cause high "ploidy"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider the following "null" test. We generated synthetic reads for a single K. Pneumoniae genome, **so there is no strain heterogeneity**, aligned and called variants using freebayes (with no filters), and ran floria. 


.. code-block:: sh

    cat floria_out_dir/contig_ploidy_info.tsv
    -------------------------------------------
    contig	average_global_ploidy	whole_contig_multiplicity	approximate_coverage_ignoring_indels	average_local_ploidy	total_vartig_bases_covered	average_global_ploidy_min1hapq	average_local_ploidy_min1hapq	avg_err
    NZ_CP081897.1	1.754	0.015	15.251	1.971453287197232	80785	0.054	1.842	0.0000

 The global ploidy is still > 1 here, but notice that ``whole_contig_multiplicity`` is very small. 

The global ploidy is computed by checking how many times each SNP is covered by vartigs. However, there are almost no variants in the vcf file, because there is only 1 strain (only about 1000 "variants", which are really miscalls, for 5 Mb genome). floria phases the sections with "variation", and within these sections, the ploidy is about 2.

However, if we look at ``whole_contig_multiplicity``, we see only 0.015 of the contig is covered by vartigs. This indicates the vartigs are very sparse, so not much phasing is being done. 

**Conclusion:** Be careful on blindly relying on ``average_global_ploidy``. Looking at ``whole_contig_multiplicity`` gives a better understanding of what is going on. 

Species cross mappings cause high "ploidy" 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose we had both E. coli and K. pneumoniae in our community. Because E. coli and K. pneumoniae are not too distant (about 78% ANI), reads from E. coli can map to K. pneumoniae, which we call "cross mapping" here. 

To simulate cross mapping, we can simulate E. coli reads and K. pneumoniae reads, map them to K. pneumoniae (which share about 78% ANI). Now let's see how phasing works:

.. code-block:: sh

    cat cross_map_example/contig_ploidy_info.tsv
    -------------------------------------------
    contig	average_global_ploidy	whole_contig_multiplicity	approximate_coverage_ignoring_indels	average_local_ploidy	total_vartig_bases_covered	average_global_ploidy_min1hapq	average_local_ploidy_min1hapq	avg_err
    NZ_CP081897.1	1.774	0.171	18.393	1.9017563501504444	895083	1.514	1.779	0.0021

As we can see, the ploidy is still about 2, since there are two species present. However, the contig multiplicity is again quite small, indicating that many of the reads do not map and hence much of the contig does not get phased. However, still almost 1 MB of vartigs are present. 

**Conclusion:** Be careful of cross mappings. It's best to use a varied pool of reference genomes to avoid this artefact. 

Short reads give smaller contig multiplicity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We ran floria with short reads on a 3 strain simulated K. pneumoniae community. Here are the results:


.. code-block:: sh

    cat short_three_kpneumoniae_strains/contig_ploidy_info.tsv
    ----------------------------------------------------------
    contig	average_global_ploidy	whole_contig_multiplicity	approximate_coverage_ignoring_indels	average_local_ploidy	total_vartig_bases_covered	average_global_ploidy_min1hapq	average_local_ploidy_min1hapq	avg_err
    NZ_CP081897.1	2.808	1.677	105.219	2.948221469711441	8797680	2.732	2.899	0.0015

Notice that the contig multiplicity is quite a bit smaller than 3, the true number of strains. This is because short-reads can not span long enough gaps, so some gaps are left unfilled and have no vartigs covering them. 

On the other hand, for **long reads** on the same community, no such issues occurs.:

.. code-block:: sh

    cat long_three_kpneumoniae_strains/contig_ploidy_info.tsv
    ----------------------------------------------------------
    contig	average_global_ploidy	whole_contig_multiplicity	approximate_coverage_ignoring_indels	average_local_ploidy	total_vartig_bases_covered	average_global_ploidy_min1hapq	average_local_ploidy_min1hapq	avg_err
    NZ_CP081897.1	2.985	2.917	96.077	2.9852791178497586	15304581	2.982	2.983	0.0429


Now, on a **real community**, we got the following **short-read** results


.. code-block:: sh

    contig	average_local_ploidy	average_global_ploidy	approximate_coverage_ignoring_indels	total_vartig_bases_covered	whole_contig_multiplicity	average_local_ploidy_min1hapq	average_global_ploidy_min1hapq	avg_err
    NC_021016.1	2.106	1.936	298.657	2511481	0.806	2.034	1.787	0.0023
    NZ_AP024085.1	2.092	1.718	98.279	1193407	0.416	2.093	1.493	0.0038

We found that there were at least 2 strains present in each of these two contigs, but the multiplicity was even less than 1. This is because the strains are more similar to each other. When the strains are more similar, much of the contig can't be phased. 

**Conclusion:** the contig multplicity is helpful, but be aware of its nuances, especially when using short reads. For long-reads, the dispcrepancy between ploidy and multiplicity is less. 

Number of strains is more accurate when the coverage is higher
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Visualization is important. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^
