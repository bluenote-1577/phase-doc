Guide: understanding output phasings and their pitfalls
=======================================

This guide gives a deeper explanation for how to interpret the phasings output by floria and common pitfalls. 

It would be beneficial to do :doc:`../tutorials/tut1` before reading this guide. 

How many strains are there in my contig?
---------------------------------------

This is the question that the ``contig_ploidy_info.tsv`` file is trying to answer. For an introduction on how to do this, see :ref:`interpret_ploidy_usage` and the associated reference for this file. Paraphrasing 



Spurious variants cause high "ploidy"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Cross mappings
^^^^^^^^^^^^^^^


Filtering by ``whole_contig_multiplicity``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``total_vartig_bases_covered metric is also critical``. Occasionally, you may observe a contig with high ploidies but a small ``total_vartig_bases_covered value``. This could suggest mismappings or the appearance of false strains due to repetitive elements. If the ``total_vartig_bases_covered`` value is large, it's likely that multiple strains are present.

For instance, consider contig1, which is a genome larger than 2,000,000 bases. Its global ploidy is 1, suggesting the presence of only one strain. Furthermore, the number of bases covered is considerably lower than the genome size. This indicates that the variants and mappings could be spurious, or there might be only slight heterogeneity. On the other hand, contig2 is more likely to be a multi-strain contig.

Visualization is important. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Noisier reads require higher coverage for accurate information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Switch errors
^^^^^^^^^^^^


