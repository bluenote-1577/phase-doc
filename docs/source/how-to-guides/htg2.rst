Guide: understanding output phasings and their pitfalls
=======================================

This guide gives a deeper explanation for how to interpret the phasings output by floria and common pitfalls. 

It would be beneficial to do :doc:`../tutorials/tut1` before reading this guide. 

The first thing you should do after running floria is to examine the ``contig_ploidy_info.tsv`` file found in the output folder. It gives very useful high-level statistics for your results. 


Filtering by ``whole_contig_multiplicity``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``total_vartig_bases_covered metric is also critical``. Occasionally, you may observe a contig with high ploidies but a small ``total_vartig_bases_covered value``. This could suggest mismappings or the appearance of false strains due to repetitive elements. If the ``total_vartig_bases_covered`` value is large, it's likely that multiple strains are present.

For instance, consider contig1, which is a genome larger than 2,000,000 bases. Its global ploidy is 1, suggesting the presence of only one strain. Furthermore, the number of bases covered is considerably lower than the genome size. This indicates that the variants and mappings could be spurious, or there might be only slight heterogeneity. On the other hand, contig2 is more likely to be a multi-strain contig.

Visualization is important. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Noisier reads require higher coverage for accurate information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spurious variants cause high "ploidy"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Cross mappings
^^^^^^^^^^^^^^^


Switch errors
^^^^^^^^^^^^


