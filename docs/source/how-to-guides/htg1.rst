Guide: optimizing parameters for phasing
=======================================

.. _htg1-crosslink:

This guide explains a few key parameters used in floria that are not obvious from their help description. 

Error rate -e
^^^^^^^^^^^^^

The error rate is arguably floria's most important parameter. It is an estimate of the fraction of the alleles within your reads that are assigned to the wrong allele.

The -e parameter controls controls for how collapsed your strains may be, among other things. In theory, the lower the -e, the more strains you could potentially phase, at the cost of perhaps more spurious haplosets appearing. 

By default, floria v0.0.1 uses a heuristic method to estimate -e. The heuristic we settled on is simple: in the input bam file, we estimate a distribution of mismatch error rates for an arbitrary contig in the bam file by considering coverage information, and then take the 66th percentile mismatch error rate to be -e. We set -e to be a minimum of 0.01 by default. 

Our recommendation is to use default -e unless you're focusing on a specific stain of interest, where then playing around with -e and visualizing the results could be beneficial. 

Block length -l
^^^^^^^^^^^^^^

The block length controls for how contiguous your haplosets are versus their quality. Shorter blocks will give more haplosets, but they will be less likely to contain switch errors within the blocks. 

In floria v0.0.1, we set -l to be the 66th percentile read length by default, except when short-reads are detected, in which case we set -l to 500. 

For paired-end reads, if you have prior information about the fragment length, then perhaps changing -l to a larger or smaller value scaling with your fragment length can be beneficial. 

Beam search number of solutions -n
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This controls the number of solutions in floria's beam search. In principle, larger -n should increase phasing quality at the cost of speed. In practice, floria's beam search heuristics already do a pretty good job. 

Ploidy sensitivity -s
^^^^^^^^^^^^^^^^^^^^

This is an experimental option that controls how stringent the heuristic is for detecting the ploidy (i.e. the number of strains). You can set -s to 1, 2, or 3. Higher ploidy sensitivity should give a larger number of strains. 

We have not tested -s 1 or -s 3 extensively, but we recommend using this option if you really believe floria is not phasing the correct number of strains. 

--no-stop-heuristic
^^^^^^^^^^^^^^^^^^

floria determines the ploidy for the blocks by iteratively phasing different ploidies until the stop heuristic (also controlled by -s) is triggered or the error rate of the block is less than -e. This flag turns off the stopping heuristic, and only stops phasing when the error rate is less than -e. 

We do not recommend this option for noisier reads (i.e. Nanopore). This option may be useful for HiFi or accurate short-reads, where the stopping heuristic is not as well tested. 


--ignore-monomorphic
^^^^^^^^^^^^^^^^^^^

This option ignores SNPs in your VCF that have very low minor allele frequency, and mask the outputs with "?". The option is intended for when the population you're phasing does not actually have many of the variants specified by your VCF file, so you're not interested in analyzing SNPs that are not actually polymorphic.

This will change floria's actual results slightly as well. 



