Tutorial 1: running floria and visualizing its outputs
================================================

This minimal tutorial goes over floria's inputs and outputs by running floria on toy data. 

We will assume that you've installed floria. See the :doc:`quick_start` page if you have not installed floria. We will assume you have the following commands working properly. 

.. code-block:: sh

   # clone the repository if you have not done so
   git clone https://github.com/bluenote-1577/floria
   cd floria

   # assuming you've installed floria
   floria -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output


Inputs
^^^^^

The example files represent a mixture of synthetic Klebsiella Pneumoniae reads over three different strains. 

``test_long.bam`` represent synthetic nanopore reads from three different strains at 70x, 30x, and 10x coverage. We then mapped these reads to the ``MN-03.fa`` reference genome to generate the bam file. We only retained reads between positions 1-120,000 for this example. 

``test.vcf`` was obtained by calling SNPs on ``test_long.bam``. We used longshot, a long-read SNP caller. 

``MN-03.fa`` is the contig we mapped onto. 

Output log
^^^^^^^^^

Depending on your version of floria, you should see something like the following after running the command. 

.. code-block:: 

    2023-06-12T23:33:40.561Z INFO  [floria::file_reader] 33,50,66 non-hard clipped read length percentiles: 9943, 13066, 16494. If -l is not set, estimated -l is set to 16494.
    2023-06-12T23:33:40.561Z INFO  [floria::file_reader] If -e is not set, estimated -e is set to 0.04081632653061224.
    2023-06-12T23:33:40.563Z INFO  [floria] Preprocessing VCF/Reference
    2023-06-12T23:33:40.568Z INFO  [floria] Finished preprocessing in 5.20466ms
    2023-06-12T23:33:40.568Z INFO  [floria] Reading and realigning inputs for contig NZ_CP081897.1 (BAM/VCF).
    2023-06-12T23:33:40.664Z INFO  [floria] Contig NZ_CP081897.1 has 954 SNPs
    2023-06-12T23:33:40.664Z INFO  [floria] Reading inputs, realigning time taken 96.312178ms
    2023-06-12T23:33:40.664Z INFO  [floria] Local phasing with 20 threads...
    2023-06-12T23:33:41.041Z INFO  [floria] Phasing time taken 376.723842ms
    2023-06-12T23:33:41.041Z INFO  [floria] Solving flow problem...
    2023-06-12T23:33:41.041Z INFO  [floria] Flow solved in time 339.672µs
    2023-06-12T23:33:41.084Z INFO  [floria::file_writer] Final SNP error rate for all haplogroups is 0.05546256145299258
    2023-06-12T23:33:41.085Z INFO  [floria] Total time taken is 522.39064ms

The first line indicates what the ``-l`` parameter, the block length, is. We take the 66th percentile read length to be the block lengths in floria by default unless it is user supplied.

The second line indicates how we estimate the ``-e`` parameter, the baseline error rate. The user can also supply this. 

The rest of the log shows how long certain steps take in each of floria's steps. The "Total time ..." line indicates floria has exited successfully. 

floria's output directory structure
^^^^^^^^^^^^^^^^^^^^^^^^

The output structure is outlined in :ref:`usage-outputs` in detail. Let's go over the important files. 

contig_ploidy_info.tsv
********************

.. code-block:: sh
    
    cat example_output/contig_ploidy_info.tsv
    contig	average_local_ploidy	average_global_ploidy	approximate_coverage_ignoring_indels	total_haplotig_bases_covered	average_local_ploidy_min1hapq	average_global_ploidy_min1hapq
    NZ_CP081897.1	3.000	3.000	79.601	354838	2.984	2.984

This file gives information about how many strains appear to be present in our sample. For each contig, a line with information will be present in this file. We only have one contig, so we only get one (non-header) line. 









