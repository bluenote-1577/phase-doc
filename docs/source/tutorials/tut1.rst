Tutorial 1: running floria and visualizing its outputs
================================================

This minimal tutorial goes over floria's inputs and outputs by running floria on toy data. 

Running floria and requirements
^^^^^^^^^^^^

We will assume that you've installed floria. See the :doc:`../quick_start` page if you have not installed floria. Once floria is installed, run the following commands.

.. code-block:: sh

   # clone the repository if you have not done so
   git clone https://github.com/bluenote-1577/floria
   cd floria

   # assuming you've installed floria
   floria -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output

For visualization, you will need:

#. `IGV <https://software.broadinstitute.org/software/igv/>`_ installed. 
#. matplotlib, numpy, scipy, python3 installed.


Inputs
^^^^^

The provided files represent **a mixture of synthetic Klebsiella Pneumoniae reads** for **three different strains**. 

``test_long.bam`` represent synthetic nanopore reads from three different strains at 70x, 30x, and 10x coverage. We then mapped these reads to the ``MN-03.fa`` reference genome to generate the bam file. We only retained reads between positions 1-120,000 for this example. 

``test.vcf`` was obtained by calling SNPs on ``test_long.bam``. We used longshot, a long-read SNP caller to call these SNPs. 

``MN-03.fa`` is the contig we mapped onto. This is an arbitrary Klebsiella Pneumoniae strain.

Output log
^^^^^^^^^

Depending on your version of floria, you should see something like the following after running the command. 

.. code-block:: sh

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

The "Total time ..." line indicates floria has exited successfully. A few lines that may be of interest:

#. The first line indicates what the ``-l`` parameter, the block length. See 
#. The second line indicates how we estimate the ``-e`` parameter, the baseline error rate. The user can also supply this. 
#. The line "Final SNP error rate ..." describes the fraction of alleles that are discordant in your final outputs. 

For understanding the ``-l`` and ``-e`` parameters, see :doc:`how-to-guides/htg1`. 

floria's main outputs
^^^^^^^^^^^^^^^^^^^^^

The output structure is outlined in :ref:`usage-outputs` in detail, and a more detailed description of each output file is also in the link. Let's go over the important files. 

contig_ploidy_info.tsv
********************

.. code-block:: sh
    
    cat example_output/contig_ploidy_info.tsv

    contig	average_local_ploidy	average_global_ploidy	approximate_coverage_ignoring_indels	total_vartig_bases_covered	average_local_ploidy_min1hapq	average_global_ploidy_min1hapq  avg_err
    NZ_CP081897.1	3.000	3.000	79.601	354838	2.984	2.984   0.055

This file gives information about how many strains appear to be present in our sample. For each contig, a line with information will be present in this file. We only have one contig, so we only get one (non-header) line. 

The most important columns is the ``average_global_ploidy`` line, which gives a suggestion of how many strains there are. In our case, it looks like there are exactly 3 strains present, which is correct. 

As a sanity check, the ``total_vartig_bases_covered`` column gives 354838 bases, which is approximately 3 (the number of strains) times 120,000. Recall that we restricted our reads to lie within 1-120,000 along the contig. We'll see explain what a "vartig" is below. 

NZ_CP081897.1.vartigs
**********************

For the rest of the files, they are stored in ``example_output/NZ_CP081897.1``. floria outputs results for each individual contig in the bam file. Because we only have one contig, there is only one output folder. Let's start with the ``*.vartigs`` file. 









