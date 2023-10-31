Utility scripts
==============

floria comes with some utillity scripts that are useful for visualizing results. These scripts can be found in ``floria/scripts`` in floria's github repository: https://github.com/bluenote-1577/floria.

Requirements
------------

All scripts below use python3 and require the following standard dependencies (can be installed by ``pip install ...``)

#. matplotlib
#. numpy
#. natsort (https://pypi.org/project/natsort/)

visualize_vartigs.py
-------------------

Very useful script for quickly visualizing output vartigs. 

.. image:: img/visualize-vartig-example.png
  :width: 400
  :alt: floria visualization example.

Each line is a vartig passing certain thresholds. The top plot is colored by the alternate allele ratio of the vartig, i.e. what fraction of the vartig's alleles are alternate? For a good phasing, you would expect the coverage and allele ratios of vartigs from the same strain to be similar. 

.. code-block:: sh

    #visualize vartigs only if HAPQ >= 15, COV >= 50, LENGTH >= 10000
    python visualize_vartigs.py floria_out_dir/contig1/contig1.vartigs -q 15 -c 50 -l 10000

haplotag_output_dir.py
---------------------

If you want to visualize floria's output phasings in the IGV, use this script. This adds a ``HP:i`` tag to the same BAM file that was input to floria for all reads in contigs that were successfully phased by floria. 

.. code-block:: sh

    floria -b my_bam_file.bam -v some.vcf -r ref.fa 
    # floria_out_dir is output directory

    python haplotag_output_dir.py -d floria_out_dir/ -b my_bam_file.bam -o output_prefix

    #haplotagged bam file with HP:i flags. 
    ls output_prefix.bam

For an example of how to visualize this output bam file, see :doc:`tutorials/tut1`. 
    
haplotag_bam.py  
----------------

This script does the same thing as haplotag_output_dir.py, but you can specify exactly which contig to haplotag. 

