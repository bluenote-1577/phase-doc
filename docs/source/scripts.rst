Utility scripts
==============

floria comes with some utillity scripts that are useful for visualizing results. These scripts can be found in ``floria/scripts`` in floria's github repository: https://github.com/bluenote-1577/floria.

Requirements
^^^^^^^^^^^

All scripts below use python3 and require the following standard dependencies (can be installed by ``pip install ...``)

#. matplotlib
#. numpy
#. natsort (https://pypi.org/project/natsort/)

visualize_vartigs.py
-------------------

Very useful script for quickly visualizing output vartigs. 

.. image:: ../img/tut1_vartigs.png
  :width: 400
  :alt: floria visualization example.

.. code-block:: sh

    #visualize vartigs only if HAPQ >= 15, COV >= 50, LENGTH >= 10000
    python visualize_vartigs.py floria_out_dir/contig1/contig1.vartigs -q 15 -c 50 -l 10000
    
haplotag_bam.py  
----------------

haplotag_output_dir.py
---------------------

write_contig_headers_vcf.py
-------------
