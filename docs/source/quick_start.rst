Quick start
=====

Install
-------

See the github repository https://github.com/bluenote-1577/floria for installation instructions. 

In summary:

#. You can compile floria on linux and mac. 

   *  Only basic software like rust, gcc, cmake are required. 
   
#. We offer a precompiled executable for linux 

   *  Should work for > 99% of linux users 

#. Conda 

  * Floria is also available through bioconda. 

Basic run
---------
If you have `rust <https://www.rust-lang.org/tools/install>`_ installed (version > 1.63), basic linux tools present (gcc, g++, cmake), and are running linux, the following should work. Otherwise, check the github repository for more instructions. 

.. code-block:: sh

   # example files are inside the floria repository; clone it
   git clone https://github.com/bluenote-1577/floria
   cd floria
   
   # see github for installation instructions if you're not on linux/don't have rust installed
   cargo install --path . --root ~/.cargo
   floria -h

   # you can also download a binary if you're on a linux x86-64 architecture
   wget https://github.com/bluenote-1577/floria/releases/download/latest/floria
   chmod +x floria 
   ./floria -h
   
If you've cloned the github repository, try the following command to run floria. 

.. code-block:: sh

   # clone the repository if you have not done so
   #git clone https://github.com/bluenote-1577/floria
   #cd floria

   floria -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output

   #visualize phasings, see Utility scripts on the sidebar for installation information
   #pip install natsort matplotlib numpy
   python scripts/visualize_vartigs.py example_output/NZ_CP081897.1/NZ_CP081897.1.vartigs
   
Once floria runs successfully, it should exit without any errors and you should see results in the ``example_output`` folder. 

To learn more about how floria works, start with the :doc:`introduction` page and the :doc:`tutorials`, in particular, the first tutorial :doc:`tutorials/tut1` picks up from this page. For a deeper look at using floria effectively, see :doc:`how-to-guides`.


