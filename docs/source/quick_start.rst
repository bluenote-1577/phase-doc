Quick start
=====

Install
-------

See the github repository TODO for installation instructions. In summary:

#. You can compile floria on linux and mac. 

   *  Only basic software like rust, gcc, cmake are required. 
   
#. We offer a precompiled executable for linux 

   *  Should work for > 95% of linux users (requires AVX-2 SIMD instructions that most modern CPUs have)

#. Conda installation forthcoming 

Basic run
---------
Assuming floria is installed and is in PATH (i.e. typing ``floria -h`` works):

.. code-block:: sh

   # example files are inside the floria repository; clone it
   git clone https://github.com/bluenote-1577/floria
   cd floria
   
   # to install on linux, uncomment and run the command below 
   #cargo install --path . --root ~/.cargo
   
   floria -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output
   
Once floria runs successfully, it should exit without any errors and you should see results in the ``example_output`` folder. See the :doc:`tutorials`, :doc:`how-to-guides`, and :doc:`introduction` pages to learn about how to use floria effectively.  
