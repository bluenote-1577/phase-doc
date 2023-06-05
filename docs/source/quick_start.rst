Quick start
=====

Install
-------

See the github repository TODO for installation instructions. In summary:

#. You can compile floria on linux and mac. 

   *  Only basic software like rust, gcc, cmake are required. 
   
#. We offer a precompiled binary for linux 

   *  Should work for > 95% of linux users (requires AVX-2 SIMD instructions that most modern CPUs have)

#. Conda installation forthcoming 

Basic run
---------
.. code-block:: sh

   #clone repository for example files
   git clone https://github.com/bluenote-1577/glopp
   glopp -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output
   
Once floria runs successfully, it should exit without any errors and you should see results in the ``example_output`` folder. See the tutorials, guides, and importan
