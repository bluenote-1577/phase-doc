Tutorial 1: running floria and visualizing its outputs
================================================

This minimal tutorial goes over floria's inputs and outputs by running floria on toy data. 

We will assume that you've installed floria. See the :doc:`quick_start` page if you have not installed floria. 

If you've cloned the github repository, try the following command to run floria. 

.. code-block:: sh

   # clone the repository if you have not done so
   git clone https://github.com/bluenote-1577/floria
   cd floria

   # assuming you've installed floria
   floria -b tests/test_long.bam -v tests/test.vcf -r tests/MN-03.fa -o example_output
   ls example_output

Let's examine floria's output in ``example_output``. 

