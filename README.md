bscall
======

Variant Caller for Bisulfite Sequencing Data.


------------
Installation
------------

Before starting the installation of bscall, you should check if your system has the GSL library already installed.

If your system does not have GSL library then you can download it from [GSL](https://www.gnu.org/software/gsl/) and follow the installation steps. 

Once GSL is already available on your system then you can compile and install bscall.

1) Change GSL library paths. In order to compile bscall you must specify the GSL headers and library directories. 
   To do that, edit Gsl.mk file with the proper paths. Just modify two lines starting with GSL_LIB and GSL_INC.

Gsl.mk:

    #1. MODIFY HERE THE GSL LIBRARY LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/lib
    GSL_LIB = -L/path/to/GSL/lib/
    #2. MODIFY HERE THE GSL HEADERS LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/include
    GSL_INC = -I/path/to/GSL/include/ 

2) After editing Gsl.mk just type make all to get the code compiled.

Compile:

    make all

3) Installation. If the compilation process has been successfully completed then a binary file should be found at bin directory. Just copy it to a directory included in your
$PATH.

Copy binary:

    cp bin/bscall /usr/bin

--------------
Running bscall
--------------

Run bscall from a BAM file to get a bcf output:

    samtools view -h my_aligned.bam chr1 | bs_call -r my_reference.fasta -p -L5 -n my_sample_name | bcftools convert -o mysample_chr1.bcf -O b

This command assumes that you have samtools and bcftools already installed on your system.

The parameters configured for this example are -p (Paired End Data) and -L5 (5 bases to trim from left of read pair).


---------
Changelog
---------
2.0.1 Fix Single End Memory Leak
2.0.1 Use of dbSNP to evaluate SNP Calling.
2.0.1 Output JSON stats to measure SNP and Methylation calls.
2.0.1 Index set of Bed files for dbSNP.

---------
Licensing
---------

bscall is licensed under GPL. See LICENSE for more information.

------
Author
------

Simon Heath at (CNAG/CRG) Centre Nacional d’Anàlisi Genòmica / Centre de Regulació Genòmica.
simon.heath@gmail.com

