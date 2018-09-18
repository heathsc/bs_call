bscall
======

Variant Caller for Bisulfite Sequencing Data.


------------
Installation
------------

Before starting the installation of bscall, you should check if your
system has the GSL library already installed.

If your system does not have GSL library then you can download it from
[GSL](https://www.gnu.org/software/gsl/) and follow the installation
steps. Once GSL is available on your system then you can compile and
install bscall.

Configure:

The compilation should be configured by typing:

	./configure
	
If GSL is not in a standard location on your system then the location
should be specified using the --with-gsl option to configure. For
example, if the installation prefix for the gsl library is /opt/local
(so the libraries can be found in /opt/local/lib and the include
directory gsl in /opt/local/include) then the configuration command
line should be:

	./configure --with-gsl=/opt/local

Compile:

    make all

Install:

If the compilation process has been successfully completed then a
binary file should be found at bin directory. Just copy it to a
directory included in your $PATH.

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

    2.0.2 Document configuration process.
    2.0.1 Fix argument -k about discarded reads that do not form proper pairs.
    2.0.1 Fix Single End Memory Leak.
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

