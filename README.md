exon-capture-phylo
==================

A target recovery pipeline for transcriptome-enabled exon-captures over novel or diverged sample taxa. 

Getting started
---------------

Clone the scripts onto your computer:

    git clone git://github.com/digitase/exon-capture-phylo.git

Take a look at our demo.config configuration template and modify it for your machine.

Download our demo dataset demo_data.tar.gz from:

    [link](http://example.com)

Documentation
-------------

Detailed usage information can be found in exon-capture-phylo_manual.pdf

For additional POD documentation for each Perl script, run:
    
    perldoc <script_name>.pl

Contacts
--------

* Ben Bai, ANU (u5205339@anu.edu.au)
* Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)

License and copyright

Copyleft 2013-2014 by the authors.
This script is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 

TODO
----

* Consider blastx --seg yes --soft_masking T for ortholog detection improvements. Reference:

    Gabriel Moreno-Hagelsieb and Kristen Latimer
    Choosing BLAST options for better detection of orthologs as reciprocal best hits
    Bioinformatics (2008) 24 (3): 319-324
    first published online November 26, 2007 doi:10.1093/bioinformatics/btm585 

* Link demo dataset with readme

* add #!/usr/bin/perl -w to scripts

* add change log
