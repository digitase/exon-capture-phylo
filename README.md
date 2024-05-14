exon-capture-phylo
==================

A target recovery pipeline for transcriptome-enabled exon-captures over novel or diverged sample taxa.

Citation
---------------

Bragg, Jason G., Sally Potter, Ana C. Afonso Silva, Conrad J. Hoskin, Benjamin Y. H. Bai, and Craig Moritz. ‘Phylogenomics of a Rapid Radiation: The Australian Rainbow Skinks’. BMC Evolutionary Biology 18, no. 1 (December 2018): 15. https://doi.org/10.1186/s12862-018-1130-4.

Getting started
---------------

Clone the scripts onto your computer:

    git clone git://github.com/digitase/exon-capture-phylo.git

Take a look at our `demo.config` configuration template.

Download our demo dataset `demo_data.tar.gz` from [link](https://docs.google.com/file/d/0B4t9UUDOrtFnRzlNQVV4bHZ5LWc/edit).

See the "Demo usage" section of `exon-capture-phylo_manual.pdf`

Documentation
-------------

Detailed usage information can be found in `ECP_MANUAL.pdf`

For additional POD documentation for each Perl script, run:
    
    perldoc <script_name>.pl

Authors
--------

* Ben Bai, ANU (u5205339@anu.edu.au)
  * Pipeline design and streamlining

* Dr. Jason Bragg, ANU (jason.bragg@anu.edu.au)
  * Core methodology and original scripts

License and copyright
---------------------

Copyleft 2013-2014 by the authors.

This package is free software; you can redistribute it and/or modify it under the same terms as Perl itself (GNU GENERAL PUBLIC LICENSE).

TODO
----

* Consider blastx --seg yes --soft_masking T for ortholog detection improvements. Reference:

    Gabriel Moreno-Hagelsieb and Kristen Latimer
    Choosing BLAST options for better detection of orthologs as reciprocal best hits
    Bioinformatics (2008) 24 (3): 319-324
    first published online November 26, 2007 doi:10.1093/bioinformatics/btm585 
