# ogrdbstats 0.5.0

* This is the first packaged version available from CRAN, and also available as a Docker image. The documentation has been reworked, and R help has been added for externally callable functions.


# ogrdbstats 0.5.1

* Fixed a problem which prevented haplotyping plots from being produced.
* Fixed various issues with report formatting

* Added support for recent versions of IgDiscover, which now uses a AIRR format for filtered.tsv
* Added automatic sequence gapping for the AIRR format column sequence_alignment, should it not be gapped (this was required for IgDIscover support)
* Fixed the creation of output files where the input file included more than one period in its name
