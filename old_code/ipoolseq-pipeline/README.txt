Installing Dependencies
-----------------------

The file "condaenv.pkgs" defines a conda environment that provides all programs
necessary for running the iPool-Seq analysis pipeline. After installing conda
(see https://conda.io), a conda environment can be created with

  conda create -n ipoolseq --file condaenv.pkgs"

Remember that (as all conda environments), this environment must than be activ-
ated before it can be used by doing

  activate ipoolseq

To prevent python and R from using possibly incompatible package version
installed in your home directory, we recommend you additionally do

  export PYTHONNOUSERSITE=1
  export R_LIBS_USER='-'

Running the Pipeline
--------------------

Download the BAM files belonging to 12 libraries from ERA[1], and store the file
named r<experiment_id>/<library>.bam as data/r.<experiment_id>.<library>/raw.bam.
Then run, for each of the folders in data/

  make data/<folder>/ngm.results.rda data/<folder>/ngm.stats.rda

[1] ftp://ftp.sra.ebi.ac.uk/vol1/ERA112/ERA1125781/bam/

More Information
----------------

See Supplementary Information S1 Text "iPool-Seq Analysis Pipeline".
