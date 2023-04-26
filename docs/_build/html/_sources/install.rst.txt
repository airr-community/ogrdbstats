.. _install:

Installation
============

OGRDBstats is available as an R package from `CRAN <https://cran.r-project.org/web/packages/ogrdbstats/index.html>`_ and also fully configured in a 
`Docker image <https://hub.docker.com/r/williamlees/ogrdbstats>`_. Stable and development versions are provided.

Installing from CRAN
--------------------

OGRDBstats requires a recent installation of `Pandoc <https://pandoc.org/>`_. If `Rstudio <https://www.rstudio.com/tags/rstudio-ide/>`_ is installed on your machine, Pandoc will already be installed. 
Otherwise please follow the installation instructions on the Pandoc website, or install Rstudio.

Once Pandoc is installed, please start the R interpreter and install ogrdbstats from CRAN:

.. code-block:: R

    install.packages('ogrdbstats')

Alternatively, you can install the development version from Github:

.. code-block:: R

    install.packages('devtools')
    devtools::install_github('https://github.com/airr-community/ogrdbstats')

Once the package is installed, please download `ogrdbstats.R <https://raw.githubusercontent.com/airr-community/ogrdbstats/master/ogrdbstats.R>`_ to a directory of your choice. OGRDBstats is used from the command line via this script, for example:

.. code-block:: bash

    $ Rscript ogrdbstats.R --help
	
You can put `ogrdbstats.R` wherever you want, and refer to it with its full path when you want to run it. In the documentation we will always refer to it as though
it is in the current directory. Please note that Rscript will exit silently if it cannot find `ogrdbstats.R`.



Docker Image
----------------

The image is available in Docker Hub at `williamlees/ogrdbstats <https://hub.docker.com/r/williamlees/ogrdbstats>`_.

.. code-block:: bash

    # Pull the latest stable version
    $ docker pull williamlees/ogrdbstats:stable
    
    # Pull the latest development build
    $ docker pull williamlees/ogrdbstats:latest
	
To use OGRDBstats from the command line, please use the following command in Linux/Mac:

.. code-block:: bash

    $ docker run -v $(pwd):/scratch williamlees/ogrdbstats:stable run_ogrdbstats
	
Or at the Windows command prompt:
	
.. code-block:: bat

    > docker run -v %cd%:/scratch williamlees/ogrdbstats:stable run_ogrdbstats

This will run ogrdbstats.R in the current directory, in exactly the same way that `Rscript ogrdbstats.R` functions in the above CRAN-based installation. For example you can try:

.. code-block:: bat

    > docker run -v %cd%:/scratch williamlees/ogrdbstats:stable run_ogrdbstats --help



