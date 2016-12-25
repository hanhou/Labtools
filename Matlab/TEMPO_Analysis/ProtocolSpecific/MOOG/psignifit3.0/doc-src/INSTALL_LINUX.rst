Debian/Ubuntu
=============

If you are using `Debian <http://www.debian.org/>`_, the following packages need to be installed:

* ``make``
* ``gcc``
* ``python``
* ``python-dev``
* ``python-numpy (provides python-numpy-dev)``
* ``python-scipy``
* ``python-matplotlib``
* ``python-sphinx``
* ``doxygen``
* ``python-nose``
* ``swig``

In order to check whether or not you have the packages already installed, type::

    aptitude search make gcc python python-dev python-numpy python-scipy python-matplotlib python-nose

Packages that are installed on your machine are listed with a leading <i>

In order to install missing packages, type::

    sudo aptitude install make gcc python python-dev python-numpy python-scipy python-matplotlib python-nose

You will want to download the most recent version of psignifit from::

`Psignifit3 Downloads <http://sourceforge.net/.projects/psignifit/files/>`_

Extract the file by typing::

    unzip psignifit3.0_beta_<date of the snapshot>.zip
    cd psignifit3.0_beta_<date of the snapshot>

where you replace <date of the snapshot> by the date string in the file name.

System-wide installation
------------------------
On the command line, navigate to the root directory of the psignifit distribution. Now, you can simply type::

    sudo python setup.py install

as root and everything will be installed to the right place.

If you want a special flavor of the Python installation and are familiar with using Python
setup-scripts, you can also use special options for the installation, by
executing the ``setup.py`` script explicitly. An example can be found in
the section `Install into users home directory`_.


Install into users home directory
---------------------------------
If you do not have root/admin rights on your computer the setup routine allows installation into your home-directory.
You may install psignifit locally by typing::

    python setup.py install --home=$HOME

where ``$HOME`` should be automatically replaced by the name of your own home-directory.
This command will install psignifit into ``$HOME/lib/python/psignifit``.
To use psignifit from this path, you will also have to set the ``$PYTHONPATH``
variable. Either you invoke the Python interpreter from the commandline by
calling::

    PYTHONPATH=$HOME/lib/python python

or you set the ``$PYTHONPATH`` variable in your ``.bashrc`` (or equivalent) file
by adding the line::

    export PYTHONPATH=$HOME/lib/python

Yet another option is to set the ``$PYTHONPATH`` variable directly from the
Python interpreter using the ``os`` module.


Installing the command line interface (optional)
------------------------------------------------

Download psignifit from `sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and
extract the compressed file to a folder in your home directory. Navigate into the folder.
You have two installation options. By default, the command line interface will be installed to a
folder called ``bin`` in your home directory. You can change this behavior by editing the
``Makefile``. At the beginning of the ``Makefile``, you find a line::

    CLI_INSTALL=$(HOME)/bin

replace this by e.g. ``/usr/bin/`` for system wide installation.

Once you have the Makefile in your desired shape type::

    make cli-install

If the installation directory is not on your system search path, you may have to add it.
To do so, add::

    export PATH=$PATH:$HOME/bin

to your ``.bashrc`` (if you use bash). If you use zsh, the same line should be in your
``.zshrc.local`` file.

Now, you should be able to call::

    psignifit-mcmc -h
    psignifit-diagnostics -h
    psignifit-bootstrap -h
    psignifit-mapestimate -h

And see some usage messages after each call.


Testing your installation (optional)
------------------------------------

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

NOTE: Currently a couple of tests report failures although they actually pass. This will change in the future.
But for now, don't be too alarmed if tests fail. As long as you can start the tests at all, everything is probably ok.

