Mac OSX
=======

This is the installation instruction for the Python version, not the Matlab version.

The easiest way to use psignifit is to use the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.

You will also need gcc. You can check whether your machine already has gcc installed by typing::

>>> which gcc

If this gives you the output::

	gcc not found

you have to download gcc through the Apple Developer Tools. Register for a  developer account (you can use your normal apple account for this and it's free, you don't have to join the developer program) which will allow you to access the developer tools where you want to download Xcode (this is a very large file but as far as we know is the only way of downloading gcc) at the time of writing Xcode 3 is free (and has everything you need) so there is no need to pay for Xcode 4. If you are not running Snow Leopard, you will have to find an older version of Xcode such as 3.1.

For OSX Lion Xcode 4 is free, too. Further information can be found `here <http://jessenoller.com/2011/07/30/quick-pythondeveloper-tips-for-osx-lion/>`.

You will want to download the most recent version of psignifit from::

`Psignifit3 Downloads <http://sourceforge.net/.projects/psignifit/files/>`_


Extract the file by typing::

    unzip psignifit3.0_beta_<date of the snapshot>.zip
    cd psignifit3.0_beta_<date of the snapshot>

where you replace <date of the snapshot> by the date string in the file name. Now simply run::

    python setup.py install


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


Testing your installation
-------------------------

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

