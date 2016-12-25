Windows
=======

This is the installation instruction for the python version, not the Matlab version.

The easiest way to use psignifit is to use the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.
In that case, you might want to download psignifit from

`Psignifit3 Downloads <http://sourceforge.net/projects/psignifit/files/>`_

Choose one of the files ``psignifit3.0_beta_<date of the snapshot>_win32-py2.6.exe``. Executing the
installer should install psignifit3.0 on your system.

Installing the command line interface (optional)
-------------------------------------------------

Download the file ``psignifit_cli_3_beta_installer_<date of the snapshot>.exe``
from `sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and run
it.  This is a windows installer and can be installed as usual.  Follow the
instructions on the screen. At the end of the installation, you will be asked
whether you want to add psignifit-cli to your environment path. You should leave
this box checked. You will not be able to use psignifit from within Matlab if
you uncheck this box!


Testing your installation
-------------------------

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

