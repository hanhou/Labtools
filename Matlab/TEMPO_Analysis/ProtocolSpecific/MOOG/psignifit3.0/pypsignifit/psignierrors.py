#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

class NosamplesError ( Exception ):
    """An exception that is raised whenever we try to use samples but there are none"""
    def __init__ ( self, msg ):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class SamplingError ( Exception ):
    """An exception that is raised if some thing is wrong with the sampling"""
    def __init__ ( self, msg ):
        self.msg = msg
    def __str__ ( self ):
        return repr(self.msg)
