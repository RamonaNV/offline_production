.. 
.. copyright  (C) 2010
.. The Icecube Collaboration
.. 
.. $Id: faq.rst 94948 2012-11-04 16:21:52Z nwhitehorn $
.. 
.. @version $Revision: 94948 $
.. @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
.. @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: nwhitehorn $


Frequently asked questions and common pitfalls
===============================================

1. I get weird errors like "error: macro "BOOST_PP_SEQ_ELEM_III" requires 2
arguments, but only 1 given" when compiling tableio.

    This happens when you use versions of dataclasses older than tableio itself. Use
    at least version 10-06-00 (part of offline-software 10-06-00). As a general rule,
    you should use tableio with the offline-software release it's part of.