I3RemoveLargeDT - Reviewed by H. Dembinski 05/29/15
=================================================== 

URL: http://code.icecube.wisc.edu/svn/projects/sim-services/trunk/private/I3RemoveLargeDT.cxx

Repository UUID: 16731396-06f5-0310-8873-f7f720988828

Revision: 133136

Last Changed Author: hdembinski

Last Changed Rev: 133136

Last Changed Date: 2015-05-29 18:48:28 -0400 (Fri, 29 May 2015)

Documentation
*************
The auto-generated page on the web is missing:
http://software.icecube.wisc.edu/simulation_trunk/projects/sim-services/index.html

Basic documentation is available under resources/docs, but not registered by the "i3_project" call in CMakeLists.txt, the latter needs to be changed to::

    i3_project(sim-services
      PYTHON_DIR python
      PYTHON_DEST icecube/sim_services
      DOCS_DIR resources/docs)

Internal code documentation of functions/classes is ok.

I am missing an overall brief summary of what I3RemoveLargeDT is supposed to do.

Code
*************

I3RemoveLargeDT.cxx
+++++++++++++++++++

A central task in this module is to compute the minimum, maximum, and median of a range. I recommend to use the accumulator framework in boost for these tasks:
http://www.boost.org/doc/libs/1_58_0/doc/html/accumulators.html

Accumulators are descriptive, because they are high-level, and they do not require a temporary allocation of a vector for storage (in case of the median that is surprising, but they have an algorithm which approximates the median very well without storing a sorted range, I just glanced over the relevant paper).

Lines 120, 123, 181, 182

  These lines are very long. That makes it hard to read the code side-by-side with this code review. The latter two lines may become more readable, if they are broken at appropriate places.

Lines 129, 130, 167, 168

  Unusual initialization of doubles in "copy constructor" style. For PODs, the C-style assignment initialization is more intuitive, which is equivalent in terms of compiler instructions (and would use the copy constructor anyway if this were a class).

Lines 140, 170

  `BOOST_FOREACH` would also make these loop heads more readable.

Lines 132-136

  This comment is not correct. The median is eventually computed, not the mean.

Tests
+++++
http://software.icecube.wisc.edu/coverage/00_LATEST/sim-services/private/sim-services/I3RemoveLargeDT.cxx.gcov.html

According to the coverage report, this module is barely tested. This is odd, because resources/tests/remove_large_dt_pes.py seems to do reasonable in-situ testing, but it is not listed under "i3_test_scripts" in the CMakeLists.txt. This needs to be changed. I recommend to use wildcards for "i3_test_scripts". Kudos for the "All your base are belong to us" reference :).

- Lines - 1.6 %
- Functions - 16.7 %
- Branches - 0.6 %

Review Response 
+++++++++++++++
**A. Olivas 8/19/2015**

* All suggestions in the Documentation section have been addressed.  The code is integrated into the meta-project docs and there's a page about the module, with a good intro.
* Nearly all of the suggestions in Code have been addressed, with the exception of using boost::accumulators.  This is more of a wishlist item, but worth looking into, so it's a ticket for now and maybe it'll become a feature in a future release.
* The Test coverage is not pretty good. 90.6% line coverage and 91.7% function coverage.
