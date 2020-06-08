#
# $Id: cdr.cmake 169709 2012-11-04 16:17:04Z nwhitehorn $
#

# This utility macro extracts all of the arguments given except the
# first, and places them into the variable named var.
#
#   car(var arg1 arg2 ...)
macro(cdr var junk)
  set(${var} ${ARGN})
endmacro(cdr)
