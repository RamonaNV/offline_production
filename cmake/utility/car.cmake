#
# $Id: car.cmake 169709 2012-11-04 16:17:04Z nwhitehorn $
#

# This utility macro extracts the first argument from the list of
# arguments given, and places it into the variable named var.
#
#   car(var arg1 arg2 ...)
macro(car var)
  set(${var} ${ARGV1})
endmacro(car)
