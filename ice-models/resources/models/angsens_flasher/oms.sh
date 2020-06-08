#!/bin/awk -f

BEGIN {
    A=0.68/2; B=0.25; // A*(1+1.5*x-x^3/2)+B*x*(x^2-1)^3
    B=0.3; if(ARGC>1) B=ARGV[1]
    print 2*A
    print A
    print 1.5*A-B
    print 0
    print -A/2+3*B
    print 0
    print -3*B
    print 0
    print B
    print 0
    print 0
    print 0
}

