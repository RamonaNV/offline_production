Fits with 85-string flasher data (all horizontal LEDs at full brightness)
numbers in brackets are sat. llh values

1 (3257.39): from h2-50cm

2 (3272.44): from no hole ice

3 (3254.97): from (3+7)/2, nominal DOM eff.

4 (3261.22): from 0.335*(1+[c]), nominal DOM eff.

5 (3255.05): from (1)*0.3+(2)*0.1+(3)*0.4+(4)*0.2

6 (3252.09): from (1)*0.1+(3)*0.6+(5.0)*0.2+(5)*0.1
   fit to 0.335*(1+x)+P3(x^2)*x^3*(x^2-1)

7 (3255.37): from h2-50cm (7.5), then (6)*0.85+(7.5)*0.15
   fit to 0.335*(1+1.5*x-x^3/2)+P0*x*(x^2-1)^3

8 (3256.25): from (6)*0.85+(7)*0.15, new parameterization (oms.sh)

9 (best): new flasher angular sensitivity at p=0.25


oms.sh [p]: create new flasher ang. sens curve with parameter p.
       Best p=0.30. Error band: 0.20...0.40
