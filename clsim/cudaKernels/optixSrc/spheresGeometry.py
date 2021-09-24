''' The MIT License (MIT)

Copyright (c) 2021, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
 
#### import the simple module from the paraview
from paraview.simple import *

import numpy as np
import csv

path  = " ../../../build/" 
outpath = " ../../../test-data/clsim"


doms = [[0,0,0,0,0]]

with open(path +'doms.csv', 'r' ) as csvfile:
    read = csv.reader(csvfile, delimiter=',')
    for row in read:
        doms.append(row)
 
print(len(doms))
print(doms[1])

geom = []

for dom in doms[2:]:
   
    sphere1 = Sphere()
    sphere1.ThetaResolution = 32
    sphere1.PhiResolution = 32
    sphere1.Radius = 0.1651
    sphere1.Center = [float(dom[0] ),float(dom[1] ),float(dom[2] )]
    geom.append(sphere1)
 
appendGeometry1 = AppendGeometry(Input=geom)
    
#show data from appendGeometry1
#appendGeometry1Display = Show(appendGeometry1, renderView1, 'GeometryRepresentation')

objwriter = POBJWriter(FileName = outpath+'doms.obj' , Input = appendGeometry1)
objwriter.UpdatePipeline()