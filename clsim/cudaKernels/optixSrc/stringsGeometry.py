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


domstrings = [[0,0,0]]


with open(path +'domstrings.csv', ) as csvfile:
    read = csv.reader(csvfile, delimiter=',')
    for row in read:
        domstrings.append(row)

print(len(domstrings))

geom = []
factor = 10.0
resolution = 32

for xyid in domstrings[2:]:
    
    # create a new 'Cylinder'
    cylinder1 = Cylinder()
    cylinder1.Resolution = resolution
    cylinder1.Height =  2470.0
    cylinder1.Radius = 0.0235*factor
    cylinder1.Center = [0,0,0]
 
    transform1 = Transform(Input=cylinder1)
    transform1.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform1.Transform.Rotate = [90.0, 0.0, 0.0]

    # create a new 'Transform'
    transform2 = Transform(Input=transform1)
    transform2.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform2.Transform.Translate = [float(xyid[0] ), float(xyid[1] ), 735.0]
    geom.append(transform2)
 

appendGeometry1 = AppendGeometry(Input=geom)

 
appendGeometry1Display = Show(appendGeometry1, renderView1, 'GeometryRepresentation')


objwriter = POBJWriter(FileName = outpath +'stringsrad10.obj' , Input = appendGeometry1)
objwriter.UpdatePipeline()