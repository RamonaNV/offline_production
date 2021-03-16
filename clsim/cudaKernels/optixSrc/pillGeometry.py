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


outpath = " ../../../test-data/clsim"
path  = " ../../../build/" 

domstrings = [[0,0,0]]

with open(path +'domstrings.csv', ) as csvfile:
    read = csv.reader(csvfile, delimiter=',')
    for row in read:
        domstrings.append(row)

print(len(domstrings))
print(domstrings[1])

doms = [[0,0,0,0,0]]

with open(path +'doms.csv', 'r' ) as csvfile:
    read = csv.reader(csvfile, delimiter=',')
    for row in read:
        doms.append(row)
 
print(len(doms))
print(doms[1])


factor = 1.0
resolution = 32

geom = []

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 2*resolution
cylinder1.Height = 0.234*factor
cylinder1.Radius = 0.15*factor

# create a new 'Sphere'
sphere2 = Sphere()
sphere2.Center = [0.0, 0.0, 0.117*factor]
sphere2.Radius = 0.15*factor
sphere2.ThetaResolution = resolution
sphere2.PhiResolution = resolution
sphere2.EndPhi = 92.0

# create a new 'Sphere'
sphere1 = Sphere()
sphere1.Center = [0.0, 0.0, -0.117*factor]
sphere1.Radius = 0.15*factor
sphere1.ThetaResolution = resolution
sphere1.PhiResolution = resolution
sphere1.StartPhi = 88.0

# create a new 'Transform'
transform1 = Transform(Input=cylinder1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Rotate = [90.0, 0.0, 0.0]
 
# create a new 'Append Geometry'
appendGeometry1 = AppendGeometry(Input=[sphere2, transform1, sphere1])

    # create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=appendGeometry1)

tessellate1 = Tessellate(Input=extractSurface1)
tessellate1.OutputDimension = 2

extractSurface2 = ExtractSurface(Input=tessellate1)

for dom in doms[2:]:
     
    # create a new 'Transform'
    transform2 = Transform(Input=extractSurface2)
    transform2.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform2.Transform.Translate =  [float(dom[0] ),float(dom[1] ),float(dom[2] )]

    geom.append(transform2)

appendGeometry2 = AppendGeometry(Input=geom)

objwriter = POBJWriter(FileName = outpath+'pillDOMs.obj' , Input = appendGeometry2)
objwriter.UpdatePipeline()
 
 
    
 