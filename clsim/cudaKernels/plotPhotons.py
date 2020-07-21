from matplotlib import pyplot as plt
import csv

x1 = []
y1 = []
z1 = []
with open('photonsCudaSh.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))
        z1.append(float(row[2]))


x2 = []
y2 = []
z2 = []

with open('photonsCLSh.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))
        z2.append(float(row[2]))



fig, ax = plt.subplots(2,2)
#marker size indicates number of scattering before hit
ax[0,0].scatter(x1,y1,s=z1, color='coral')
ax[1,0].scatter(x2,y2,s=z2, color='mediumaquamarine')

ax[1,1].scatter(x1,y1 ,marker="X" , color='coral')
ax[1,1].scatter(x2,y2,  s=80, facecolors='none', color='mediumaquamarine')

ax[0,1].scatter(x1,y1,s=z1, color='coral')
ax[0,1].scatter(x2,y2,s=z2, color='white')

ax[0,0].set(xlabel='x ', ylabel='y')
ax[1,0].set(xlabel='x', ylabel='y')

ax[0,0].title.set_text('CUDA')
ax[1,0].title.set_text('CL')
ax[1,1].title.set_text('comparison ')
ax[0,1].title.set_text('comparison should be white ')
 
plt.show()
#plt.savefig('photondist.png')
