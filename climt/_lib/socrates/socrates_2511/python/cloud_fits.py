# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots the fitted cloud optical properties against the scaled data points
from a monitoring file created by the scatter_average program.
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if (len(sys.argv) > 2):
        filename = sys.argv[1]
        ndists = int(sys.argv[2])
    elif (len(sys.argv) > 1):
        filename = sys.argv[1]
        ndists = 56
    else:
        raise RuntimeError('please enter a mon file name')

size=np.zeros(ndists)
actual=np.zeros(ndists)
fitted=np.zeros(ndists)
linenum=0
plotnum=0
plotrow=0
with open(filename) as f:
    for line in f:
        linenum += 1
        if (line.find("Fitted") == 0):
            title=line
        if (line.find("Scaled Size") == 4):
            f.readline()
            for n in range(0, ndists):
                cols=f.readline().strip().split()
                size[n]=float(cols[0])
                actual[n]=float(cols[1])
                fitted[n]=float(cols[2])
            if (plotnum == 0):
                if (plotrow == 0):
                    fig, ax = plt.subplots(3,3, sharex=True, figsize=(13, 10))
            ax[plotrow,plotnum].plot(size,actual, 's', label='Actual')
            ax[plotrow,plotnum].plot(size,fitted, '+', label='Fitted')
            ax[plotrow,plotnum].set_title(title)
            plotnum += 1
            if (plotnum == 3):
                plotnum=0
                plotrow += 1
            if (plotrow == 3):
                plotnum=0
                plotrow=0
                plt.legend()
                plt.tight_layout()
                plt.show()
    if (plotrow > 0):
        plt.tight_layout()
        plt.show()
