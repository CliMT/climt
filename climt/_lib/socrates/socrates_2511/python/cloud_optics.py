# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots cloud optical properties against characteristic dimension
from a fitting file created by the scatter_average program.
'''

import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        raise RuntimeError('please enter a fit file name')

re=np.zeros(4096)
extinction=np.zeros(4096)
scattering=np.zeros(4096)
ss_albedo=np.zeros(4096)
asymmetry=np.zeros(4096)
plotrow=0
n_param=0
scheme=0
with open(filename) as f:
    for line in f:
        if (line.find("Min and max characteristic dimension") == 0):
            cols=line.strip().split()
            min_d=float(cols[5])
            max_d=float(cols[6])
            for n in range(0, 4096):
                re[n]=min_d + float(n)*(max_d-min_d)/4096.0
        if (line.find("Index of parametrization scheme") == 0):
            cols=line.strip().split()
            scheme=int(cols[5])
        if (line.find("Number of fitting parameters") == 0):
            cols=line.strip().split()
            n_param=int(cols[5])
        if (line.find("Band") == 0):
            cols=line.strip().split()
            band=int(cols[2])
            param=np.zeros(n_param)
            for n in range(0, int(n_param/4)):
                cols=f.readline().strip().split()
                param[n*4]=float(cols[0])
                param[n*4+1]=float(cols[1])
                param[n*4+2]=float(cols[2])
                param[n*4+3]=float(cols[3])
            if (plotrow == 0):
                fig, ax = plt.subplots(3,3, sharex=True, figsize=(13, 10))
            if (scheme == 20):
                for n in range(0, 4096):
                    extinction[n] = (
                      ( param[0]+re[n]*(param[1]+re[n]*param[2]) ) /
                      ( 1.0+re[n]*(param[3]+re[n]*(param[4]+re[n]*param[5])) ) )
                    ss_albedo[n] = ( 1.0 -
                      ( param[6]+re[n]*(param[7]+re[n]*param[8]) ) /
                      ( 1.0+re[n]*(param[9]+re[n]*param[10]) ) )
                    scattering[n] = extinction[n] * ss_albedo[n]
                    asymmetry[n] = (
                      ( param[11]+re[n]*(param[12]+re[n]*param[13]) ) /
                      ( 1.0+re[n]*(param[14]+re[n]*param[15]) ) )
                ax[plotrow,0].plot(re*1e6, extinction,
                    linestyle='solid', linewidth=1.0, color='black',
                    label='Extinction (m2/kg)')
                ax[plotrow,0].plot(re*1e6, scattering,
                    linestyle='dotted', linewidth=1.0, color='blue',
                    label='Scattering (m2/kg)')
                ax[plotrow,0].legend()
                ax[plotrow,0].set_title('Band '+str(band))
                ax[plotrow,1].plot(re*1e6, ss_albedo,
                    linestyle='solid', linewidth=1.0, color='blue',
                    label='Single Scattering Albedo')
                ax[plotrow,1].legend()
                ax[plotrow,2].plot(re*1e6, asymmetry,
                    linestyle='solid', linewidth=1.0, color='green',
                    label='Asymmetry')
                ax[plotrow,2].legend()
                plotrow += 1
                if (plotrow == 3):
                    ax[2,0].set_xlabel('Characteristic Dimension (micron)')
                    ax[2,1].set_xlabel('Characteristic Dimension (micron)')
                    ax[2,2].set_xlabel('Characteristic Dimension (micron)')
                    plotrow=0
                    plt.tight_layout()
                    plt.show()
    if (plotrow > 0):
        plt.tight_layout()
        plt.show()
