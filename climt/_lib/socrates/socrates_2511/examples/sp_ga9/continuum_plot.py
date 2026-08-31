import numpy as np
import matplotlib.pyplot as plt

wn_limit=914

fig=plt.figure()

caviar = np.loadtxt('caviar_s296_plot')
elsey_shine = np.loadtxt('elsey_shine_s296_plot')
mt_ckd = np.loadtxt('mt_ckd3p2_s296_plot')

ax1=fig.add_subplot(211)
ax1.plot(caviar[1:wn_limit,0], caviar[1:wn_limit,1], color='tab:green', label='CAVIAR')
ax1.plot(elsey_shine[1:wn_limit,0], elsey_shine[1:wn_limit,1], color='tab:blue', label='Elsey-Shine')
ax1.plot(mt_ckd[1:wn_limit,0], mt_ckd[1:wn_limit,1], color='tab:red', label='MT_CKD 3.2')
ax1.set_yscale('log')
ax1.set_xlabel('Wavenumber (cm$^{-1}$)')
ax1.set_ylabel('Cs (cm$^2$ molec$^{-1}$ atm$^{-1}$)')
ax1.set_title('Self-continuum 296K')
plt.legend()

caviar = np.loadtxt('caviar_s260_plot')
elsey_shine = np.loadtxt('elsey_shine_s260_plot')
mt_ckd = np.loadtxt('mt_ckd3p2_s260_plot')

ax2=fig.add_subplot(212)
ax2.plot(caviar[1:wn_limit,0], caviar[1:wn_limit,1], color='tab:green', label='CAVIAR')
ax2.plot(elsey_shine[1:wn_limit,0], elsey_shine[1:wn_limit,1], color='tab:blue', label='Elsey-Shine')
ax2.plot(mt_ckd[1:wn_limit,0], mt_ckd[1:wn_limit,1], color='tab:red', label='MT_CKD 3.2')
ax2.set_yscale('log')
ax2.set_xlabel('Wavenumber (cm$^{-1}$)')
ax2.set_ylabel('Cs (cm$^2$ molec$^{-1}$ atm$^{-1}$)')
ax2.set_title('Self-continuum 260K')
plt.legend()

plt.tight_layout()
plt.show()
