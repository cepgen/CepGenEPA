import matplotlib.pyplot as plt
import numpy as np
import sys
import math


plt.rcParams["axes.linewidth"] = 1.8
plt.rcParams["xtick.major.width"] = 1.8
plt.rcParams["xtick.minor.width"] = 1.8
plt.rcParams["ytick.major.width"] = 1.8
plt.rcParams["ytick.minor.width"] = 1.8

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["ytick.labelsize"] = 15

plt.rcParams["legend.fontsize"] = 15

plt.rcParams['legend.title_fontsize'] = 'x-large'



##################################################################

def cs_muon_w_condition_Hamzeh(wvalue):  # Eq.62 of Physics Reports 364 (2002) 359-450
    re = 2.8179403262e-15 * 137.0 / 128.0
    me = 0.510998950e-3
    mmuon = 0.105658
    hbarc2 =  0.389
    alpha2 = (1.0/137.0)*(1.0/137.0)

    # Element-wise calculation of beta using np.where
    beta = np.sqrt(np.where(1.0 - 4.0 * mmuon * mmuon / wvalue**2.0 >= 0.0, 1.0 - 4.0 * mmuon * mmuon / wvalue**2.0, np.nan))

    # Element-wise calculation of cs using np.where
    cs = np.where(wvalue > mmuon, ( 4.0 * np.pi * alpha2 * hbarc2 ) / wvalue**2.0 * (beta) * \
             ( (3.0 - (beta**4.0))/(2.0 * beta) * np.log((1.0 + beta)/(1.0 - beta)) - 2.0 + beta**2.0 ), 0.0) * 1e9

    return cs



##################################################################


def cs_muon_w_condition_muon(wvalue):
    re = 2.8179403262e-15 * 137.0 / 128.0
    me = 0.510998950e-3
    mmuon = 0.105658
    hbarc2 =  0.389
    alpha2 = (1.0/137.0)*(1.0/137.0)

    # Element-wise calculation of beta using np.where
    beta = np.sqrt(np.where(1.0 - 4.0 * mmuon * mmuon / wvalue**2.0 >= 0, 1.0 - 4.0 * mmuon * mmuon / wvalue**2.0, np.nan))

    # Element-wise calculation of cs using np.where
    cs = np.where(wvalue > mmuon, (4.0 * np.pi * alpha2 * hbarc2 ) / wvalue**2.0 *
               ( (1.0 + 4.0*mmuon**2.0/wvalue**2.0 - 8.0*mmuon**4.0/wvalue**4.0) * np.log((1.0+beta)/(1.0-beta)) - beta* (1.0 + 4.0*mmuon**2.0/wvalue**2.0)), 0.0) * 1e9

    return cs


##################################################################


def cs_muon_w_condition_Krzysztof(wvalue):
    re = 2.8179403262e-15 * 137.0 / 128.0
    me = 0.510998950e-3
    mmuon = 0.105658
    hbarc2 = 0.389
    alpha2 = (1.0/137.0)*(1.0/137.0)

    # Element-wise calculation of beta using np.where
    beta = np.sqrt(np.where(1.0 - 4.0 * mmuon * mmuon / wvalue**2.0 >= 0, 1.0 - 4.0 * mmuon * mmuon / wvalue**2.0, np.nan))

    # Element-wise calculation of cs using np.where
    cs = (4.0 * np.pi * hbarc2 * alpha2 ) / wvalue**2.0 * \
         ( (1.0 + 4.0 * mmuon**2.0 / wvalue**2.0 - 8.0 * mmuon**4.0 / wvalue**4.0) * np.log((1.0 + beta) / (1.0 - beta)) -
          beta * (1.0 + 4.0 * mmuon**2.0 / wvalue**2.0) ) * 1e9


    return cs


##################################################################




"""
def cs_DM_w(wvalue):
    re = 2.8179403262e-15 * 137.0 / 128.0
    me = 0.510998950e-3
    mDM = 10.0
    hbarc2 =  0.389
    alpha2 = (1.0/137.0)*(1.0/137.0)

    # Element-wise comparison
    condition = (1.0 - 4.0 * mDM * mDM / wvalue) >= 0
    beta = np.sqrt(np.where(condition, 1.0 - 4.0 * mDM * mDM / wvalue, np.nan))

    # Element-wise calculation of cs
    cs = np.where(wvalue > mDM, (4.0 * np.pi * alpha2 * hbarc2 ) / wvalue* (beta) * \
             (2.0 - beta**2.0 - (1-beta**4.0)/(2.0 * beta)*np.log((1.0+beta)/(1.0-beta))), 0.)

    return cs
"""


def trap_integ(wv, fluxv):
    wmin = np.zeros(len(wv) - 1)
    integ = np.zeros(len(wv) - 1)

    for i in range(len(wv) - 2, -1, -1):
        wvwid = wv[i + 1] - wv[i]
        cs_0 = cs_muon_w_condition_Hamzeh(wv[i])
        cs_1 = cs_muon_w_condition_Hamzeh(wv[i + 1])
#        cs_0 = cs_DM_w(wv[i])
#        cs_1 = cs_DM_w(wv[i + 1])
#        cs_0 = cs_DM_w_old(wv[i])
#        cs_1 = cs_DM_w_old(wv[i + 1])
        traparea = wvwid * 0.5 * (fluxv[i] * cs_0 + fluxv[i + 1] * cs_1)
        wmin[i] = wv[i]
        if i == len(wv) - 2:
            integ[i] = traparea
        else:
            integ[i] = integ[i + 1] + traparea

    nanobarn = 1.e+40

    return wmin, integ   # * 1000000000.0


sys.path.append('./values')




from wgrid_10_100_100 import *

wv = np.array(wvalues[3])
ie = np.array(inel[3])
el = np.array(elas[3])

wv1, int_inel = trap_integ(wv, ie)
wv2, int_el = trap_integ(wv, el)

fig, ax = plt.subplots(figsize = (9.0, 8.0))
ax.set_xlim(10.0, 1000.0)
ax.set_ylim(1.e-3, 10.e2)


inel_label = ('$M_N<$ ${{{:g}}}$ GeV').format(inel[0]) + (' ($Q^2_p<$ ${{{:g}}}$ GeV$^2$)').format(inel[2])
title_label = ('$Q^2_e<$ ${{{:g}}}^{{{:g}}}$ GeV$^2$').format(10,np.log10(inel[1]))
plt.loglog(wv2[:303], int_el[:303], linestyle = 'solid',  linewidth=2,  label = 'Elastic')
plt.loglog(wv1[:303], int_inel[:303], linestyle = 'dotted',  linewidth=2, label = inel_label)

#plt.grid()

plt.legend(title = title_label)

#plt.grid()




# Save the output values in a text file
output_data = np.column_stack((wv2[:303], int_el[:303], int_inel[:303]))
header = 'W_Value Elastic Inelastic'
np.savetxt('output_values_mu.txt', output_data, header=header, fmt='%0.8e', delimiter='\t')



font1 = {'family':'serif','color':'black','size':24}
font2 = {'family':'serif','color':'black','size':24}

plt.xlabel("W$_0$ [GeV]", fontdict=font2)
#plt.ylabel("$\sigma_{\\tau^+\\tau^-}$ (W > W$_0$) [pb]", fontdict=font2)
plt.ylabel(r"$\sigma_{{\rm ep}\to {\rm e}(\gamma\gamma\to\mu^+\mu^-){\rm p}^{(\ast)}}$ (W > W$_0$) [pb]", fontdict = font2)


plt.savefig("cs_muon_MN2_mMin2_q2min_Final_25April_10_100_100.pdf")
plt.savefig("cs_muon_MN2_mMin2_q2min_Final_25April_10_100_100.jpg")



plt.show()


"""
from wgrid_2_4_4_0908 import *

wv = np.array(wvalues[3])
ie = np.array(inel[3])
wv1, int_inel = trap_integ(wv, ie)

inel_label = 'M_N < ' + str(inel[0])
plt.loglog(wv1[:101], int_inel[:101], '-', label = inel_label)
plt.legend(title = title_label)

from wgrid_3_4_4_0908 import *

wv = np.array(wvalues[3])
ie = np.array(inel[3])
wv1, int_inel = trap_integ(wv, ie)

inel_label = 'M_N < ' + str(inel[0])
plt.loglog(wv2[:101], int_inel[:101], '-', label = inel_label)
plt.legend(title = title_label)
"""


