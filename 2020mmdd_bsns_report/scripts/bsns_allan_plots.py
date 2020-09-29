import allantools as at
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

file_dict = {
    "30 GHz": 'C:\\Users\\James\\Documents\\GitHub\\postings\\2020mmdd_bsns_report\\data\\data\\BSNS_stab_test_30GHz_final.csv',
    "90 GHz": 'C:\\Users\\James\\Documents\\GitHub\\postings\\2020mmdd_bsns_report\\data\\data\\BSNS_stab_test_90GHz_28aug20201_fixed.csv',
    "220 GHz": 'C:\\Users\\James\\Documents\\GitHub\\postings\\2020mmdd_bsns_report\\data\\data\\BSNS_stab_test_220GHz1.csv',
             }

src_index = [1,1,1]
src_fs = [1.0,0.5,0.5]
keys = list(file_dict.keys())
start_ind = [0,5000*0,0]
T_cal = [0,1000,1000]

for srcind in np.arange(len(keys)):
    #print(file_dict[keys[srcind]])
    i = src_index[srcind]
    data = np.genfromtxt(file_dict[keys[srcind]],
                     skip_header=1,delimiter=',')
    data = data[start_ind[srcind]::,:]

    P = data[:,1+i]/np.mean(data[:,1+i])-1
    t = data[:,0+i]-data[0,0+i]
    T = data[:,2+i]*1000

    inc = src_fs[srcind]
    tint = np.arange(0,t[-2],inc)
    Pint = np.interp(tint,t,P)
    Tint = np.interp(tint, t, T)


    Ttau, Tad, Tade, Tns = at.oadev(Tint/np.mean(Tint)-1,1/inc,'phase','decade')
    Ptau, Pad, Pade, Pns = at.oadev(Pint,1/inc,'phase','decade')

    savedict = {
        'Ttau': Ttau,
        'Tad': Tad,
        'Tade': Tade,
        'Tns': Tns,
        'Ptau': Ptau,
        'Pad': Pad,
        'Pade': Pade,
        'Pns': Pns,
        'Pint': Pint,
        'Tint': Tint,
        'tint': tint,
    }

    savefile = '..\\data\\'+keys[srcind]+'\\bsns_src_stability_final.mat'
    sio.savemat(savefile,savedict)

    plt.figure(1)
    plt.plot(tint, (Tint/np.mean(Tint)-1))

    plt.figure(2)
    ax = plt.gca()
    plt.errorbar(Ttau,Tad,Tade,fmt='--o')
    #plt.errorbar(tau,(ad),(ade),fmt='--o')
    ax.set_yscale('log')
    ax.set_xscale('log')
    #plt.title(keys[srcind])
    plt.grid()

    plt.figure(3)
    plt.plot(tint, Pint)

    plt.figure(4)
    ax = plt.gca()
    plt.errorbar(Ptau,Pad,Pade,fmt='--o')
    #plt.errorbar(tau,(ad),(ade),fmt='--o')
    ax.set_yscale('log')
    ax.set_xscale('log')
    #plt.title(keys[srcind])
    plt.grid()

    plt.figure(5)
    plt.plot(Tint/np.mean(Tint)-1,Pint,'.')

plt.show()
