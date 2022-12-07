# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 19:43:18 2021

@author: GabiK
"""

import lmfit as fit
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.stats import norm

def fit_func(x, intensity1, shift1, std1,intensity2, shift2, std2 ):
    return intensity1*np.exp(-(1/2)*(x-shift1)**2/std1**2) + intensity2*np.exp(-(1/2)*(x-shift2)**2/std2**2)

def fit_func2(x, intensity1, shift1, std1):
    return intensity1*np.exp(-(1/2)*(x-shift1)**2/std1**2)

#set working directory
#path='D:/CountingDark_R/emergence/emergence_split_earlysummer/'
path='C:/Users/krive/Dropbox/PhD/Publications/CountingDark_R/emergence/emergence_split_earlysummer/'
os.chdir(path) #sets the current working directory to path
files=[]
for file in glob("*.txt"):
    files.append(file) #creates a list of txt files in working dir
#print(files)

fh = open(path+'/output/earlysummer_start_days.txt', 'w') #write text file with result
fh2 = open(path+'/output/Eldena_2021_smoothened_y.txt', 'w') #write text file with result

for cfile in files:
    Site=cfile.split("_")[0]
    print(Site)
    year=cfile.split("_")[1].split(".")[0]
#    print(Site+year)
#    print(cfile)
    try:
        data = np.loadtxt(cfile)
    except:
        print("could not load", cfile)
        continue
    #plt.plot(data[:,1])
    
    x=data[:,0]
    y=data[:,1]
    
    plt.close('all')
    
    fig, ax = plt.subplots()
    plt.title(Site + " " + year)
    plt.xlabel("days since January 1")
    plt.ylabel("number of passes / night")
    ax.plot(x,y, 'ko', ms=2)
    y = savgol_filter(y, 51, 2)
    #print(y)
    fh2.write(str(list(y.flatten())))
    print("file saved!")
    pr=np.max(y[0:200])*0.1 #ratio should be revised
    peaks,_ = find_peaks(y,prominence=pr,distance=10)
    ax.plot(x,y, '-', label = 'smoothened data')
#    print(peaks)
    if len(peaks) < 2:
        print(cfile, "less than 2 peaks found")
        
        gmodel = fit.Model(fit_func2)
        
        parameters = fit.Parameters()
        parameters.add('intensity1', value = 2000, min = 1, max = 3000)
        parameters.add('shift1', value = peaks[0], min = peaks[0]-4, max = peaks[0]+4)
        parameters.add('std1', value = 11, min = 5, max = 21)
        
        fitres = gmodel.fit(y[0:200], parameters, x = x[0:200])
        I = fitres.best_values['intensity1']
        s = fitres.best_values['shift1']
        d = fitres.best_values['std1']
        
        ax.plot(x,I*np.exp(-(1/2)*(x-s)**2/d**2),label = 'emergence',color='orange')
        cutoff=norm.ppf(0.99, loc=s, scale=d)
        
        ax.axvline(cutoff, ls = '--')
        
        plt.tight_layout()
        plt.savefig('./plots/'+Site+'_'+year+'.png')
        
        cutoff="{:.0f}".format(cutoff)
        
        fh.write( Site + '\t'+  year + '\t'+ cutoff +'\n')
        continue
#    for i in peaks:
#        #print(i)
#        if i <200:
#            ax.axvline(i, ls = '--') #plot dashed lines where peaks are
    if peaks[1]>200:
        print('no early summer activity', cfile)
        fh.write( Site + '\t'+  year + '\t'+ 'NA' +'\n')
        continue
    

    gmodel = fit.Model(fit_func)
    
     
    parameters = fit.Parameters()
    parameters.add('intensity1', value = 2000, min = 1, max = 3000)
    parameters.add('shift1', value = peaks[0], min = peaks[0]-4, max = peaks[0]+4)
    parameters.add('std1', value = 11, min = 5, max = 21)
    
    parameters.add('intensity2', value = 2000, min = 1, max = 3000)
    parameters.add('shift2', value = peaks[1],min = peaks[1]-4, max = peaks[1]+4)
    parameters.add('std2', value = 11, min = 5, max = 21)
    
    #fitres = gmodel.fit(y[0:120], parameters, x = x[0:120])
    fitres = gmodel.fit(y[0:200], parameters, x = x[0:200])
    I = fitres.best_values['intensity1']
    s = fitres.best_values['shift1']
    d = fitres.best_values['std1']
    
    I2 = fitres.best_values['intensity2']
    s2 = fitres.best_values['shift2']
    d2 = fitres.best_values['std2']
    print(I,s,d, I2, s2, d2)

    
    ax.plot(x,I*np.exp(-(1/2)*(x-s)**2/d**2),label = 'emergence',color='orange')
    ax.plot(x,I2*np.exp(-(1/2)*(x-s2)**2/d2**2), label = 'early summer',color='green')
    ax.legend()
    
    def solve_gaussian_crossing(a1,a2,s1,s2,x1,x2):
        a = s2**2 -s1**2
        b = 2*(x2*s1**2-x1*s2**2)
        c = s2**2*x1**2-s1**2*x2**2  -2*s1**2*s2**2*(np.log(a1)-np.log(a2))
        D = b**2-4*a*c
        if D < 0:
            print('no solution')
            return []
        if D == 0:
            print('only one solution')
            return -b/(2*a)
        sol1 = (-b + np.sqrt(D))/(2*a)
        sol2 = (-b - np.sqrt(D))/(2*a)
        #print(sol1)
        #print(sol2)
        return np.array([sol1, sol2])

    result=solve_gaussian_crossing(I,I2,d,d2,s,s2)
    print(result)
    
    ind = np.where((result>0) & (result<200)) #arbitrary, should be revised
    print(result[ind][0])
    ax.axvline(result[ind][0], ls = '--')

    if len(ind) == 0:
        print(cfile, "no solutions between 0 and 200")
        fh.write( Site + '\t'+  year + '\t'+ 'NA' +'\n')
        continue
    if len(ind)>1:
        print(cfile, "more than one solution between 0 and 200")
        print(result)
#        continue
    try:
        date="{:.0f}".format(result[ind][0])
    except:
        print(cfile, ' date string exception')
        fh.write( Site + '\t'+  year + '\t'+ 'NA' +'\n')
        continue
    
    plt.tight_layout()
    plt.savefig('./plots/'+Site+'_'+year+'_emergence_split.png')
    
    fh.write( Site + '\t'+  year + '\t'+ date +'\n')
fh.close()