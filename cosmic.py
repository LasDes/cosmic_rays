#Analysis of muon detection from cosmic rays for nuclear physics practical course.
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal
from scipy import optimize
from scipy.optimize import curve_fit
import astropy.constants as const

#open data

start_up1=np.loadtxt('start_up',usecols=(0))
start_up2=np.loadtxt('start_up',usecols=(1))

stop_up1=np.loadtxt('stop_up',usecols=(0))
stop_up2=np.loadtxt('stop_up',usecols=(1))

cal1=np.loadtxt('cal',usecols=(0))
cal2=np.loadtxt('cal',usecols=(1))

plt.figure(0)
plt.plot(start_up1,start_up2)

plt.figure(1)
plt.plot(stop_up1,stop_up2)

plt.figure(2)
plt.plot(cal1,cal2)
plt.xlabel('channel')
plt.ylabel('counts')
#Use above channel vs and known times to get a channel time correlation
peaks=scipy.signal.find_peaks(cal2,height=None,threshold=1,distance=100)
print(peaks[0])
#print(peaks[0][0])
delays=[16,24,32,40,48]
cal1_peaks=cal1[peaks[0]]
#linear fit these data
def lin_func(x,m,c):
    y=m*x+c
    return y
'''
popt,pcov=scipy.optimize.curve_fit(lin_func,delays,cal1_peaks)
print(popt)

time_fit=np.linspace(0,np.max(delays)+10)
fit=lin_func(time_fit,popt[0],popt[1])

plt.figure(3)
plt.scatter(delays,cal1_peaks,s=5)
plt.plot(time_fit,fit,'--')
plt.xlabel('time')
plt.ylabel('channel')
plt.ylim(0,1500)
plt.xlim(0,55)
'''
popt,pcov=scipy.optimize.curve_fit(lin_func,cal1_peaks,delays)
print(popt)

plt.figure(3)
plt.scatter(cal1_peaks,delays,s=5)
plt.xlabel('channel')
plt.ylabel('time')
#plt.ylim(0,1500)
#plt.xlim(0,55)


time_start_up1=start_up1*popt[0]+popt[1]
#I think we used 16ns delay
plt.figure(4)
plt.plot(start_up1,start_up2,color='black')
plt.plot(stop_up1,stop_up2,color='r')

plt.figure(5)
plt.plot(time_start_up1,start_up2)
plt.xlabel('time ns')

time_stop_up1=stop_up1*popt[0]+popt[1]
plt.figure(5)
plt.plot(time_stop_up1,stop_up2)
plt.xlabel('time ns')

plt.show()
