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

'''
Angular distribution and east-west effect
'''

angles_east=np.array([0,10,30,50,70,90]) #degrees
angles_west=-1*angles_east    #degrees
total_angles=np.concatenate([angles_east[1:],angles_west[::-1]])

print(angles_east[1:])
print(total_angles)

counts_east=np.array([207,215,170,71,51,30])
err_counts_east=np.sqrt(counts_east)

counts_west=np.array([207,197,195,99,36,16])
err_counts_west=np.sqrt(counts_west)

total_counts=np.concatenate([counts_east[1:],counts_west[::-1]])
err_total_counts=np.concatenate([err_counts_east[1:],err_counts_west[::-1]])

delta_theta=6   #degrees error to same signifigance as angle measurements.

time=15*60 #time of each measurement in seconds

'''
Random fit
'''

def angular_func(x,a,b,c):
    y=a*(np.cos(b*x))**2 +c
    return y
popt, pcov = curve_fit(angular_func,total_angles,total_counts,bounds=([207,0,0],[207.1,0.05,50]))

fit_angles=np.linspace(-90,90,1000)
fit_counts=angular_func(fit_angles,popt[0],popt[1],popt[2])
print('para',popt)

plt.figure(6)
#plt.errorbar(angles_east,counts_east,err_counts_east,delta_theta,fmt='o',color='blue',ms=4,label='counts east')
#plt.errorbar(angles_west,counts_west,err_counts_west,delta_theta,fmt='o',color='blue',ms=4,label='counts west')
plt.errorbar(total_angles,total_counts,err_total_counts,delta_theta,fmt='o',ms=4,label='counts')
plt.plot(fit_angles,fit_counts,'--')
plt.xlabel('Zenith angle [$^{\circ}$]')
plt.ylabel('Counts')
#plt.legend()
plt.savefig('plots/angular_distribution.png',dpi=400,bbox_inches='tight')


'''
East-west asymmetry coefficient and error analysis of it
'''
eastonly_angles=angles_east[1:]
westonly_angles=-1*angles_west[1:]

eastonly_counts=counts_east[1:]
westonly_counts=counts_west[1:]

err_eastonly_counts=err_counts_east[1:]
err_westonly_counts=err_counts_west[1:]

#print(err_eastonly_counts,err_westonly_counts)
#print(eastonly_angles,westonly_angles)
x=0.0
y=0.0
err_x_2=0.0
err_y_2=0.0
for i in range(len(eastonly_angles)):
    x+=(westonly_counts[i]-eastonly_counts[i])
    y+=(westonly_counts[i]+eastonly_counts[i])
    err_x_2+=((westonly_counts[i])**2 + (eastonly_counts[i])**2)
    err_y_2+=((westonly_counts[i])**2 + (eastonly_counts[i])**2)

err_x=np.sqrt(err_x_2)
err_y=np.sqrt(err_y_2)

epsilon=x/y
err_epsilon=epsilon*(np.sqrt(((err_x/x)**2)+((err_y/y)**2)))


print('The east-west asymmetry coefficient is: ',(chr(949)),' = ',epsilon,(chr(177)),err_epsilon)




plt.show()
