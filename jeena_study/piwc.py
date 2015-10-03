####TASK3
None
##INITIALIZING
###load('signal5_t3.mat');%contains the variables R and t
##load('Butter.mat');%A,B
fs=44100#Hz
fc=10e+3#Hz
Nsc=128#number of subcarriers in OFDM. Ts is the transmissiontime of Nsc data symbols
Ncp=20 #length of the cyclic prefix
N=Nsc+Ncp #length of OFDM block with cp
Tofdm=58e-3  #symbol time fo 1 OFDM symbol without cp 
#deltaf=1/Tofdm;
#bw=Nsc*deltaf;
SRofdm=1/(Tofdm/Nsc)#Hz sample rate of the OFDM 2206.9 'bw'
#R signal loader, the received signal
#t time interval for signal

#CHANNEL
###Down conversion
import math
cI=math.cos(2*math.pi*fc*t)*math.sqrt(2)  #I carrier
cQ=math.sin(2*math.pi*fc*t)*math.sqrt(2)   #Q carrier
rI=R*cI
rQ=R*(-cQ)
rS=rI+i*rQ
