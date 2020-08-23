#!/usr/bin/env python
import os,sys
import string
import numpy as np
import scipy as scp
import cmath

 
#   ,gggggggggggggg                                                                                      
#  dP""""""88""""""                                      I8                                              
#  Yb,_    88                                            I8                                              
#   `""    88                                         88888888  gg                                       
#       ggg88gggg                                        I8     ""                                       
#          88   8  gg      gg   ,ggg,,ggg,     ,gggg,    I8     gg     ,ggggg,     ,ggg,,ggg,     ,g,    
#          88      I8      8I  ,8" "8P" "8,   dP"  "Yb   I8     88    dP"  "Y8ggg ,8" "8P" "8,   ,8'8,   
#    gg,   88      I8,    ,8I  I8   8I   8I  i8'        ,I8,    88   i8'    ,8I   I8   8I   8I  ,8'  Yb  
#     "Yb,,8P     ,d8b,  ,d8b,,dP   8I   Yb,,d8,_    _ ,d88b, _,88,_,d8,   ,d8'  ,dP   8I   Yb,,8'_   8) 
#       "Y8P'     8P'"Y88P"`Y88P'   8I   `Y8P""Y8888PP88P""Y888P""Y8P"Y8888P"    8P'   8I   `Y8P' "YY8P8P
#                                                                                                                                                                                          


def mean_end(a,p): 
	# a - array, p - percent of the used points from the end of the array; returns float, bool
	# """Function returning a mean of the end part of an array"""
	# It is going to return mean of the last p proportion of the array and a "it did it correctly" bool
	n = int(p*len(a))
	l = len(a)
	if n > 1 and n <= len(a):
		sum = 0
		for i in a[-n:]:
		    sum = sum+i
		return sum/n
	else:
		print(f"There is something goofy about the n in function mean_end \n n =", n, "Lenght of the list =", l)
	return 0.,False


def find_nearest(a,v):
    # a - array, v - value, returns int
    #finds the index of the element of a that is closes to v
    return (np.abs(a-v)).argmin()


def get_freq(sampling_freq, n_samples):
	#"""Function which gives back an array of frequencies corresponding to results of the FFT"""
        max_freq = 0
        n_out = 0 # the number of expected 
        increment = sampling_freq/n_samples
        if (n_samples % 2) == 0:
                max_freq = sampling_freq/2
                n_out = (n_samples+2)/2
        else:
                n_out = (n_samples+1)/2
                max_freq = (n_out-1)*increment
        return  np.arange(0., max_freq + 0.5*increment, increment)


def get_peak_freq(fft_output, frequencies, min_search_freq, max_search_freq):
       #we take in two lists: the abs of the fft output, and a list of the same length of the corresponding frequencies
       #and a frequency range: [min_search_freq, max_search_freq]
       #return the local maximum of fft_output in this frequency range. 
       if not (len(fft_output) == len(frequencies)):
               print(f"get_peak_freq gets lists of different length: len(fft_output)=",fft_output, " len(frequencies)=",frequencies)
               return -1
       max_output = 0
       output_freq = 0
       for i in range(len(frequencies)):
               if frequencies[i] < min_search_freq:
                       continue
               if frequencies[i] > max_search_freq:
                       break 
               if fft_output[i] > max_output:
                       max_output = fft_output[i]
                       output_freq = frequencies[i]
       return output_freq


def diff(x,t):
	#inputs: array of numbers x, array of times of the same length t; 
	#returns the time derivative dx/dt, and an array of corresponding times. both with length = len(x) - 1
	dx_over_dt=range(0,len(x)-1)
	dt=range(0,len(x)-1)
	if len(x)==len(t):
		for i in range(0,len(x)-1):
			dx=x[i+1]-x[i]
			dt2=t[i+1]-t[i]
			dx_over_dt[i] = dx/dt2
			dt[i]=(t[i]+t[i+1])/2
	else: 
		print(f"The arrays are not the same length!")
	return dx_over_dt, dt

def get_sampling_frequencies(time_stream):
	#determins the sampling frequency of the stream in Hz
	#samp_freq_bp   = len(points_t_bp[100:200])/(points_t_bp[200]-points_t_bp[100])
	n = len(time_stream)
	return float(n)/(time_stream[n-1] - time_stream[0])

def polyprint(x):
	Len = len(x)
	if Len ==0:
		print("")
	elif Len ==1:
		print(x[0])
	elif Len ==2:
		print(x[0], x[1])
	elif Len ==3:
		print(x[0], x[1],x[2])
	elif Len ==4:
		print(x[0], x[1],x[2],x[3])
	else:
		print(f"reqrite polyprint to include ",Len)

def Write_Strem_to_file(stream, filename):
	file_handle = open(filename,'w')
	length = len(stream.t)
	for i in range(length):
		value = (t[i], x_norm[i], y_norm[i], z_norm[i])
		s = str(value)
		file_handle.write(s)
	file_handle.close()
