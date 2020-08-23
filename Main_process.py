#!/usr/bin/env python
import os,sys
import string
import numpy as np
import scipy as scp
import matplotlib
#matplotlib.use('PS')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import cmath
import sync
from Jana_Functions import *

#     _____      __  __  _                 
#    / ___/___  / /_/ /_(_)___  ____ ______
#    \__ \/ _ \/ __/ __/ / __ \/ __ `/ ___/
#   ___/ /  __/ /_/ /_/ / / / / /_/ (__  ) 
#  /____/\___/\__/\__/_/_/ /_/\__, /____/  
#                            /____/        
# settings 
# Naming of the files. 
tracker_data_files = ["SET1_BP.txt", "SET1_SP.txt"]
tracker_longnames  = ["Big Propeller", "Small Propeller"]#, "Base"] #This MUST have the same length as tracker_data_files
length_unit = "mm"

#time that the deceleration occurs, also the when we begin the fft.
#Deceleration = 12.3 #set2
Deceleration = 12.5 #set1


standard_tail_percent = 0.25
color_sequence = ['blueviolet','darkorange','limegreen']

#Switches to turn on and off plots and sections of the analysis
plot2 = True	#Time domain sync
plot3 = True	#X,Y vibrations
plot4 = True    #FFT sync
plot5 = False   #FFT phase
endprints = False #I forget what these do.  

if len(tracker_data_files) != len(tracker_longnames):
	print("ERROR! The lengths of tracker_longnames and tracker_data_files don't match and they have to.")

#     _____ __                               ________               
#    / ___// /_________  ____ _____ ___     / ____/ /___ ___________
#    \__ \/ __/ ___/ _ \/ __ `/ __ `__ \   / /   / / __ `/ ___/ ___/
#   ___/ / /_/ /  /  __/ /_/ / / / / / /  / /___/ / /_/ (__  |__  ) 
#  /____/\__/_/   \___/\__,_/_/ /_/ /_/   \____/_/\__,_/____/____/  
#                                                                   
# stream class

class Stream:
	def __init__(self, filename, long_name, tail_percent,color = 'magenta'):
		self.long_name = long_name
		self.color  = color
		self.i_deceleration = -1
		x = []
		y = []
		z = []
		self.t = []
		file_handle = open(filename,'r')
		n_skipped_lines = 0
		for line in file_handle.readlines(): # Read text file by lines
			words = line.strip().split() # Splitting lines to separate numbers ("words")
			if len(words) < 3:
				n_skipped_lines = n_skipped_lines+1
				continue
			x.append( float(words[0])) 
			y.append( float(words[1])) 
			z.append( float(words[2])) 
			self.t.append( float(words[3])) 
		if n_skipped_lines > 3:
			print("warning: skipping a lot of lines (",n_skipped_lines,") while loading ",filename)
		file_handle.close() # Closing the file
		#normalize coords and save those.
		self.x_norm = np.array(x) - mean_end(x, tail_percent)
		self.y_norm = np.array(y) - mean_end(y, tail_percent)
		self.z_norm = np.array(z) - z[0]

	def offset_time_series(self,offset):
		for i in range(len(self.t)):
			self.t[i] += offset

#      __  ___      __           _____ __                                
#     /  |/  /___ _/ /_____     / ___// /_________  ____ _____ ___  _____
#    / /|_/ / __ `/ //_/ _ \    \__ \/ __/ ___/ _ \/ __ `/ __ `__ \/ ___/
#   / /  / / /_/ / ,< /  __/   ___/ / /_/ /  /  __/ /_/ / / / / / (__  ) 
#  /_/  /_/\__,_/_/|_|\___/   /____/\__/_/   \___/\__,_/_/ /_/ /_/____/  
#                                                                        
print(f"Making Streams")
streams = []

for i in range(len(tracker_data_files)):
	streams.append( Stream( filename = tracker_data_files[i],  long_name = tracker_longnames[i], tail_percent = standard_tail_percent, color = color_sequence[i])) 

if len(tracker_data_files) < 2:
	print("WARNING! Less than 2 trackers entered. We need at least 2 to do this analysis")
	

#      _______           __   _______                   ____  _ ________    
#     / ____(_)___  ____/ /  /_  __(_)___ ___  ___     / __ \(_) __/ __/____
#    / /_  / / __ \/ __  /    / / / / __ `__ \/ _ \   / / / / / /_/ /_/ ___/
#   / __/ / / / / / /_/ /    / / / / / / / / /  __/  / /_/ / / __/ __(__  ) 
#  /_/   /_/_/ /_/\__,_/    /_/ /_/_/ /_/ /_/\___/  /_____/_/_/ /_/ /____/  
#                                                                           
# Time difference constants
print(f"Finding time constants")
Time_offsets = []
for i in range(1,len(streams)):
	#print(f'streams[i].t, {streams[i].t}')
	dt = sync.doSync(streams[0].z_norm, streams[0].t, streams[i].z_norm, streams[i].t) - (streams[0].t[0] - streams[i].t[0])
	Time_offsets.append( dt )
	streams[i].offset_time_series(dt)
	print(f"offset time between ",streams[0].long_name," and ",streams[i].long_name, " is %.6fs"%(dt))

#      ____  __      __     _____                      ______               __              __      __           ________          __ 
#     / __ \/ /___  / /_   /__  /       _   _______   /_  __/  _      _____/ /_  ___  _____/ /__   / /_   ____  / __/ __/_______  / /_
#    / /_/ / / __ \/ __/     / /       | | / / ___/    / /    (_)    / ___/ __ \/ _ \/ ___/ //_/  / __/  / __ \/ /_/ /_/ ___/ _ \/ __/
#   / ____/ / /_/ / /_      / /__      | |/ (__  )    / /    _      / /__/ / / /  __/ /__/ ,<    / /_   / /_/ / __/ __(__  )  __/ /_  
#  /_/   /_/\____/\__/     /____/      |___/____/    /_/    ( )     \___/_/ /_/\___/\___/_/|_|   \__/   \____/_/ /_/ /____/\___/\__/  
#                                                           |/                                                                        
#Plotting the Z components versus Time to check the shift
if plot2:
	print(f"plot z's")
	plt.figure(2, figsize=(14,5))
	for stream in streams:
		#print "Healt of stream ",stream.long_name," t size: ", len(stream.t), " z size: ", len(stream.z_norm)
		plt.plot(stream.t, stream.z_norm, color=stream.color, linewidth=2, label=stream.long_name)
	plt.title('Check of the Z-fit')
	plt.ylabel('Normalized Z ['+length_unit+']')
	plt.xlabel('Time [s]')
	plt.legend()
	plt.show()
	#plt.savefig('plot_Z_vs_t.png',dpi=600)

#      __  ___      _          ___                __           _     
#     /  |/  /___ _(_)___     /   |  ____  ____ _/ /_  _______(_)____
#    / /|_/ / __ `/ / __ \   / /| | / __ \/ __ `/ / / / / ___/ / ___/
#   / /  / / /_/ / / / / /  / ___ |/ / / / /_/ / / /_/ (__  ) (__  ) 
#  /_/  /_/\__,_/_/_/ /_/  /_/  |_/_/ /_/\__,_/_/\__, /____/_/____/  
#                                               /____/               

# finrd the Index of closest point to the deceleration time of the `op
print(f"find i's of deceleration")
for stream in streams:
	stream.i_deceleration = find_nearest(np.array(stream.t),Deceleration)

#FFT's of the measurements from the deceleration leg onward.  #probably should do something with these.
print(f"Do FFT")
for stream in streams:
	stream.res_x = np.fft.rfft(stream.x_norm[stream.i_deceleration: ])
	stream.res_y = np.fft.rfft(stream.y_norm[stream.i_deceleration: ])

print("finding freq's")
for stream in streams:
	sample_frequencies = get_sampling_frequencies(stream.t) 
	stream.freq = get_freq(sample_frequencies, len(stream.x_norm) - stream.i_deceleration)

trunk = []
for stream in streams:
	#vstream = sync.Build_Derivative(stream.z_norm, sync.get_sampling_period(stream.t))
	trunk.append(stream.i_deceleration)
    
for i in range(len(streams)):
		print(streams[i].long_name, "number of data points: ", len(streams[i].z_norm), " trunk i_deceleration: ", trunk[i])

if plot3:
	print(f"plot x,y")
	plt.figure(3, figsize=(14,5))
	plt.subplot(211)
	plt.gca().set_title('X')
	for i in range(len(streams)):
		plt.plot(streams[i].t[trunk[i]:], streams[i].x_norm[trunk[i]:], color=streams[i].color, linewidth=1, label=streams[i].long_name)
	#for stream in streams:
	#	plt.plot(stream.t, stream.x_norm, color=stream.color, linewidth=1, label=stream.long_name)
	plt.ylabel('Normalized X ['+length_unit+']')

	plt.subplot(212)
	plt.gca().set_title('Y')
	for i in range(len(streams)):
		plt.plot(streams[i].t[trunk[i]:], streams[i].y_norm[trunk[i]:], color=streams[i].color, linewidth=1, label=streams[i].long_name)
	#for stream in streams:
	#	plt.plot(stream.t, stream.y_norm, color=stream.color, linewidth=1, label=stream.long_name)
	plt.xlabel('Time [s]')
	plt.ylabel('Normalized Y ['+length_unit+']')
	plt.legend()
	#plt.savefig('plot3.png',dpi=600)
	plt.show()

#      ____  __      __     __________________
#     / __ \/ /___  / /_   / ____/ ____/_  __/
#    / /_/ / / __ \/ __/  / /_  / /_    / /   
#   / ____/ / /_/ / /_   / __/ / __/   / /    
#  /_/   /_/\____/\__/  /_/   /_/     /_/     
#                                             
if plot4:
	print(f"do plot 4")
	plt.figure(4, figsize=(14,5))
	plt.subplot(211)
	plt.gca().set_title('X')
	for stream  in streams:
		plt.semilogy(stream.freq, np.abs(stream.res_x), color=stream.color,linewidth=1, label=stream.long_name)
	plt.ylabel('Amplitude')
	#if which_system == system_3Trackers:
	#	plt.annotate('Cryogenic recovery system (4.3-4.5 Hz)', xy=(4.5, 180), xytext=(7, 300),
	#		arrowprops=dict(facecolor='darkorange', shrink=0.03),
	#		)
	#	plt.annotate(' ', xy=(4.5, 15), xytext=(7, 300),
	#		arrowprops=dict(facecolor='blueviolet', shrink=0.03),
	#		)
	#	plt.annotate('42.12 Hz', xy=(42.12, 3.2), xytext=(45, 50),
	#		arrowprops=dict(facecolor='blueviolet', shrink=0.03),
	#		)
	#	plt.annotate('6.06 Hz', xy=(6.06, 1.5), xytext=(9, 30),
	#		arrowprops=dict(facecolor='darkorange', shrink=0.03),
	#		)
	plt.subplot(212)
	plt.gca().set_title('Y')
	for stream in streams:
		plt.semilogy(stream.freq, np.abs(stream.res_y), color=stream.color,linewidth=1, label=stream.long_name)
	plt.xlabel('Frequency [Hz]')
	plt.ylabel('Amplitude')
	#if which_system == system_3Trackers:
	#	plt.annotate('Cryogenic recovery system', xy=(4.25, 150), xytext=(7, 400),
	#		arrowprops=dict(facecolor='darkorange', shrink=0.03),
	#		)
	#	plt.annotate('6.86 Hz', xy=(6.86, 2.5), xytext=(9, 30),
	#		arrowprops=dict(facecolor='blueviolet', shrink=0.03),
	#		)
	#	plt.annotate('43.20 Hz', xy=(43.2, 0.7), xytext=(45, 2),
	#		arrowprops=dict(facecolor='darkorange', shrink=0.03),
	#		)
	plt.legend()
	#plt.savefig('plot4_FFT.png',dpi=600)
	#plt.show()

#      ____  __      __     __________________         __    _ 
#     / __ \/ /___  / /_   / ____/ ____/_  __/  ____  / /_  (_)
#    / /_/ / / __ \/ __/  / /_  / /_    / /    / __ \/ __ \/ / 
#   / ____/ / /_/ / /_   / __/ / __/   / /    / /_/ / / / / /  
#  /_/   /_/\____/\__/  /_/   /_/     /_/    / .___/_/ /_/_/   
#                                           /_/                
if plot5:
	print(f"do plot 5")
	plt.figure(5, figsize=(14,5))
	plt.subplot(211)
	for stream in streams:
		plt.plot(stream.freq, [cmath.phase(z) for z in stream.res_y],color=stream.color,linewidth=1, label=stream.long_name)
	plt.xlabel('Frequency [Hz]')
	plt.ylabel('Phase')
	plt.subplot(212)
	for stream in streams:
		plt.plot(stream.freq, [cmath.phase(z) for z in stream.res_x], color=stream.color,linewidth=1, label=stream.long_name)
	#kplt.plot(freq_sp, [cmath.phase(z) for z in res_y_sp], color='darkorange',linewidth=1, label="Small propeller")
	#plt.plot(freq_bp, [cmath.phase(z) for z in res_y_bp], color='blueviolet',linewidth=1, label="Big propeller")
	#plt.plot(freq_base, [cmath.phase(z) for z in res_y_base], color='limegreen',linewidth=1, label="Base")
	plt.xlabel('Frequency [Hz]')
	plt.ylabel('Phase')
	#plt.savefig('plot5_phiFFT.png',dpi=600)
	#plt.show()


#End prints
#remind me what these do? 
if endprints:
	#These are used to locate frequency peaks within certain frequency ranges.
	#polyprint prints the resulting lists in one line. 
	#get_peak_freq(fft_output, frequencies, min_search_freq, max_search_freq)
	polyprint([get_peak_freq(stream.res_x, stream.freq, 2, 50) for stream in streams])
	polyprint([get_peak_freq(stream.res_y, stream.freq, 2, 50) for stream in streams])

	polyprint([get_peak_freq(stream.res_x, stream.freq, 0, 10) for stream in streams])
	polyprint([get_peak_freq(stream.res_y, stream.freq, 0, 10) for stream in streams])

	print(get_peak_freq(streams[0].res_y, streams[0].freq, 17, 20)) #big propeller only

	polyprint([get_peak_freq(stream.res_x, stream.freq, 20, 30) for stream in streams])
	polyprint([get_peak_freq(stream.res_y, stream.freq, 20, 30) for stream in streams])

	polyprint([get_peak_freq(stream.res_x, stream.freq, 30, 40) for stream in streams])
	polyprint([get_peak_freq(stream.res_y, stream.freq, 30, 40) for stream in streams])

	polyprint([get_peak_freq(stream.res_x, stream.freq, 40, 50) for stream in streams])
	polyprint([get_peak_freq(stream.res_y, stream.freq, 40, 50) for stream in streams])
