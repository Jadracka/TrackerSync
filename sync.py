#!/usr/bin/env python
import string
from math import sin, cos, floor
import random 
import numpy as np
import scipy as scp
import matplotlib
#matplotlib.use('PS')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cmath

"""
  Here we assume some physical function f(t) that has been measured in two streams 
  streamA and streamB, both with equal sampling period T, so that for every sample, indexed by int i,j: 
  streamA[i] = f(i*T) + small_noise 			(i = 0..sizeA-1)
  streamB[j] = f(j*T + t_offset) + small_noise 		(j = 0..sizeB-1)
  streamA[i] ~= streamB[i-t_offset/T] where the two streams overlap
  We want to the offset t_offset, which parameterizes how out of sync the stream are. 
  streamA[i] ~= streamB[i-t_offset/T] where the two streams overlap
   0123456789ABCDEF  #A
    0123456789 	   #B aligned to A with offset = 1
  0123456789         #B aligned to A with offset = -1
"""

"""	#Assumes the same sampling frequency
	float SyncStreams1(float streamA[], float streamB[], float sampling_period) 

	#allows for different sampling frequencies. 
	float SyncStreams2(float streamA[], int sizeA, float streamB[], int sizeB, float sampling_periodA, float sampling_periodB)

	#supporting functions
	float* Build_interpolated_stream(float stream_old[], int size_old, float sampling_period_old, float sampling_period_new, int& newsize)
	enum find_offset_mode{automatic, hint, exhaustive}
	int Find_offset(float streamA[], int sizeA, float streamB[], int sizeB, find_offset_mode hintmode = automatic, int hint_val = 0)
	double GetChi2perDOF(float streamA[], int sizeA, float streamB[], int sizeB, int offset)
	double GetChi2perDOF2(float streamA[], int sizeA, float streamB[], int sizeB, int offset, float STOF)
	float Find_STOF(float streamA[], int sizeA, float streamB[], int sizeB, int offset)

	#Test Drivers.
"""

#####################################################
def Build_Derivative(stream, sampling_period):
	#returns a time derivative of stream, which starts exactly 1/2 a sampling_period after steram
	#and has 1 fewer sample
	size = len(stream) - 1
	if size < 1:
		return stream_v
	stream_v = [0]*size
	for inew in range(1,size):
		stream_v[inew] = (stream[inew]-stream[inew-1])/sampling_period
	return stream_v

def Find_Start_of_Motion(stream_v, thresh):
	#returns the index of stream_v for which the derivative has been over thresh for 3 measurements.
	#if it fails to find it, return -1
	for i in range(len(stream_v)):
		if abs(stream_v[i]) > thresh:
			count_over_thresh += 1
			if count_over_thresh >= 3:
				return i
		else:
			count_over_thresh = 0
	return -1

def Find_End_of_Motion(stream_v, thresh):
	#returns the index of stream_v for which the derivative has been over thresh for 3 measurements.
	#if it fails to find it, return -1
	#for i in range(len(stream_v)):
	count_over_thresh = 0
	for i in range(len(stream_v)-1,-1,-1):
		if abs(stream_v[i]) > thresh:
			count_over_thresh += 1
			if count_over_thresh >= 3:
				return i
		else:
			count_over_thresh = 0
	return -1

def Make_Motion_Stream(stream, sampling_period, thresh):
	#returns the section of stream that includes motion over thresh, 
	#and the time between the beginning of stream and the beginning of the new stream.
	stream_v = Build_Derivative(stream, sampling_period)
	begin = 1 + Find_Start_of_Motion(stream_v, thresh)
	end   = 1 + Find_End_of_Motion(stream_v, thresh) #last one to be included 
	print(f"Make_Motion_Stream returns from %.3fs to %.3fs"%(begin*sampling_period, (end-1)*sampling_period))
	#print(f"Make_Motion_Stream returns from ",begin*sampling_period, "s to ", (end-1)*sampling_period,"s")
	return stream[begin:end], begin*sampling_period #stream begin through end-1 incluseively 

def Build_interpolated_stream(stream_old, size_old, sampling_period_old, sampling_period_new): 
	#float* Build_interpolated_stream(float stream_old[], int size_old, float sampling_period_old, float sampling_period_new, int& newsize):
	#Build a new stream interpolated from stream_old with new sampling period sampling_period_new. 
	#Returns newsize as second output
	#Rules: stream_out[0] = stream_old[0]
	#The old stream must contain at least 2 elements, 
	#usually used to interpolate finer(smaller) sampling periods into larger ones. 

	#need at least 2 elements for interpolation. and don't allow nonsensical sampling periods.	
	if size_old <= 1  or sampling_period_old <= 0. or sampling_period_new <= 0.: 
		print(f"Error Building interpolated stream")
		return stream_old, size_old

	sampling_ratio = sampling_period_new/sampling_period_old #float #expect > 1 
	newsize = int(floor(size_old/sampling_ratio ))#figure out how many elements there are going to be. and set #test
	stream_out = [0] * newsize #wad

	stream_out[0] = stream_old[0]
	i_old_prev = 0 #int

	for inew in range(1,newsize): #for(int inew = 1; inew < newsize; ++inew){ #loop over the cells to be filled in the new stream.
		t_new = inew*sampling_period_new
		i_old = int(floor(t_new/sampling_period_old)) #you moved this
		t_new_fpart = t_new - i_old*sampling_period_old #float
		i_new_fpart = t_new_fpart/sampling_period_old

		if i_old + 1 >= size_old: #protect against overflowing stream_old.
			#TODO: check this so that you never run into this problem.
			print(f"Something's wrong with the sizeing in Build_interpolated_stream, trying to access an iold+1 = ",i_old+1)
			print(f"   you overestimated newize by ",newsize - inew," = ",newsize," - ", inew )
			print(f"   Build_interpolated_stream should continue to funciton despite this error.")
			newsize = inew-1 #set the newsize to the correct value.
			return stream_out

		stream_out[inew] = (1. - i_new_fpart)*stream_old[i_old] + i_new_fpart*stream_old[i_old+1]
	#end for
	return stream_out, newsize	
#end Build_interpolated_stream

def GetChi2perDOF(streamA, sizeA, streamB, sizeB, offset): 
	#double GetChi2perDOF(float streamA[], int sizeA, float streamB[], int sizeB, int offset):
	#returns Chi2 divided by number of overlapping cells
	#if zero overlap, return double max value.

	#Guard against strings missing each other entirely by just giving a huge chi2
	if offset > sizeB or (-offset) > sizeA or sizeA <= 1 or sizeB <= 1:
		"""if offset > sizeB:
			print(f"Warning! offset > sizeB", offset, sizeB)
		elif (-offset) > sizeA:
			print(f"Warning! (-offset) > sizeA",(-offset), sizeA)
		elif sizeA <= 1:
			print(f"Warning! sizeA <=1",sizeA)
		else:
			print(f"Warning! sizeB <=1",sizeB) """
		return 1e30 #std::numeric_limits<double>::max()

	Astart = max(0,offset) #int
	Aend   = min(sizeA,sizeB+offset) #int
	chi2 = 0.
	for iA in range(int(Astart),int(Aend)): #for(int iA = Astart; iA < Aend; iA++)
#		print(f'This is a test! This is a test! {iA-offset}')#{streamA[iA]},{streamB[iA-offset]}')
		chi2 += pow(streamA[iA] - streamB[iA-offset],2)
	return chi2/ float(Aend - Astart)
#end GetChi2perDOF

def GetChi2perDOF2(streamA, sizeA, streamB, sizeB, offset, STOF): 
	#double GetChi2perDOF2(float streamA[], int sizeA, float streamB[], int sizeB, int offset, float STOF):
	#returns Chi2 divided by number of overlapping cells with interpolation for sub-bin offsets.
	#if zero overlap, return double max value or other nonsensical input
	#STOF is in [-1..1]

	#Guard against strings missing each other and out of bounds stof by giving a huge chi2
	if offset > sizeB or (-offset) > sizeA or abs(STOF)>1. or sizeA < 2 or sizeB <2: 
		return 1e30 #std::numeric_limits<double>::max()
	Astart = 0 #int
	Aend = 0  #int
	chi2 = 0. #double
	if offset>=0:# need to access one streamB one element to the left of usual position
		Astart = max(0,offset+1) 
		Aend   = min(sizeA,sizeB+offset)
		for iA in range(Astart,Aend): #for(int iA = Astart; iA < Aend; iA++){
			interp = streamB[iA-offset]*(1. - STOF) + streamB[iA -1 -offset]*(STOF) #float
			chi2 += pow(streamA[iA] - interp,2)
	else: #offset<0, need to access streamB one element to the right of the usual position.
		Astart = max(0,offset) 
		Aend   = min(sizeA,sizeB+offset-1)
		for iA in range(Astart,Aend): #for(int iA = Astart; iA < Aend; iA++){
			interp = streamB[iA-offset]*(1. + STOF) + streamB[iA +1 -offset]*(-STOF) #float
			chi2 += pow(streamA[iA] - interp,2)
	return chi2/float(Aend - Astart)
#end GetChi2perDOF2

def plot_diff(streamA, sizeA, streamB, sizeB, offset_, sampling_period, chi2):
	#makes a plot of the difference between streamA and streamB after 
	#duck1
	offset_ /= sampling_period
	offset = int(floor(offset_))
	STOF = float(offset_ - offset)

	Astart = int(max(0,offset+1)) #int 
	Aend   = int(min(sizeA,sizeB+offset))
	tA = [0]
	for i in range(Aend -1 -Astart): 
		tA.append(tA[i]+sampling_period)

	diff = [0]*(Aend-Astart)
	for iA in range(Aend-Astart): 
		interp = streamB[iA + Astart -offset]*(1. - STOF) + streamB[iA + Astart -1 -offset]*(STOF) #float
		diff[iA] = streamA[iA + Astart] - interp

#	"""
	plt.figure(4,figsize=(14,5))
	plt.plot(tA, diff, color='darkorange',linewidth=2, label="")
	plt.title('Difference between synced streams')
	plt.ylabel('Normalized Z (mm)')
	plt.xlabel('Time (s)')
	#plt.legend()
	plt.show()
#	"""

automatic, hint, exhaustive = (1,2,3) #enum find_offset_mode{automatic, hint, exhaustive}
def Find_offset(streamA, sizeA, streamB, sizeB, hintmode = automatic, hint_val = 0): #tested, passes test3
	#int Find_offset(float streamA[], int sizeA, float streamB[], int sizeB, find_offset_mode hintmode, int hint_val):
	#returns the integer offset needed to alingn the streams to the nearest sample. 
	#This is a very dumb function that runs super slow. 
	#hintmode = automatic, ignore hint and run on automatic, assuming there's some vague attempt to align the streams
	#hintmode = hint, 		manually set a hint that should be within 10% of the correct offset (10% of the shorter stream)
	#hintmode = exhaustive, scan absolute entire range. 

	#default range: have 1/3 overlap at extreme ends of the scan.
	#under automatic
	if sizeA <= 1 or sizeB <= 1: 
		print(f"Error! Bad sizes given to Find_offset. A,B:",sizeA,sizeB)
	offset_scan_start = -sizeA + int(sizeB/3) #was int in python 2, have to force it to be an int in python3
	offset_scan_end = int((2*sizeB)/3) #was int in python 2, have to force it to be an int in python3

	if hintmode == hint:
		aSmidgeon = min(sizeA,sizeB)/10 #10% #int
		offset_scan_start = max(hint_val - aSmidgeon, offset_scan_start)
		offset_scan_end = min(hint_val + aSmidgeon, offset_scan_end)
	elif hintmode == exhaustive:
		offset_scan_start = -sizeA + 1
		offset_scan_end = sizeB-1

	min_chi2perDOF = 1e30 # double #std::numeric_limits<double>::max()
	min_offset = 0 #int
	#duck this loop runs a lot, greatest needed place for improvement. 
	length = offset_scan_end + 1 - offset_scan_start
	#for offset in range(offset_scan_start,offset_scan_end+1):  
	offset = offset_scan_start
	ticks_since_last_change = 0
	last_inc = 1
	mod_1000 = floor(float(offset)/1000.)
	while offset <= offset_scan_end:
		if floor(float(offset)/1000.) > mod_1000: 
			mod_1000 = floor(offset/1000)
			print(f"Offset loop at %.1f%s. Running min chi2: %.2f"%(100.*float(offset - offset_scan_start)/float(length),'%', min_chi2perDOF))
		#do the deed	
		this_chi2perDOF = GetChi2perDOF(streamA, sizeA, streamB, sizeB, offset) #double #receives sizeA==0!
		if this_chi2perDOF < min_chi2perDOF: 
			min_chi2perDOF = this_chi2perDOF
			min_offset = offset
			ticks_since_last_change = 0
		else:
			ticks_since_last_change += last_inc
		#fast increment:
		if this_chi2perDOF > 30:
				last_inc = 20
		elif this_chi2perDOF > 25:
				last_inc = 10
		elif this_chi2perDOF > 20:
				last_inc = 5
		elif this_chi2perDOF < 0.00001:
			if ticks_since_last_change > 200:
				last_inc = 20
			elif ticks_since_last_change > 50:
				last_inc = 10
			elif ticks_since_last_change > 10:
				last_inc = 2
		else:
			last_inc = 1
		offset += last_inc

	if min_offset == offset_scan_start or min_offset == offset_scan_end:
		print(f"Warning! min_offset lands on the edge of of overlap range. Try running with hint = -2 \n")

	chi2perDOFm1 = GetChi2perDOF(streamA, sizeA, streamB, sizeB, min_offset-1) 
	chi2perDOF = GetChi2perDOF(streamA, sizeA, streamB, sizeB, min_offset)
	chi2perDOFp1 = GetChi2perDOF(streamA, sizeA, streamB, sizeB, min_offset+1)
	return min_offset
#end Find_offset

def Find_STOF(streamA, sizeA, streamB, sizeB, offset): 
	#float Find_STOF(float streamA[], int sizeA, float streamB[], int sizeB, int offset):
	#locatest sub-time-step offset fraction. will lie on (-1,1) by binary search by semi-binary search.
	#bounds of the STOF scan, starting at their maximum values of +-1.
	print(f"Doing Find_STOF")
	STOFmax = 2. #float
	STOFmin = -2. #float
	N_interations = 30 #how far we're going with this optomization.  #int

	while N_interations > 0: #loop N_interations times.
		N_interations -= 1
		#query the chi2 for points 1/4 and 3/4 between STOFmax and STOFmin
		probe_left  = max(-1.,0.75*STOFmin + 0.25*STOFmax)  #float
		probe_right = min( 1.,0.75*STOFmax + 0.25*STOFmin) #float
		chi2_left  = GetChi2perDOF2(streamA, sizeA, streamB, sizeB, offset, probe_left) #float
		chi2_right = GetChi2perDOF2(streamA, sizeA, streamB, sizeB, offset, probe_right) #float
		#print "i=",N_interations, "probe left: ", probe_left, " chi2_left: ", chi2_left, " probe_right: ", probe_right, " chi2_right: ", chi2_right

		if chi2_left < chi2_right: #go left
			#print "go left"
			STOFmax = probe_right 
		else: #go right
			#print "go right"
			STOFmin = probe_left  
	#end while
	
	chi2_left  = GetChi2perDOF2(streamA, sizeA, streamB, sizeB, offset, STOFmin) #float
	chi2_right = GetChi2perDOF2(streamA, sizeA, streamB, sizeB, offset, STOFmax) #float

	if chi2_left < chi2_right:
		return STOFmin, chi2_left
	else:
		return STOFmax, chi2_right
#end FIND_STOF

#def SyncStreams1(streamA, sizeA, streamB, sizeB, sampling_period): #wad.
def SyncStreams1(streamA, streamB, sampling_period, zero_suppress=True): #wad.
	#float SyncStreams(float streamA[], int sizeA, float streamB[], int sizeB, float sampling_period): #converted.
	#returns dt such that A(t) = B(t+dt), 
	#where A(t) and B(t) are their time functions with their own clock, so the streams both start at t=0
	#so if streamA starts at t = 0, the start of streamB is at dt
	#returns deltat, chi2
	sizeA = 0
	sizeB = 0
	streamAm=[]
	streamBm=[]
	dt_Am = 0
	dt_Bm = 0
	if zero_suppress:
		streamAm, dt_Am = Make_Motion_Stream(streamA,sampling_period,0.1)
		streamBm, dt_Bm = Make_Motion_Stream(streamB,sampling_period,0.1)
		sizeA = len(streamAm)
		sizeB = len(streamBm)
	else:
		streamAm = streamA
		streamBm = streamB
		sizeA = len(streamA)
		sizeB = len(streamB)

	print(f"Syncing streams: find_offset")
	offset = Find_offset(streamAm, sizeA, streamBm, sizeB) #int
	print(f"course offset", offset, "found, now find fine offset")
	STOF, chi2 = Find_STOF(streamAm, sizeA, streamBm, sizeB, offset) #STOF is sub-time-step offset fraction. will lie on (-1,1). #float
	#full_time_offset = sampling_period*(STOF + offset) + (dt_Am - dt_Bm) 
	full_time_offset = sampling_period*(STOF + offset + 1) + (dt_Am - dt_Bm)  #duck1 sometimes this works better. Seems like there's an intermittent off-by-1 error.
	plot_diff(streamA, len(streamA), streamB, len(streamB), full_time_offset, sampling_period, chi2)
	return full_time_offset, chi2
	#end SyncStreams1

def SyncStreams1_v0(streamA, sizeA, streamB, sizeB, sampling_period): #wad.
	#float SyncStreams(float streamA[], int sizeA, float streamB[], int sizeB, float sampling_period): #converted.
	#returns dt such that A(t) = B(t+dt), 
	#where A(t) and B(t) are their time functions with their own clock, so the streams both start at t=0
	#so if streamA starts at t = 0, the start of streamB is at dt
	#returns deltat, chi2
	print(f"Syncing streams: find_offset")
	offset = Find_offset(streamA, sizeA, streamB, sizeB) #int
	print(f"course offset", offset, "found, now find fine offset")
	STOF, chi2 = Find_STOF(streamA, sizeA, streamB, sizeB, offset) #STOF is sub-time-step offset fraction. will lie on (-1,1). #float
	plot_diff(streamA, sizeA, streamB, sizeB, sampling_period*(STOF + offset), sampling_period, chi2)
	return sampling_period*(STOF + offset), chi2
	#end SyncStreams1

def SyncStreams2(streamA, streamB, sampling_periodA, sampling_periodB):  #wad
	#float SyncStreams(float streamA[], int sizeA, float streamB[], int sizeB, float sampling_periodA, float sampling_periodB):
	#interpolate the stream with the lower sampling period, then run SyncStream on the interpolated streams.
	#returns deltat, chi2
	sizeA = len(streamA)
	sizeB = len(streamB)
	if sizeA ==0 or sizeB ==0:
		print(f"Error! in SyncStreams2, given bad sizes. A,B=",sizeA,sizeB)

	if abs(sampling_periodA - sampling_periodB) < 0.001*sampling_periodA: #if same period.
		return SyncStreams1(streamA, streamB, sampling_periodA)
	elif sampling_periodA > sampling_periodB: 
		#interpolate streamB since it has the finer sampling frequency
		print(f"interpolate stream B")
		Interpolated_stream, newsize = Build_interpolated_stream(streamB, sizeB, sampling_periodB, sampling_periodA) #float array
		result, chi2 = SyncStreams1(streamA, Interpolated_stream, sampling_periodA) #float
		return result, chi2
	else: 
		#interpolate streamA since it has the finer sampling frequency
		print(f"interpolate stream A")
		Interpolated_stream,newsize  = Build_interpolated_stream(streamA, sizeA, sampling_periodA, sampling_periodB) #float array
		"""
			#here, have a look at what's going on here. 
			tA = [0]
			tB = [0]
			tAi = [0]
			for i in range(1,sizeA):
				tA.append(tA[i-1]+sampling_periodA)
			for i in range(1,sizeB):
				tB.append(tB[i-1]+sampling_periodB)
			for i in range(1,newsize):
				tAi.append(tAi[i-1]+sampling_periodB)
			#print(f"StreamA time series:  ", tA[0], tA[1], tA[2])	
			print(f"StreamB time series:  ", tB[0], tB[1], tB[2])
			print(f"StreamAi time series: ", tAi[0], tAi[1], tAi[2])	

			plt.figure(3,figsize=(14,5))
			plt.plot(tA, streamA, color='darkorange',linewidth=2, label="Stream A, uninterpoltaed.")
			plt.plot(tAi, Interpolated_stream, color='blue',linewidth=2, label="Stream A, interpoltaed.")
			plt.plot(tB, streamB, color='limegreen',linewidth=2, label="Stream B")
			plt.title('Check interpolation sanity')
			plt.ylabel('Normalized Z (mm)')
			plt.xlabel('Time (s)')
			plt.legend()
			plt.show()"""

		result, chi2 = SyncStreams1(Interpolated_stream, streamB, sampling_periodB) #float
		return result, chi2
#end SyncStreams2

def test12(noffset, nB):
	#def test12(float noffset, int nB):
	#used to execute tests 1 and 2
	#returns the difference between the input offset and teh found offset
	nA = 1000 #int
	streamA = [sin(6.283*i/float(nA)) for i in range(nA)]
	streamB = [sin(6.283*(noffset + i)/float(nA)) for i in range(nB)]
	ret,chi2 = SyncStreams1(streamA, streamB, 1.,False)
	print(f"test1, input n: ",nB, "vs 1000. input offset", noffset, "returns ",ret)
	return ret - noffset
#end test12

def test12_plot(toffset, nB):
	#def test12(float noffset, int nB):
	#used to execute tests 1 and 2
	#returns the difference between the input offset and teh found offset
	#passes well.
	nA = 1000 #int
	freq = 3. #sampling frequency
	t0A = 200.3 #start of clock A in seconds
	t0B = 500.5 #start of clock B in seconds
	streamA = [sin(6.283*i/float(nA)) for i in range(nA)]
	streamB = [sin(6.283*(toffset*freq + i)/float(nA)) for i in range(nB)]
	ret,chi2 = SyncStreams1(streamA, streamB, 1./freq,False) #ah, ret is already measured in s. It's a measure of freq.
	print(f"test1, input n: ",nB, "vs 1000. input offset", toffset, "returns ",ret)

	tA = [(x/freq) + t0A for x in range(nA)]
	tB = [(x/freq) + t0B for x in range(nB)]
	absolute_time_diff = ret + tA[0] - tB[0] 
	TB = [t + absolute_time_diff for t in tB]

	plt.figure(2,figsize=(14,5))
	plt.plot(tA, streamA, color='darkorange',linewidth=2, label="Stream A")
	plt.plot(TB, streamB, color='limegreen',linewidth=2, label="Stream B")
	plt.title('Check of the Z-fitting')
	plt.ylabel('Normalized Z (mm)')
	plt.xlabel('Time (s)')
	plt.legend()
	#plt.show()
	plt.savefig('test12_plot.png')
	return ret - toffset

#test12_plot(45.6, 800)

def test3(noffset):
	#Test 3: make a pair of streams that differ by an integer number of steps. Try out Find_offset.
	nA = 1000 #int
	streamA = [sin(6.283*i/float(nA)) for i in range(nA)]
	streamB = [sin(6.283*(noffset + i)/float(nA)) for i in range(nA)]
	auto_int = Find_offset(streamA, len(streamA), streamB, len(streamB), automatic)
	hint_int = Find_offset(streamA, len(streamA), streamB, len(streamB), hint, noffset + 2)
	exhv_int = Find_offset(streamA, len(streamA), streamB, len(streamB), exhaustive)
	print(f"Test 3: with offset = ", noffset, "Automatic returns ", auto_int, "hint +2 returns ", hint_int, "and exhaustive gives ",exhv_int)

def test4(noffset, STOF):
	#takes in int noffset, float STOF
	#Test 4: make a pair of streams that differ by less than one full time step. Try out Find_STOF. 
	nA = 1000 #int
	streamA = [sin(6.283*i/float(nA)) for i in range(nA)]
	streamB = [sin(6.283*(STOF + noffset + i)/float(nA)) for i in range(nA)]
	print(f"Test4: with int offset = ",noffset, " sub sample offset inputed: ", STOF)
	result = Find_STOF(streamA, len(streamA), streamB, len(streamB), noffset) 
	print(f"Test4: inputed: ", STOF, "returns ",result)
	print(f"")
	return result

def test_set():
	#This puts the sinc code through it's paces, asking it to sinc sine waves under various conditions.
	#Here, we want the average to come out low, and for nothing to segfault or some such bad bahavior. 
	nplus = 1200
	nref = 1000
	nminu = 800
	ofs_p = 23.8
	ofs_m = -23.8
	#Test 1
	print(f"Test 1a: smoke test")
	result1 = test12(ofs_m, nplus)
	print(f"Test 1a passes. Yay, it didn't crash")
	print(f"Test 1b: do we get 0 offset when applying identical distributions?")
	res_identical= test12(0, nref)
	if res_identical == 0:
		print(f"Test 1b passes. Ientical streams produce 0 offset.")
	else: 
		print(f"Test 1b fails!. Ientical streams produce offset = ", res_identical)
		
	print(f"Test2: operation under various conditions")	
	print(f"Resulting differences from test1 suite")
	print(f"         nB =  same,   20%% more,    20%% less")
	print(f"offset  ",ofs_m,": ", test12( ofs_m, nref), result1, test12(ofs_m, nminu))
	print(f"offset   0  : ", res_identical, test12(0, nplus), test12(0, nminu))
	print(f"offset  ",ofs_p,": ", test12( ofs_p, nref), test12(ofs_p, nplus), test12(ofs_p, nminu))
	print(f"")
	print(f"Test3 suite")
	test3(-45)
	test3(0)
	test3(23)
	print(f"")
	print(f"Test4 suite")
	test4(23,0)
	test4(23,0.355)
	test4(23,-0.255)
	print(f"now things I expect to not really work")
	test4(23,0.555)
	test4(23,0.855)
	test4(23,-0.855)
	
#end test1set
#test_set()

"""
 void test(){
	#TODO: write a utility to visualize these streams.
	#Should produce a plot with two pannels: one unsynced and one synced.  

	#Test 1a: make a pair of streams with various sizes and lengths and offsets. Test that you can get chi2 without segfault from GetChi2perDOF
	#Test 1b: make sure that chi2 is 0 when streams are identical, and small when nearly identical. 

	#Test 2: make a pare of streams with various sizes and offsets (now floating point) and try out the floating point version of GetChi2perDOF

	#Test 3: make a pair of streams that differ by an integer number of steps. Try out Find_offset.

	#Test 4: make a pair of streams that differ by less than one full time step. Try out Find_STOF. passes

	#Test 5**: try out syncstream for various combinations of stream length and offset. Make sure that it's able to reproduce an offset you give it.

	#Test 6: Test the interpolator. Make a stream and up-sample it. Then display both overlaid. 
	#Test 7: try making streams with different length and different sampling frequencies. Print them. 
"""

"""def read_z_and_t_from_file(filename):
	#read in a text files with the format x\t y\t z\t t
	#converts that into streams of x,t
	#xA = [];
	#yA = [];
	zA = [];
	tA = [];
	sA = open(filename,'r')
	for line in sA.readlines():
		words = line.split()
		#xA.append(float(words[0].strip()))
		#yA.append(float(words[0].strip()))
		zA.append(float(words[0].strip()))
		tA.append(float(words[0].strip()))
	sA.close()
	return zA,tA
"""

def read_z_and_t_from_file(filename):
	#read in a text files with the format 
	#x\t y\t z\t t
	#or
	#X,2.4016438086718566e+000,Y,-3.3961773305913123e+000,Z,5.4790474559541195e+001,Time(sec),0.000000	
	#converts that into streams of x,t
	#xA = [];
	#yA = [];
	print(f"Reading from ",filename)
	zA = [];
	tA = [];
	sA = open(filename,'r')
	for line in sA.readlines():
		words = []
		for x in line.split():
       			 for y in x.split(','):
		                new.append(y)
		if words[0] == "X":
			#format 1
			if len(words) < 8:
				continue
			#xA.append(float(words[1].strip()))
			#yA.append(float(words[3].strip()))
			zA.append(float(words[5].strip()))
			tA.append(float(words[7].strip()))
		else:
			zA.append(float(words[0].strip()))
			tA.append(float(words[0].strip()))
	sA.close()
	print(f"done read")
	return zA,tA	
		
def zero_start(stream):
	#this subtracts the inital value from the stream. 
	#this expects the first few hundred samples to be about the same up to noise. 
	#first look through the first 10 samples to get a measure of the variance, 
	#then ask if sample 200 is consistent with this. if it is, subtract off the average of the first 200 samples. 
	#if it's not consistent, try with something lower.
	first = stream[0]
	return [x - first for x in stream]	
	#get the mean and variance of the first 10	
"""
	if len(stream) < 201:
		return [x - first for x in stream]	
	mean = 0
	for i in range(10):
		mean += stream[i]
	mean /= 10.
	varriance = 0
	for i in range(10):
		varriance += ( mean - stream[i])**2
	stdev = sqrt(varriance / 10. ) 
	n = 1
	if abs(stream[200] - mean) < 2.*stdev:
		n = 200
	elif abs(stream[150] - mean) < 2.*stdev:
		n = 150
	elif abs(stream[100] - mean) < 2.*stdev:
		n = 100
	elif abs(stream[50] - mean) < 2.*stdev:
		n = 50
	else: 
		n = 10
	mean = 0
	for i in range(n)
		mean += stream[i]	
	mean =/ n
	return [x - mean for x in stream]	
"""

def get_sampling_period(time_stream):
	#determins the sampling frequency of the stream in Hz
	n = len(time_stream)
	return (time_stream[n-1] - time_stream[0])/float(n)

def doSync(zA, tA, zB, tB, verbose = False):
	#this does the full suite of synchronization procedures, including normalizing the z's, and handeling the time offsets.
	#it takes in tuples of the z and t data for two data streams, A, and B, and returns the absolute time difference between them. 
	#It assumes that the beginning of the z's are zero, with very low noise, as is the case with the mu2e field mapping system measurements.
	#normalize the Z's:
	ZA = zero_start(zA)
	ZB = zero_start(zB)

	#figure out the sampling frequencies 
	TA = get_sampling_period(tA)
	TB = get_sampling_period(tB)
	
	#Do the synchronization
	if abs(TA - TB) < 0.001*TA:
		if verbose:
			print(f"Running sync for identical frequencies. Sampling period:", TA ," [s] :: ",1./TA," [Hz]")
		result,chi2 = SyncStreams1(ZA, ZB, TA ) #WAD
	else:
		if verbose:
			print(f"Running sync for different frequencies",TA, TB)
			print(f"    Sampling period of stream A: ", TA, " [s] :: ",1./TA," [Hz]")
			print(f"    Sampling period of stream B: ", TB, " [s] :: ",1./TB," [Hz]")
		result,chi2 = SyncStreams2(ZA, ZB, TA , TB ) #wad 

	absolute_time_diff = result + tA[0] - tB[0] 

	#report results. 
	"""if verbose:
		print "done with the main work of syncing streams. Minimized chi2: ", chi2
		print "Time offset between the start of the two streams", result
		print "Time offset between the two clocks", result + tA[0] - tB[0]
		#now build a time shifted version of B for the sake of plotting 
		TB = [t + absolute_time_diff for t in tB]
		plt.figure(2,figsize=(14,5))
		plt.plot(tA, ZA, color='darkorange',linewidth=2, label="Stream A")
		plt.plot(TB, ZB, color='limegreen',linewidth=2, label="Stream B")
		plt.title('Check of the Z-fitting')
		plt.ylabel('Normalized Z (mm)')
		plt.xlabel('Time (s)')
		plt.legend()
		plt.show()"""

	return absolute_time_diff
#end doSync

def read_and_sync(streamA_filename, streamB_filename):
	# void read_and_sync(string streamA, string streamB, float sampling_frequencyA, float sampling_frequencyB){
	#read in two text files with the format x\t y\t z\t t
	#then print and return the calculated sync time constant between them.
	
	#first read
	zA, tA = read_z_and_t_from_file(streamA_filename)
	zB, tB = read_z_and_t_from_file(streamB_filename)
	return doSync(zA, tA, zB, tB)
#end read_and_sync
