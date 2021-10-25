#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from BayesISOLA.helpers import my_filter


def filter_ita_data(self):
   
   #in itaca
    for st in self.data:    
        st.detrend('linear')
        for tr in st:
            tr.data=tr.data-tr.data[0]
#        if fmin==None or fmax==None:
#            print('filter not set, data not filtered')
#        else:
        stats=st[0].stats
        stn = self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]
        fmin=stn['fmin']
        fmax=stn['fmax']
#        print(fmin,fmax)
        my_filter(st,fmin,fmax)
        #???self.d.data_are_corrected=True

def trim_ita_data(self, noise_slice=False, noise_starttime=None, noise_length=None):
	"""
	Filter ``self.data`` using function :func:`prefilter_data`.
	Decimate ``self.data`` to common sampling rate ``self.max_samprate``.
	Optionally, copy a time window for the noise analysis.
	Copy a slice to ``self.data``.
	
	:type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
	:param starttime: Specify the start time of trimmed data
	:type length: float
	:param length: Length in seconds of trimmed data.
	:type noise_slice: bool, optional
	:param noise_slice: If set to ``True``, copy a time window of the length ``lenght`` for later noise analysis. Copied noise is in ``self.noise``.
	:type noise_starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`, optional
	:param noise_starttime: Set the starttime of the noise time window. If ``None``, the time window starts in time ``starttime``-``length`` (in other words, it lies just before trimmed data time window).
	:type noise_length: :class:`~obspy.core.utcdatetime.UTCDateTime`, optional
	:param noise_length: Length of the noise time window (in seconds).
	"""

	starttime = self.d.event['t']+self.grid.shift_min+self.t_min
	length = self.t_max - self.t_min + self.grid.shift_max + 10
	endtime = starttime + length
	if noise_slice:
		if not noise_length:
			noise_length = length*4
		if not noise_starttime:
			noise_starttime = starttime - noise_length
			noise_endtime = starttime
		else:
			noise_endtime = noise_starttime + noise_length
		DECIMATE = int(round(self.max_samprate / self.samprate))
	for st in self.d.data:
	        stats = st[0].stats
	        fmax = self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]['fmax']
	        self.data.append(st.copy())

	for st in self.data:
		stats = st[0].stats
		fmin = self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]['fmin']
		fmax = self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]['fmax']
		decimate = int(round(st[0].stats.sampling_rate / self.max_samprate))
		if noise_slice:
			self.noise.append(st.slice(noise_starttime, noise_endtime))
			#print self.noise[-1][0].stats.endtime-self.noise[-1][0].stats.starttime, '<', length*1.1 # DEBUG
			if (len(self.noise[-1])!=3 or (self.noise[-1][0].stats.endtime-self.noise[-1][0].stats.starttime < length*1.1)) and self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]['use'+stats.channel[2]]:
				self.log('Noise slice too short to generate covariance matrix (station '+st[0].stats.station+'). Stopping generating noise slices.')
				noise_slice = False
				self.noise = []
			elif len(self.noise[-1]):
				my_filter(self.noise[-1], fmin/2, fmax*2)
				self.noise[-1].decimate(int(decimate*DECIMATE/2), no_filter=True) # noise has 2-times higher sampling than data
		self.prefilter_data(st)
		st.decimate(decimate, no_filter=True)
		st.trim(starttime, endtime)
			#print(stats.station,stats.starttime,stats['endtime'],endtime)
		# TODO: kontrola, jestli neorezavame mimo puvodni zaznam


def prefilter_data(self, st):
	"""
	FILTER OUT frequencies above Green's function computation high limit using obspy filter lowpass
	
	:param st: stream to be filtered
	:type st: :class:`~obspy.core.stream`
	"""
	f = self.freq / self.tl
	st.filter('lowpass',freq=f,corners=4) #note: corners=4 is default

 
def decimate2_shift(self):
	"""
	Generate ``self.data_shifts`` where are multiple copies of ``self.data`` (needed for plotting).
	Decimate ``self.data_shifts`` to sampling rate for inversion ``self.samprate``.
	Filter ``self.data_shifts`` by :func:`my_filter`.
	Generate ``self.d_shifts`` where are multiple vectors :math:`d`, each of them shifted according to ``self.SHIFT_min``, ``self.SHIFT_max``, and ``self.SHIFT_step``
		"""
	self.d_shifts = []
	self.data_shifts = []
	self.shifts = []
	starttime = self.d.event['t']# + self.t_min
	length = self.t_max-self.t_min
	endtime = starttime + length
	decimate = int(round(self.max_samprate / self.samprate))
	for SHIFT in range(self.grid.SHIFT_min, self.grid.SHIFT_max+1, self.grid.SHIFT_step):
		#data = deepcopy(self.data)
		shift = SHIFT / self.max_samprate
		self.shifts.append(shift)
		data = []
		#plot=True
		for st in self.data:
			st2 = st.copy()#slice(starttime+shift-self.elemse_start_origin, endtime+shift+1) # we add 1 s to be sure, that no index will point outside the range
#			st2.trim(starttime+shift-self.elemse_start_origin, endtime+shift+1, pad=True, fill_value=0.) # short records are not inverted, but they should by padded because of plotting
			stats = st2[0].stats
			stn = self.d.stations_index['_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])]
			if self.invert_displacement:
#				st2.detrend('linear')
				st2.integrate()
			if stn['accelerograph']:
				st2.integrate()
			st2.decimate(decimate, no_filter=True)
			fmin = stn['fmin']
			fmax = stn['fmax']
#			my_filter(st2, fmin, fmax) #ALREADY FILTERED
			st2.trim(starttime+shift, endtime+shift+1, pad=True, fill_value=0.) # we add 1 s to be sure, that no index will point outside the range
			data.append(st2)
		self.data_shifts.append(data)
		c = 0
		d_shift = np.empty((self.components*self.npts_slice, 1))
		for r in range(self.d.nr):
			for comp in range(3):
				if self.d.stations[r][{0:'useZ', 1:'useN', 2:'useE'}[comp]]: # this component has flag 'use in inversion'
					weight = self.d.stations[r][{0:'weightZ', 1:'weightN', 2:'weightE'}[comp]]
					for i in range(self.npts_slice):
						try:
							d_shift[c*self.npts_slice+i] = data[r][comp].data[i] * weight
						except:
							self.log('Index out of range while generating shifted data vectors. Waveform file probably too short.', printcopy=True)
							print('values for debugging: ', r, comp, c, self.npts_slice, i, c*self.npts_slice+i, len(d_shift), len(data[r][comp].data), SHIFT)
							raise Exception('Index out of range while generating shifted data vectors. Waveform file probably too short.')
					c += 1
		self.d_shifts.append(d_shift)

