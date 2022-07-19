#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth	
from BayesISOLA.fileformats import read_elemse
from BayesISOLA.helpers import my_filter

def VR_of_components(self, n=1):
	"""
	Calculates the variance reduction from each component and the variance reduction from a subset of stations.
	
	:param n: minimal number of components used
	:type n: integer, optional
	:return: maximal variance reduction from a subset of stations
	
	Add the variance reduction of each component to ``self.stations`` with keys ``VR_Z``, ``VR_N``, and ``VR_Z``.
	Calculate the variance reduction from a subset of the closest stations (with minimal ``n`` components used) leading to the highest variance reduction and save it to ``self.max_VR``.
	"""
	npts = self.d.npts_slice
	data = self.d.data_shifts[self.centroid['shift_idx']]
	elemse = read_elemse(self.inp.nr, self.d.npts_elemse, 'green/elemse'+self.centroid['id']+'.dat', self.inp.stations, self.d.invert_displacement) # read elemse
	for r in range(self.inp.nr):
		for e in range(6):
			my_filter(elemse[r][e], self.inp.stations[r]['fmin'], self.inp.stations[r]['fmax'])
			elemse[r][e].trim(UTCDateTime(0)+self.d.elemse_start_origin)
	MISFIT = 0
	NORM_D = 0
	COMPS_USED = 0
	max_VR = -99
	self.VRcomp = {}
	for sta in range(self.inp.nr):
		SYNT = {}
		for comp in range(3):
			SYNT[comp] = np.zeros(npts)
			for e in range(6):
				SYNT[comp] += elemse[sta][e][comp].data[0:npts] * self.centroid['a'][e,0]
		comps_used = 0
		for comp in range(3):
			if self.cova.Cd_inv and not self.inp.stations[sta][{0:'useZ', 1:'useN', 2:'useE'}[comp]]:
				self.inp.stations[sta][{0:'VR_Z', 1:'VR_N', 2:'VR_E'}[comp]] = None
				continue
			synt = SYNT[comp]
			try:
				d = data[sta][comp][0:npts]
			except IndexError:
				continue	
			if self.cova.LT3:
				d    = np.zeros(npts)
				synt = np.zeros(npts)
				x1 = -npts
				for COMP in range(3):
					if not self.inp.stations[sta][{0:'useZ', 1:'useN', 2:'useE'}[COMP]]:
						continue
					x1 += npts; x2 = x1+npts
					y1 = comps_used*npts; y2 = y1+npts
					d    += np.dot(self.cova.LT3[sta][y1:y2, x1:x2], data[sta][COMP].data[0:npts])
					synt += np.dot(self.cova.LT3[sta][y1:y2, x1:x2], SYNT[COMP])
				
			elif self.cova.Cd_inv:
				d    = np.dot(self.cova.LT[sta][comp], d)
				synt = np.dot(self.cova.LT[sta][comp], synt)
			elif self.cova.Cd_inv_shifts:
				d=np.dot(self.cova.LT_shifts[self.centroid['shift_idx']][sta][comp],d)
				synt=np.dot(self.cova.LT_shifts[self.centroid['shift_idx']][sta][comp],synt)
			else:
				pass
			comps_used += 1
			misfit = np.sum(np.square(d - synt))
			norm_d = np.sum(np.square(d))
			VR = 1 - misfit / norm_d
			self.inp.stations[sta][{0:'VR_Z', 1:'VR_N', 2:'VR_E'}[comp]] = VR
			if self.inp.stations[sta][{0:'useZ', 1:'useN', 2:'useE'}[comp]]:
				MISFIT += misfit
				NORM_D += norm_d
				VR_sum = 1 - MISFIT / NORM_D
				COMPS_USED += 1
				#print sta, comp, VR, VR_sum # DEBUG

		if COMPS_USED >= n:
			if COMPS_USED > 1:
				self.VRcomp[COMPS_USED] = VR_sum
			if VR_sum >= max_VR:
				max_VR = VR_sum
				self.max_VR = (VR_sum, COMPS_USED)
	write_VR_by_comp(self)
	return max_VR

def write_VR_by_comp(self):
    #calculate VR without matrix weighing:
    npts = self.d.npts_slice
    data = self.d.data_shifts[self.centroid['shift_idx']]
    elemse = read_elemse(self.inp.nr, self.d.npts_elemse, 'green/elemse'+self.centroid['id']+'.dat', self.inp.stations, self.d.invert_displacement) # read elemse
    for r in range(self.inp.nr):
        for e in range(6):
            my_filter(elemse[r][e], self.inp.stations[r]['fmin'], self.inp.stations[r]['fmax'])
            elemse[r][e].trim(UTCDateTime(0)+self.d.elemse_start_origin)
    MISFIT = 0
    NORM_D = 0
    COMPS_USED = 0
    max_VR = -99
    self.VRcomp = {}
    for sta in range(self.inp.nr):
        SYNT = {}
        for comp in range(3):
            SYNT[comp] = np.zeros(npts)
            for e in range(6):
                SYNT[comp] += elemse[sta][e][comp].data[0:npts] * self.centroid['a'][e,0]
        comps_used = 0
        for comp in range(3):
            if not self.inp.stations[sta][{0:'useZ', 1:'useN', 2:'useE'}[comp]]:
                self.inp.stations[sta][{0:'VR_Z', 1:'VR_N', 2:'VR_E'}[comp]] = None
                continue
            synt = SYNT[comp]
            try:
                    d = data[sta][comp][0:npts]
                    comps_used += 1
                    misfit = np.sum(np.square(d - synt))
                    norm_d = np.sum(np.square(d))
                    VR = 1 - misfit / norm_d
            except IndexError:
                    pass	
            self.inp.stations[sta][{0:'VR_Z', 1:'VR_N', 2:'VR_E'}[comp]] = VR

     #recalculate distance
    for r in range(self.inp.nr):
        if obspy.__version__[0] == '0':
               az,baz,dist=g.inv(self.centroid['lon'],self.centroid['lat'],self.inp.stations[r]['lon'],self.inp.stations[r]['lat'])
        else:
               dist,az,baz=gps2dist_azimuth(self.centroid['lat'],self.centroid['lon'],self.inp.stations[r]['lat'],self.inp.stations[r]['lon']) 
        self.inp.stations[r]['dist2']=dist
     #write output
     #station code, distance, VR_Z,VR_N,VR_E
    wrf=open(self.inp.outdir+'/VR_by_comp.dat','w')
    wrf.write('Station \t distance \t VR_Z \t VR_N \t VR_E \n')
    for r in range(self.inp.nr):
      if self.inp.stations[r]['VR_N'] and self.inp.stations[r]['VR_E'] and self.inp.stations[r]['VR_Z']:
        wrf.write('{0:s}:{1:s}  {2:10.1f}  {3:5.0f}%  {4:5.0f}%  {5:5.0f}% \n'.format(self.inp.stations[r]['network'],self.inp.stations[r]['code'],self.inp.stations[r]['dist2']/1000.,self.inp.stations[r]['VR_Z']*100,self.inp.stations[r]['VR_N']*100,self.inp.stations[r]['VR_E']*100))
    wrf.close()
