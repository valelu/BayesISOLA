#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os.path
import matplotlib.pyplot as plt

from BayesISOLA.helpers import next_power_of_2

def suspect_mouse(self,fcrit1=0.05,fmax=20.,ampcrit1=0.1,fcrit2=0.05,ampcrit2=0.4,figures=None, multichannel=True):
        """
        Simple routine to substitute mousetrap when raw data are unavailable:
           compares the maximal amplitude of signal spectrum below and above critical frequency fcrit
           if maximal amplitude below fcrit is larger than above fcrit -> mouse is suspected
        Note: if the mouse signal is week (its spectrum is smaller), it will not work...But hopefully the mouse signal will not be fitted
        :param fcrit: critical frequency to dissect the spectrum, where the mouse signal is present
	"""
#some parts copied from detect_mouse...
        self.log('\nSimple spectral mouse suspicion:')
        if multichannel : self.log('Assuming different parameter for HH and HN/HL spectral mouse check.')
        out = ''
        for st0 in self.data:
                st = st0.copy()
		#demean(st)
		# calculate spectra and their abs 
                sta = st[0].stats.station
                for comp in range(len(st)):
                    tr=st[comp]
                    stats = tr.stats
                    stkey='_'.join([stats.network, sta, stats.location, stats.channel[0:2]])
                    #check mouse only if use component is true:
                    if self.stations_index[stkey]['use'+stats.channel[2]]==True:
                        nfft=next_power_of_2(tr.stats.npts)
                        nft2=int(nfft/2)
                        df=1./(nfft*stats.delta)
                        if multichannel:
                          if stats.channel[0:2]=="HH": #different parameters for broadband data and also integrate spectra to velocities!
                            icrit=int(fcrit2/df)
                            ampcrit=ampcrit2
                            sp=np.fft.fft(np.cumsum(tr.data)*stats.delta,nfft)
                          else:
			#try:
                        #except:
			#	out += '  ' + sta + ' ' + stats.channel + ': MOUSE detecting problem (record too short?), ignoring component in inversion\n'
		        #	self.stations_index['_'.join([stats.network, sta, stats.location, stats.channel[0:2]])]['use'+stats.channel[2]] = False
			#else:
                            sp=np.fft.fft(tr.data,nfft)
                            icrit=int(fcrit1/df)
                            ampcrit=ampcrit1



                        imax=int(fmax/df)
                        sp=abs(sp[:nft2])
#osetri chybu ak imax<icrit
#....
                        amp1=np.max(sp[:icrit])
                        amp2=np.max(sp[icrit:imax])
				#onset, amp, dummy, dummy, fit = m1.params(degrees=True)
				#amp = abs(amp)
                        detected=False
				#if (amp > 50e-8 and fit > 0.6) or (amp > 10e-8 and fit > 0.8) or (amp > 7e-8 and fit > 0.9) or (amp > 5e-9 and fit > 0.94) or (fit > 0.985): # DEBUGGING: fit > 0.95 in the before-last parentheses?
                        if amp1/amp2>ampcrit:
                                out += '  ' + sta + ' ' + stats.channel + ': MOUSE suspected, ignoring component in inversion \n'#(time of onset: {o:6.1f} s, amplitude: {a:10.2e} m s^-2, fit: {f:7.2f})\n'.format(o=onset-t_start_origin, a=amp, f=fit)
                                self.stations_index['_'.join([stats.network, sta, stats.location, stats.channel[0:2]])]['use'+stats.channel[2]] = False
                                detected=True
                                if figures:
                                        if figures is True:
                                            figures = os.path.join(self.outdir, 'mouse')
                                        if not os.path.exists(figures): #and figures_mkdir:
                                            os.mkdir(figures)
					#m1.plot(st[comp], outfile=os.path.join(figures, 'mouse_'+('no','YES')[detected]+'_'+sta+str(comp)+'.png'), xmin=t_start_origin-60, xmax=t_start_origin+240, ylabel='raw displacement [counts]', title="{{net:s}}:{{sta:s}} {{ch:s}}, fit: {fit:4.2f}".format(fit=fit))
                                        fr=[i*df for i in range(nft2)]
                                        colors=['r','g','b']
                                        plt.xscale('log')
                                        plt.yscale('log')
                                        plt.xlabel('frequency [Hz]')
                                        plt.ylabel('Acceleration amplitude spectrum')
                                        #if (np.amin(sp)<1.e-10):
                                        #   plt.ylim(bottom=1.e-8,top=5*np.amax(sp))
                                        p=plt.plot(fr,sp,colors[comp])
                                        plt.savefig(os.path.join(figures,'mouse_'+('no','YES')[detected]+'_'+sta+stats.channel+'.png'))
                                        plt.clf()
        self.logtext['mouse'] = out
        self.log(out, newline=False)


