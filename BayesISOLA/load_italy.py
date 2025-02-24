#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path

import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings('ignore', '.*Conversion of the second argument of issubdtype from.*')
	# Avoids following warning:
		# /usr/lib/python3.7/site-packages/obspy/signal/detrend.py:31: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.


class load_italy:
	"""
     Load all necessary data for MT inversion.

    :type logfile: string, optional
    :param logfile: path to the logfile (default '$outdir/log.txt')
    :type outdir: string, optional
    :param outdir: a directory, where the outputs are saved (default 'output')
    :type output_mkdir: bool, optional
    :param output_mkdir: if ``True`` (default), creates a dir for output in case it does not exists
	
    .. rubric:: _`Variables`
    
    ``data`` : list of :class:`~obspy.core.stream`
        Data for the inversion. They are loaded by 
         load_itaca
        or load_asdf
       ??? possibility to access remotely using curl?
        The list is ordered ascending by epicentral distance of the station. After processing, the list references to the same streams as ``data``.
    ``data_deltas`` : list of floats
        List of ``stats.delta`` from ``data_raw`` and ``data``.
    ``data_are_corrected`` : bool
        Flag whether the instrument response correction was performed.
    ``event`` : dictionary
        Information about event location and magnitude.
    ``stations`` : list of dictionaries
        Information about stations used in inversion. The list is ordered ascending by epicentral distance.
    ``stations_index`` : dictionary referencing ``station`` items.
        You can also access station information like ``self.stations_index['NETWORK_STATION_LOCATION_CHANNELCODE']['dist']`` insted of ``self.stations[5]['dist']``. `CHANNELCODE` are the first two letters of channel, e.g. `HH`.
    ``nr`` : integer
        Number of stations used, i.e. length of ``stations`` and ``data``.
    ``logtext`` : dictionary
        Status messages for :func:`html_log`
    ``models`` : dictionary
        Crust models used for calculating synthetic seismograms
    ``rupture_length`` : float
        Estimated length of the rupture in meters
    ``stf_description`` : string
        Text description of source time function
	"""


	from BayesISOLA._input_crust import read_crust
	from BayesISOLA._input_event import set_source_time_function
	from BayesISOLA._input_network import create_station_index, write_stations
	from BayesISOLA._input_itaca import load_itaca#, extract_event, extract_network_coordinates, load_itadata
	from BayesISOLA._input_asdf import load_asdf	
	
	def __init__(self, outdir='output', logfile='$outdir/log.txt', output_mkdir=True):
		self.outdir = outdir
		if not os.path.exists(outdir) and output_mkdir:
			os.mkdir(outdir)
		self.logfile = open(logfile.replace('$outdir', self.outdir), 'w', 1)
		self.data = []
		self.data_deltas = [] # list of ``stats.delta`` values of traces in ``self.data`` or ``self.data_raw``
		self.logtext = {}
		self.models = {}

	def __exit__(self, exc_type, exc_value, traceback):
		self.__del__()
		
	def __del__(self):
		self.logfile.close()
		del self.data

	def log(self, s, newline=True, printcopy=False):
		"""
		Write text into log file
		
		:param s: Text to write into log
		:type s: string
		:param newline: if is ``True``, add LF symbol (\\\\n) at the end
		:type newline: bool, optional
		:param printcopy: if is ``True`` prints copy of ``s`` also to stdout
		:type printcopy: bool, optional
		"""
		self.logfile.write(s)
		if newline:
			self.logfile.write('\n')
		if printcopy:
			print(s)


