import os
import math
import obspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core import Stream,Stats,UTCDateTime,Trace

def load_itaca(self,dir='./query/',min_distance=None,max_distance=None,invert_Z_only=False):
#from obspy.core import Stream,Trace
    print(dir)
#    dir=dir+"/query/"
    itadata=Stream()
    for filename in os.listdir(dir):
        ff=dir+filename
        if os.path.isfile(ff):
         data1=dynaconvert(ff)
         for tr1 in data1:
            itadata.append(tr1)
 
    extract_event(self,itadata)

    extract_network_coordinates(self,itadata,min_distance,max_distance,invert_Z_only)
    load_itadata(self,itadata)#Reading all data in given directory. Default dir=query

def extract_event(self,itadata):

    nitadat=len(itadata)
    k=0
    #? check if all data belong to same event?
    lat=[] ; long = [] ; dpth = [] ; date = [] ; time=[] ; magn =[]; agency =[]
    lat.append(itadata[0].stats.dyna.EVENT_LATITUDE_DEGREE)
    long.append(itadata[0].stats.dyna.EVENT_LONGITUDE_DEGREE)
    dpth.append(itadata[0].stats.dyna.EVENT_DEPTH_KM)
    date.append(itadata[0].stats.dyna.EVENT_DATE_YYYYMMDD)
    time.append(itadata[0].stats.dyna.EVENT_TIME_HHMMSS)
    if (itadata[0].stats.dyna.MAGNITUDE_W):
        magn.append(itadata[0].stats.dyna.MAGNITUDE_W)
    elif (itadata[0].stats.dyna.MAGNITUDE_L):
        magn.append(itadata[0].stats.dyna.MAGNITUDE_L)
    else:
        raise ValueError('No magnitude value found')
    agency.append(itadata[0].stats.dyna.HYPOCENTER_REFERENCE)
    for i in range(1,nitadat):
        if (itadata[i].stats.dyna.EVENT_LATITUDE_DEGREE==lat[k] and \
            itadata[i].stats.dyna.EVENT_LONGITUDE_DEGREE==long[k] and\
            itadata[i].stats.dyna.EVENT_DEPTH_KM==dpth[k] and \
            itadata[i].stats.dyna.EVENT_DATE_YYYYMMDD==date[k] and \
            itadata[i].stats.dyna.EVENT_TIME_HHMMSS==time[k]):
            pass
        else:
            k=k+1
            lat.append(itadata[i].stats.dyna.EVENT_LATITUDE_DEGREE)
            long.append(itadata[i].stats.dyna.EVENT_LONGITUDE_DEGREE)
            dpth.append(itadata[i].stats.dyna.EVENT_DEPTH_KM)
            date.append(itadata[i].stats.dyna.EVENT_DATE_YYYYMMDD)
            time.append(itadata[i].stats.dyna.EVENT_TIME_HHMMSS)
            magn.append(itadata[i].stats.dyna.MAGNITUDE_W)
            agency.append(itadata[i].stats.dyna.HYPOCENTER_REFERENCE)
    if(len(lat)>1):
        print('Data from different events, using only first event!')
    #convert to event info:
    itadat_sel=Stream()
    for i in range(0,1):
        for k in range(len(itadata)):
            if (itadata[k].stats.dyna.EVENT_LATITUDE_DEGREE==lat[i] and \
            itadata[k].stats.dyna.EVENT_LONGITUDE_DEGREE==long[i] and\
            itadata[k].stats.dyna.EVENT_DEPTH_KM==dpth[i] and \
            itadata[k].stats.dyna.EVENT_DATE_YYYYMMDD==date[i] and \
            itadata[k].stats.dyna.EVENT_TIME_HHMMSS==time[i]):
                itadat_sel.append(itadata[k])

    #overwrite data by selection:
        itadata=itadat_sel

        t=UTCDateTime(date[i]+time[i])
        self.event={'lat' : lat[i], 'lon' : long[i], 'depth' : float(dpth[i])*1e3, 'mag' : float(magn[i]),\
                    't' : t, 'agency' : agency[i]}
        self.log('\nHypocenter location:\n  Agency: {agency:s}\n  Origin time: {t:s}\n  Lat {lat:8.3f}   Lon {lon:8.3f} \
        Depth{d:4.1f} km'.format(t=t.strftime('%Y-%m-%d %H:%M:%S'), lat=float(lat[i]), lon=float(long[i]), d=float(dpth[i]),\
                                 agency=agency[i]))
        self.rupture_length = math.sqrt(111 * 10**self.event['mag'])		# M6 ~ 111 km2, M5 ~ 11 km2 		REFERENCE NEEDED
        #get event name
        self.event['name']=itadata[0].stats.dyna.EVENT_NAME
#generate soutype.dat using scaling relation for rupture duration:
 #write STF file soutype:
        m0=10**(self.event['mag']*3./2.+9.1)
        dur=10.**(0.31*math.log10(m0)-4.9) #Courboulex 2016 SRL not-sub
        self.set_source_time_function('triangle',dur)
#this was replaced by new set_source_time_function
#        gwd='green/'
#       if not os.path.exists(gwd):
#            os.mkdir(gwd)
#        sfile=gwd+'soutype.dat'
#        sf=open(sfile,'w')
#        m0=10**(isola.event['mag']*3./2.+9.1)
#        dur=10.**(0.31*log10(m0)-4.9) #Courboulex 2016 SRL not-sub
#        sf.write('4\n%5.1f\n1.\n1\n' %dur) #triangular (4) stf of duration dur in velocity (1), 3rd line is not used...
#        sf.close()


def extract_network_coordinates(self,itadata,min_distance,max_distance,invert_Z_only):
#def extract_network_coordinates(itadata):
    nitadat=len(itadata)
    lat=[] ; long = [] ; stnid = [] ; netw = [] ; chann=[]; loc=[]
#    dist1=[];az1=[]
    lat.append(itadata[0].stats.dyna.STATION_LATITUDE_DEGREE)
    long.append(itadata[0].stats.dyna.STATION_LONGITUDE_DEGREE)
    stnid.append(itadata[0].stats.station)
    netw.append(itadata[0].stats.network)
    chann.append(itadata[0].stats.channel)
    loc.append(itadata[0].stats.location)
    kmax=1
#    dist1.append(itadata[0].stats.dyna.EPICENTRAL_DISTANCE_KM)
#    az1.append(itadata[0].stats.dyna.EARTHQUAKE_BACKAZIMUTH_DEGREE)
    for i in range(1,nitadat):
        skip=False
        for k in range(0,kmax):
            if (itadata[i].stats.dyna.STATION_LATITUDE_DEGREE==lat[k] and \
                itadata[i].stats.dyna.STATION_LONGITUDE_DEGREE==long[k] and\
                itadata[i].stats.station==stnid[k] and \
                itadata[i].stats.location==loc[k] and \
                itadata[i].stats.network==netw[k]):
                skip=True
        if itadata[i].stats.starttime== toUTCDateTime('19700101_000000'):
                print('Omitting station ',itadata[i].stats.station,' due to invalid timestamp')
                skip=True
        if not skip: 
            kmax=kmax+1
            lat.append(itadata[i].stats.dyna.STATION_LATITUDE_DEGREE)
            long.append(itadata[i].stats.dyna.STATION_LONGITUDE_DEGREE)
            stnid.append(itadata[i].stats.station)
            netw.append(itadata[i].stats.network)
            chann.append(itadata[i].stats.channel)
            loc.append(itadata[i].stats.location)
#            dist1.append(itadata[i].stats.dyna.EPICENTRAL_DISTANCE_KM)
#            az1.append(itadata[i].stats.dyna.EARTHQUAKE_BACKAZIMUTH_DEGREE)
        
    #cmajz od Jirky: tuto nam uz bude treba isolu....
    if min_distance==None:
        min_distance = 2*self.rupture_length
    if max_distance==None:
        max_distance = 1000 * 2**(self.event['mag']*2.)

    stats=[]
 
    for i in range(0,len(lat)):
        stn={'code':stnid[i],'lat':lat[i],'lon':long[i],'network':netw[i],\
             'location':loc[i],'channelcode':chann[i][0:2],'model':''}
        if obspy.__version__[0] == '0':
            g = Geod(ellps='WGS84')
            az,baz,dist = g.inv(self.event['lon'], self.event['lat'], long[i], lat[i])
        else:
            dist,az,baz = gps2dist_azimuth(self.event['lat'], self.event['lon'], lat[i], long[i])
        model=''
        if model not in self.models:
            self.models[model]=0
        stn['az'] = az; stn['dist'] = dist
#        stn['az'] = az1[i]; stn['dist'] = dist1[i]
        stn['useN'] = stn['useE'] = stn['useZ'] = False
        stn['accelerograph'] = False
        if invert_Z_only:
           stn['weightN'] = stn['weightE'] = 0.
        else:
           stn['weightN'] = stn['weightE'] = 1.
        stn['weightZ'] = 1.
        if dist > min_distance and dist < max_distance:   
#        if dist1[i] > 0. and dist1[i] < 1000.:
            stats.append(stn)

    stats = sorted(stats, key=lambda stn: stn['dist']) # seradit podle vzdalenosti        
#    return stats
## number of stations:
    if len(stats) > 1000:
        print('Warning: Number of stations > 1000. Using only first 1000 stations. If you need all of them, change nrp in param.inc and elemse.for in green function calculation and then return here and correct!')
        stats = stats[0:999] # BECAUSE OF GREENS FUNCTIONS CALCULATION

    self.stations=stats
#these are external
    self.create_station_index()
    self.write_stations(filename='green/station.dat')

def load_itadata(self,itadata):
#def load_itadata(itada):
    nitadat=len(itadata)
    stns=self.stations#extract_network_coordinates(itadata)
    nstats=len(stns)

    for i in range(0,nstats):
        sta=stns[i]['code']
        lat=stns[i]['lat']
        long=stns[i]['lon']
#    ch=stations[i]['channelcode']
        f={}
        for k in range(0,nitadat):
            if (itadata[k].stats.station==sta and \
                itadata[k].stats.dyna.STATION_LATITUDE_DEGREE==lat and \
                itadata[k].stats.dyna.STATION_LONGITUDE_DEGREE==long):
                comp=itadata[k].stats.channel[2]
                f[comp]=itadata[k]
                if comp=='Z':
                    self.stations[i]['useZ']=True
                elif comp=='N':
                    self.stations[i]['useN']=True
                elif comp=='E':
                    self.stations[i]['useE']=True
                if itadata[k].stats.dyna.UNITS=='cm/s^2':
                    f[comp].data=f[comp].data*1.e-2
                if itadata[k].stats.dyna.DATA_TYPE=='ACCELERATION':
                    self.stations[i]['accelerograph']=True
                #append deltas 
                if not itadata[k].stats.delta in self.data_deltas:
                    self.data_deltas.append(itadata[k].stats.delta)
#        print(i,isola.stations[i]['useZ'], isola.stations[i]['useN'],isola.stations[i]['useE'])
        #what if one trace does not exist?
        st=Stream()
        if (self.stations[i]['useZ']):
            st.append(f['Z'])
        if (self.stations[i]['useN']):
            st.append(f['N'])
        if (self.stations[i]['useE']):
            st.append(f['E'])
        # 
#        st=Stream(traces=[f['Z'],f['N'],f['E']])
#        st.plot()
#        print(st[0].stats.station,st[0].stats.starttime,st[0].stats.endtime)
        self.data.append(st)

def toUTCDateTime(value):
    #from obspy.core import UTCDateTime
    try:
        date, time = value.split('_')
    except ValueError:
        date = value
        
    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])
    
    hour = int(time[0:2])
    mins = int(time[2:4])
    secs = float(time[4:])
    
    return UTCDateTime(year, month, day, hour, mins) + secs

def strtofloat(sf):
    try:
        x = float(sf)
    except:
        return None
    return x

def strtoint(sf):
    try:
        x = int(sf)
    except:
        return None
    return x


def dynaconvert(filename_in):
    """
LICENSE:

dyna-convert.py is released under the following BSD-style license:

Copyright (c) 2014, INGV Milano-Pavia
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
    """
    import numpy as np
    
    headers = {}
#    data = StringIO()
        
    fh = open(filename_in, 'rt')
    for i in range(64): 
        key, value = fh.readline().strip().split(':', 1)
        headers[key.strip()] = value.strip()
    # create ObsPy stream object
    stream = Stream()
    header = Stats()
    header['dyna'] = {}

    header['network'] = headers['NETWORK']
    header['station'] = headers['STATION_CODE']
    header['location'] = headers['LOCATION'] 
    header['channel'] = headers['STREAM']
    try:
        header['starttime'] = toUTCDateTime(headers['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS']) # use toUTCDateTime to convert from DYNA format
    except:
        header['starttime'] = toUTCDateTime('19700101_000000')
    header['sampling_rate'] = 1/float(headers['SAMPLING_INTERVAL_S'])
    header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
    header['npts'] = int(headers['NDATA'])
    header['calib'] = 1 # not in file header

##DYNA dict float data
    header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(headers['EVENT_LATITUDE_DEGREE'])
    header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(headers['EVENT_LONGITUDE_DEGREE'])
    header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
    header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
    header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
    header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
    header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(headers['STATION_LATITUDE_DEGREE'])
    header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(headers['STATION_LONGITUDE_DEGREE'])
    header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
    header['dyna']['SITE_CLASSIFICATION_EC8']=headers['SITE_CLASSIFICATION_EC8']
    header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(headers['EPICENTRAL_DISTANCE_KM'])
    header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
    header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
    header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(headers['INSTRUMENTAL_FREQUENCY_HZ'])
    header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(headers['INSTRUMENTAL_DAMPING'])
    header['dyna']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

# data type is acceleration
    if headers['DATA_TYPE'] == "ACCELERATION" \
    or headers['DATA_TYPE'] == "ACCELERATION RESPONSE SPECTRUM":
        header['dyna']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
        header['dyna']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
# data type is velocity
    if headers['DATA_TYPE'] == "VELOCITY" \
    or headers['DATA_TYPE'] == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
        header['dyna']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
        header['dyna']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
# data type is displacement
    if headers['DATA_TYPE'] == "DISPLACEMENT" \
    or headers['DATA_TYPE'] == "DISPLACEMENT RESPONSE SPECTRUM":
        header['dyna']['PGD_CM'] = strtofloat(headers['PGD_CM'])
        header['dyna']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])
        
    header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(headers['LOW_CUT_FREQUENCY_HZ'])
    header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(headers['HIGH_CUT_FREQUENCY_HZ'])

##DYNA dict int data
    header['dyna']['STATION_ELEVATION_M'] = strtoint(headers['STATION_ELEVATION_M'])
    header['dyna']['SENSOR_DEPTH_M'] = strtoint(headers['SENSOR_DEPTH_M'])
    header['dyna']['N_BIT_DIGITAL_CONVERTER'] =  strtoint(headers['N_BIT_DIGITAL_CONVERTER'])
    header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

##DYNA dict string data
    header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
    header['dyna']['EVENT_ID'] = headers['EVENT_ID']
    header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
    header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
    header['dyna']['MAGNITUDE_W_REFERENCE'] = headers['MAGNITUDE_W_REFERENCE']
    header['dyna']['MAGNITUDE_L_REFERENCE'] = headers['MAGNITUDE_L_REFERENCE']
    header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
    header['dyna']['STATION_NAME'] = headers['STATION_NAME']
    header['dyna']['SITE_CLASSIFICATION_EC8'] = headers['SITE_CLASSIFICATION_EC8']
    header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers['MORPHOLOGIC_CLASSIFICATION']
    header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers['DATE_TIME_FIRST_SAMPLE_PRECISION']
    header['dyna']['UNITS'] = headers['UNITS']
    header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
    header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers['INSTRUMENT_ANALOG/DIGITAL']
    header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
    header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
    header['dyna']['LATE_NORMAL_TRIGGERED'] = headers['LATE/NORMAL_TRIGGERED']
    header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
    header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
    header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
    header['dyna']['PROCESSING'] = headers['PROCESSING']
    header['dyna']['DATA_LICENSE'] = headers['DATA_LICENSE']
    header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers['DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
    header['dyna']['DATA_CITATION'] = headers['DATA_CITATION']
    header['dyna']['DATA_CREATOR'] = headers['DATA_CREATOR']
    header['dyna']['ORIGINAL_DATA_MEDIATOR_CITATION'] = headers['ORIGINAL_DATA_MEDIATOR_CITATION']
    header['dyna']['ORIGINAL_DATA_MEDIATOR'] = headers['ORIGINAL_DATA_MEDIATOR']
    header['dyna']['ORIGINAL_DATA_CREATOR_CITATION'] = headers['ORIGINAL_DATA_CREATOR_CITATION']
    header['dyna']['ORIGINAL_DATA_CREATOR'] = headers['ORIGINAL_DATA_CREATOR']
    header['dyna']['USER1'] = headers['USER1']
    header['dyna']['USER2'] = headers['USER2']
    header['dyna']['USER3'] = headers['USER3']
    header['dyna']['USER4'] = headers['USER4']
    header['dyna']['USER5'] = headers['USER5']

# read data
    data = np.loadtxt(fh, dtype='float32')
    if headers['DATA_TYPE'][-8:] == "SPECTRUM":
        data_1 = np.array([], dtype=np.float32)
        data_2 = np.array([], dtype=np.float32)
        for j in range(len(data)):
            for i in range(2):
                if i == 0:
                    data_1 = np.append(data_1,data[j][i])
                elif i == 1:
                    data_2 = np.append(data_2,data[j][i])
        stream.append(Trace(data=data_1, header=header))
        stream.append(Trace(data=data_2, header=header))
    else:
        stream.append(Trace(data=data, header=header))

    fh.close()
#    print(header)
    return stream


