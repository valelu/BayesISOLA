import obspy
import pyasdf #import *
import os
import math
from obspy.geodetics.base import gps2dist_azimuth

def load_asdf(self,file='query.h5',min_distance=None,max_distance=None,invert_Z_only=False):

    ds=pyasdf.ASDFDataSet(file) #read file
    #load all to inputs:
    extract_event(self,ds)
    extract_network_coordinates(self,ds,min_distance,max_distance,invert_Z_only)
    extract_data(self,ds)


def extract_event(self,ds):
    events=ds.events
    ev=events[0]
    orig=ev.origins[0]
    mag=ev.magnitudes[0].mag
    if len(events)>1:
        self.log('More than 1 event present in asdf file, working only with first event',ev.origins[0])
    #create event
    self.event={'lat' : orig.latitude, 'lon' : orig.longitude, 'depth' : orig.depth*1e3, 'mag' : mag ,\
                    't' : orig.time, 'agency' : orig.creation_info.author}
    self.log('\nHypocenter location:\n  Agency: {agency:s}\n  Origin time: {t:s}\n  Lat {lat:8.3f}   Lon {lon:8.3f} \
        Depth{d:4.1f} km'.format(t=orig.time.strftime('%Y-%m-%d %H:%M:%S'), lat=float(orig.latitude), lon=float(orig.longitude), d=float(orig.depth),\
                                 agency=orig.creation_info.author))
    self.rupture_length = math.sqrt(111 * 10**self.event['mag'])
    #get event name:
    for desc in ev.event_descriptions:
        if desc.type=='earthquake name':
            evname=desc.text
    self.event['name']=evname
    m0=10**(self.event['mag']*3./2.+9.1)
    dur=10.**(0.31*math.log10(m0)-4.9) #scaling relation: Courboulex 2016 SRL not-sub
    self.set_source_time_function('triangle',dur)

def extract_network_coordinates(self,ds,min_distance,max_distance,invert_Z_only):
#get station coordinates
    
    nstats=len(ds.waveforms)

    if min_distance==None:
        min_distance = 2*self.rupture_length
    if max_distance==None:
        max_distance = 1000 * 2**(self.event['mag']*2.)

    stats=[]

    for station in ds.waveforms:
        lat=station.coordinates['latitude']
        long=station.coordinates['longitude']
        for tag in station.get_waveform_tags():
            str1=station[tag]
            sts=str1[0].stats
            stn={'code': sts.station,'lat': lat,'lon':long,'network': sts.network,\
             'location':sts.location,'channelcode':sts.channel[0:2],'model':''}
            if obspy.__version__[0] == '0':
               g = Geod(ellps='WGS84')
               az,baz,dist = g.inv(self.event['lon'], self.event['lat'], long, lat)
            else:
               dist,az,baz = gps2dist_azimuth(self.event['lat'], self.event['lon'], lat, long)
            model=''
            if model not in self.models:
                self.models[model]=0
            stn['az'] = az; stn['dist'] = dist
#        stn['az'] = az1[i]; stn['dist'] = dist1[i]
            stn['useN'] = stn['useE'] = stn['useZ'] = False
#            stn['accelerograph'] = False
            stn['accelerograph'] = "acc" in tag
            if invert_Z_only:
                stn['weightN'] = stn['weightE'] = 0.
            else:
                stn['weightN'] = stn['weightE'] = 1.
            stn['weightZ'] = 1.
            if dist > min_distance and dist < max_distance:   
#        if dist1[i] > 0. and dist1[i] < 1000.:
               if not stn in stats:
                   stats.append(stn)

    stats = sorted(stats, key=lambda stn: stn['dist']) # seradit podle vzdalenosti        
#    return stats
## number of stations:
    if len(stats) > 1000:
        print('Warning: Number of stations > 1000. Using only first 1000 stations. If you need all of them, change nrp in param.inc and elemse.for in green function calculation, recompile and then return here and correct!')
        stats = stats[0:999] # BECAUSE OF GREENS FUNCTIONS CALCULATION

    self.stations=stats
    self.create_station_index()
    self.write_stations(filename='green/station.dat')

def extract_data(self,ds):
#prepare data
    for stn in self.stations:
       f={}
       stkey=None
       head='_'.join([stn['network'],stn['code']])
       station=ds.waveforms[head]
       for tag in station.get_waveform_tags():
              str1=station[tag]
              stats=str1[0].stats
              #add auxiliary data as stream some stream flag
              #header #tag from tag
#              head='_'.join([stats.network,stats.station])
              str1[0].stats.auxpar=ds.auxiliary_data.Headers[head][tag].parameters
              stkey='_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])
              try:
                 if stn==self.stations_index[stkey]:
                    comp=str1[0].stats.channel[2]
                    f[comp]=str1[0]
                    if comp=='Z':
                         stn['useZ']=True
                    elif comp=='N':
                         stn['useN']=True
                    elif comp=='E':
                         stn['useE']=True
                    if not stats.delta in self.data_deltas:
                         self.data_deltas.append(stats.delta)
                    if stats.auxpar.units=='cm/s^2':
                         f[comp].data=f[comp].data*1.e-2
              except KeyError:
                 #stkey not in station index probably due to distance constraints
                 continue
       st=obspy.Stream()
       if f:
              if stn['useZ']:
                  st.append(f['Z'])
              if stn['useN']:
                  st.append(f['N'])
              if stn['useE']:
                  st.append(f['E'])
              self.data.append(st)

    
