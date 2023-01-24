import obspy
import pyasdf #import *
import os
import math
import urllib
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core import UTCDateTime

def load_asdf(self,file='query.h5',min_distance=None,max_distance=None,invert_Z_only=False,invert_BB_only=False):

    ds=pyasdf.ASDFDataSet(file) #read file
    #check and add event and stations xml
    check_metadataset(ds)
    #load all to inputs:
    extract_event(self,ds)
    extract_network_coordinates(self,ds,min_distance,max_distance,invert_Z_only,invert_BB_only)
    extract_data(self,ds)



def check_metadataset(self,ds,tokenfile='token.txt'):
    from pyasdf.exceptions import WaveformNotInFileException
    events=ds.events
    if len(events)==0:
        #extract query:
        qq=''
        self.log('Header with event xml not found, downloading quake.ml from ESM.')
        for sta in ds.waveforms.list(): 
          #if not qq:
            for tag in ds.waveforms[sta].get_waveform_tags():
              if not qq:
                wf=ds.waveforms[sta][tag]
                evid=wf[0].stats.asdf.event_ids[0].id
                evid=evid.replace("event_id","eventid")
                qq='https://'+evid.split(':')[1] #to remove smi?
                #qq=qq.replace("_","") ??? 
        #message:
        if not os.path.exists(tokenfile):
                 raise Exception('A file with token is missing. Expected:'+tokenfile)
        tfh=open(tokenfile,'r')
        token=tfh.read()
        tfh.close()
        mess=urllib.parse.urlencode({'message':token}).encode()
              #urllib.request.urlretrieve(url=qq,data=mess,filename='quake.xml')
        urllib.request.urlretrieve(url=qq,data=mess,filename='quake.xml')
        ds.add_quakeml('quake.xml')
    for sta in ds.waveforms:
        #check if sta has station xml, how?
        try:
            staxml=sta.StationXML
        except AttributeError:
            noxml=True
            self.log('XML info for station not found, downloading from ESM')
            for tag in sta.get_waveform_tags():
              if noxml:
                network=sta[tag][0].stats.network
                station=sta[tag][0].stats.station
                qq='https://esm-db.eu/fdsnws/station/1/query?network={network}&station={station}&level=channel'.format(network=network,station=station)
                tfh=open(tokenfile,'r')
                token=tfh.read()
                tfh.close()
                mess=urllib.parse.urlencode({'message':token}).encode()
                urllib.request.urlretrieve(url=qq,data=mess,filename='station.xml')
                #urllib.request.urlretrieve(url=qq,filename='station.xml')
                ds.add_stationxml('station.xml')
                noxml=False
                os.remove('station.xml')
    #how to download only auxiliary data?
    #check MP presence
    if any(['mp' in tag[-2:] for tag in ds.waveform_tags]):
           self.log('selecting data with MP counterpart')
        #data contains manually processed waveforms:
           stalist=ds.waveforms.list()
           for sta in stalist:
              taglist=ds.waveforms[sta].get_waveform_tags() 
              for tag in taglist:
                 if tag[-2:]=='cv': #it is converted data
                    tt=tag[:-2]+'mp'
                    try:
                        del ds.waveforms[sta][tt]
                    except WaveformNotInFileException:
                        del ds.waveforms[sta][tag]


def extract_event(self,ds):
    events=ds.events
    ev=events[0]
    orig=ev.origins[0]
    mag=ev.magnitudes[0].mag
    if len(events)>1:
        self.log('More than 1 event present in asdf file, working only with first event'+str(ev.origins[0]))
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
    self.rupture_duration=dur

def extract_network_coordinates(self,ds,min_distance,max_distance,invert_Z_only,invert_BB_only):
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
          if sts.starttime==UTCDateTime('19700101_000000'):
            self.log('Omitting station '+sts.station+' due to invalid timestamp')
          else:
           if sts.channel[2] in ['N','E','Z']:
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
            if invert_BB_only and not sts.channel[0:2]=="HH":
                stn['weightN'] = stn['weightE'] = stn['weightZ']= 0.
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
    if len(stats)==0:
        raise Exception('no stations with valid data')
    self.stations=stats
    self.create_station_index()
    self.write_stations(filename='green/station.dat')

def extract_data(self,ds):
    import numpy as np
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
              if len(ds.auxiliary_data.list())>0:
                 str1[0].stats.auxpar=ds.auxiliary_data.Headers[head][tag].parameters
              stkey='_'.join([stats.network, stats.station, stats.location, stats.channel[0:2]])
              try:
                 if stn==self.stations_index[stkey]:
                    comp=str1[0].stats.channel[2]
                    f[comp]=str1[0]
                    if comp=='Z':
                         stn['useZ']=True
                         if stn['weightZ']==0.: stn['useZ']=False
                    elif comp=='N':
                         stn['useN']=True
                         if stn['weightN']==0.: stn['useN']=False
                    elif comp=='E':
                         stn['useE']=True
                         if stn['weightE']==0.: stn['useE']=False
                    if not stats.delta in self.data_deltas:
                         self.data_deltas.append(stats.delta)
                    try: 
                       units=stats.auxpar.units
                       if units=='cm/s^2':
                         f[comp].data=f[comp].data*1.e-2
                    except AttributeError:
                        #no auxiliary data, was downloaded probably downloaded without them
                        #if downloaded from esm, suppose its cm/s2
                        if "esm-db" in f[comp].stats.asdf.event_ids[0].resource_id:
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
           #check if exists
