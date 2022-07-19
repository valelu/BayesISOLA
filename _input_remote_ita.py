from urllib.request import urlopen,urlretrieve
from urllib import parse
from obspy.core import Stream
import io
import os

def load_from_esm(self,eventID,tokenfile='token.txt',min_distance=None,max_distance=None,invert_Z_only=False,invert_BB_only=False,user_email='',user_password='',proc_type='CV,MP',data_type='ACC',form='hdf5',save_to=None):
    """
    Downloads waveforms from esm_db and prepares data
    
    :param eventID: ID of event to download from ESM-db. See also: https://esm-db.eu/esmws/eventdata/1/query-options.html.
    :type eventID: string
    :param tokenfile: path to file with esm token. Default: token.txt.
    :type tokenfile: string or other path-like object, optional
    :param min_distance: minimal distance criterion of stations. Default: estimated from magnitude scaling relations.
    :type min_distance: float,optional
    :param max_distance: maximal distance criterion for stations. Default: none, but carefult with the fmax!
    :type max_distance: float, optional
    :param invert_Z_only: If true, omit horizontal components in the inversion. Default: False.
    :type invert_Z_only: boolean, optional
    :param user_email: User email to log in to ESM. This will be used if tokenfile is not found.
    :type user_email: string,optional
    :param user_password: Password to log in to ESM with user_email. Will be used instead of tokenfile.
    :type user_password: string, optional
    :param proc_type: processing type of data to be downloaded. MUST CONTAIN CV!!! because only converted data is used for the inversion. If MP present, will check presence of MP against CV as "bad-quality" flag. Default: 'CV,MP'.
    :type proc_type: string, optional
    :param data_type: Data type for data to download. For ESM data, only accelerations are available as not-processed. Default: 'ACC'.
    :type data_type: string, optional
    :param form: Format of file of downloaded data. Zipfile containing Dyna 1.2 ascii waveforms as 'ascii', or ASDF file as 'hdf5'. Default: 'hdf5'.
    :type form: string, optional
    :param save_to: save the downloaded file on harddrive. Note: hdf5 is always saved as 'query.h5'.
    :type save_to: string, optional
    """

    if os.path.exists(tokenfile):
        tf=open(tokenfile,'r')
        dat2=tf.read()
        tf.close()
    else:
        if len(user_email)==0 and len(user_password)==0:
            raise Exception('Tokenfile does not exist and no user credentials were given.Input either path to file with token or user_email+user_passowrd!')
        else:
            dat='{"user_email":"'+user_email+'", "user_password":"'+user_password+'"}'
#data=parse.urlencode({"message":'{"user_email":"admin.isola@ingv.it","user_password":"oflodor"}'}).encode()
            data=parse.urlencode({"message":dat}).encode()
            url="https://esm-db.eu/esmws/generate-signed-message/1/query"
            resp=urlopen(url,data=data)# urlretrieve(url,data=data,filename='token.txt')
            dat2=resp.read()  #.decode('utf-8')
            f1=open(tokenfile,'w')
            f1.write(dat2.decode())
            f1.close()


#    prepare query
    url2="https://esm-db.eu/esmws/eventdata/1/query?eventid={evID}&processing-type={proctype}&data-type={dattype}&format={form}&add-xml=True&add-auxiliary-data=True".format(evID=eventID,proctype=proc_type,dattype=data_type,form=form)
    data2=parse.urlencode({'message':dat2}).encode()

    if form=="ascii":
        load_ascii(self,url2,data2,min_distance,max_distance,invert_Z_only,invert_BB_only,save_to)
    if form=="hdf5":
        load_hdf5(self,url2,data2,min_distance,max_distance,invert_Z_only,invert_BB_only,save_to)

def load_ascii(self,url,data,min_distance,max_distance,invert_Z_only,invert_BB_only,save_to): 
    import zipfile
    from BayesISOLA._input_itaca import dynaconvert,extract_event,extract_network_coordinates,load_itadata,check_quality
    if save_to:
#	    urlretrieve(url,data=data,filename=save_to)
	    zf=zipfile.ZipFile(save_to)
    else:
	    resp2=urlopen(url,data=data)
	    zzf=resp2.read()
	    zf=zipfile.ZipFile(io.BytesIO(zzf),'r')

	#pf=zipzfile.Path(zf)
	#for ipf in pf.iterdir():
    data1=Stream()
    for izf in zf.namelist():
          zf.extract(izf)
          st=dynaconvert(izf)
          for tr in st:
                    data1.append(tr)
          os.remove(izf) #remove file not to take too much space
    itadata=Stream()
    itadata=check_quality_ascii(data1)
    
    extract_event(self,itadata)
    extract_network_coordinates(self,itadata,min_distance,max_distance,invert_Z_only,invert_BB_only)
    load_itadata(self,itadata)

def load_hdf5(self,url,data,min_distance,max_distance,invert_Z_only,invert_BB_only,save_to): 
    import pyasdf #h5py
    from BayesISOLA._input_asdf import extract_event,extract_data,extract_network_coordinates
#WILL ALWAYS SAVE TO FILE
#    if save_to:
    if save_to:
        urlretrieve(url,data=data,filename=save_to)
        ds=pyasdf.ASDFDataSet(save_to)
    else:
        if os.path.exists('query.h5'):
            os.remove('query.h5')
        urlretrieve(url,data=data,filename='query.h5')
        ds=pyasdf.ASDFDataSet('query.h5')
        
    #if dataset contains auxiliary data:
    if len(ds.auxiliary_data.list())>0:
        for sta in ds.waveforms.list():
            sta=sta.replace(".","_")
            for tag in ds.waveforms[sta].get_waveform_tags():
                params=ds.auxiliary_data.Headers[sta][tag].parameters
                if 'converted' in params['processing']:
                    chann=params['stream']
                    tt=[tag2 for tag2 in ds.waveforms[sta].get_waveform_tags() if chann in ds.auxiliary_data.Headers[sta][tag2].parameters['stream'] and 'manual' in ds.auxiliary_data.Headers[sta][tag2].parameters['processing']]
                    if len(tt)==0: #there is no manually processed file
                        del ds.waveforms[sta][tag] #so remove the also the converted because it is baad
                    else:
                        for tag2 in tt:
                            del ds.waveforms[sta][tag2]
    else:
        if any(['mp' in tag for tag in ds.waveform_tags]):
        #data contains manually processed waveforms:
           for sta in ds.waveforms.list():
              for tag in ds.waveforms[sta].get_waveform_tags():
                 if tag[-2:]=='cv': #it is converted data
                    tt=tag[:-2]+'mp'
                    try:
                        del ds.waveforms[sta][tt]
                    except WaveformNotInFileException:
                        del ds.waveforms[sta][tag]

    extract_event(self,ds)
    extract_network_coordinates(self,ds,min_distance,max_distance,invert_Z_only,invert_BB_only)
    extract_data(self,ds)
    #here we can remove also remove the asdf file
    #os.remove('query.h5')

#    else:
#        resp2=urlopen(url2,data=data2)
#        rr=resp2.read()
#        hf=h5py.File(io.BytesIO(rr),'r')
#.....   hf['Waveforms']
#        hf['Provenance']
#        hf['AuxiliaryData']


