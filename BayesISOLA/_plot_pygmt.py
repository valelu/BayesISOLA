import pygmt

def plot_stations_pygmt(self,outfile='$outdir/stations.png',network=True,location=False,channelcode=False):
    """
	Plot a map of stations used in the inversion.
	
	:param outfile: path to file for plot output; if ``None`` plots to the screen
	:type outfile: string, optional
	:param network: include network code into station label
	:type network: bool, optional
	:param location: include location code into station label
	:type location: bool, optional
	:param channelcode: include channel code into station label
	:type channelcode: bool, optional
	:param fontsize: font size for all texts in the plot; if zero, the size is chosen automatically
	:type fontsize: scalar, optional
	
	
	The stations are marked according to components used in the inversion.
    """
    fig=pygmt.Figure()
    #define region
    #get station coordinates
    lats= [sta['lat'] for sta in self.inp.stations]
    longs=[sta['lon'] for sta in self.inp.stations]
    #epicenter
    lats.append(self.MT.centroid['lat'])
    longs.append(self.MT.centroid['lon'])
    region=[min(longs)-.5,max(longs)+.5,min(lats)-.5,max(lats)+.5]
    projection='M12i' #Mercator projection with 12inch size 
    #focal mechanism?
    if self.MT.decompose:
        focal=dict(strike=self.MT.centroid['s1'],dip=self.MT.centroid['d1'],rake=self.MT.centroid['r1'],magnitude=self.MT.centroid['Mw'])
    #get labels for colors...
    cols=[];stalabs=[]
    for sta in self.inp.stations:
        if sta['useN'] and sta['useE'] and sta['useZ']:
            col=2#'red' #red
        elif not sta['useN'] and not sta['useE'] and not sta['useZ']:
            col=0#'white' #white
        else:
            col=1#'gray' #gray
        cols.append(col)
        if network and sta['network']: l=sta['network']+':'
        else: l = ''
        l +=sta['code']
        if location and sta['location']: l+= ':'+sta['location']
        if channelcode: l+=' '+sta['channelcode']
        stalabs.append(l)

    fig.basemap(region=region,projection=projection,frame='a1')
    #get topo data:
    grid=pygmt.datasets.load_earth_relief(resolution='01m',region=region)
    fig.grdimage(grid=grid,cmap='geo',transparency=40)
    fig.coast(shorelines='1/2p',lakes='azure',borders='1p')
    #plot beachball:
    if self.MT.decompose:
        fig.meca(focal,scale='1i',longitude=longs[-1],latitude=lats[-1],depth=self.MT.centroid['z']/1.e3)
    else:
        fig.plot(x=longs[-1],y=lats[-1],style='a',color='yellow',pen='1p,black',size='1i')
    #plot stations:
    pygmt.makecpt(cmap='white,orange,red',series='0/2/1',color_model="+cnot used,some components,all used")
    fig.plot(x=longs[:-1],y=lats[:-1],style='t0.7i',color=cols,cmap=True,pen='1p,black')
    fig.text(text=stalabs,x=longs[:-1],y=lats[:-1],offset='0./1.',font='20p')
    fig.colorbar()
    if outfile:
        outfile=outfile.replace('$outdir',self.outdir)
        fig.savefig(outfile)

