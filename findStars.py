#This is a script that reads in the starlist catalog
#It will eventually do the following:
#1) Start at a given time, latitude, longitude, elevation, azimuth, field of view
#2) Find the next star that will pass through the elevation,azimuth, as well as the time it will pass through
#3) Output the star, and the time it will pass through

#Written by Charles Copley 28/7/2015
import pandas,numpy,ephem
#define the column layout

elevation = 90
elevationRadians = numpy.deg2rad(elevation)
azimuth = 180
latitude = '5.75'
longitude = '-0.30'
fieldOfView = 2#size of field of view in degrees



colspecs=[(0,4),(6,11),(16,26),(28,40),(41,53),(54,65),(66,71)]
#read in the data frame from the starlist catalog
starList=pandas.read_fwf('starList.cat',colspecs=colspecs)

starList['DECd']=""
starList['SkyAngle']=""
starList['TransitTime']=""

#Set up the observing site
Kuntunse=ephem.Observer()
Kuntunse.lat='5.75'
Kuntunse.lon='-0.30'
Kuntunse.date = '2015/8/21'

#and loop through and calculate everything
source = ephem.FixedBody()
for i in range(0,len(starList)):
	source._ra = starList['RA'][i]
	source._dec = starList['DEC'][i]
	source._epoch=ephem.J2000
	source.compute(Kuntunse)
	#starList['TransitTime'][i] = ephem.Date.datetime(Kuntunse.next_transit(source))
	
	aa = '%s %s %s %s %s' %(starList['Name'][i],starList['RA'][i],starList['DEC'][i],source.transit_time,source.transit_alt.norm)
	if(source.transit_alt.norm>numpy.radians(elevation-fieldOfView) and  (source.transit_alt.norm<numpy.radians(elevation+fieldOfView))):
		print aa

	#print(Kuntunse.next_transit(source))
	#print(source.transit_alt.norm)
	#rint(source._dec)




