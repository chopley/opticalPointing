#This is a script that reads in the starlist catalog
#It will eventually do the following:
#1) Start at a given time, latitude, longitude, elevation, azimuth, field of view
#2) Find the next star that will pass through the elevation,azimuth, as well as the time it will pass through
#3) Output the star, and the time it will pass through

#Rev 1.0 29/7/2015 
#Call as follows
# python findStars.py 90 2 (where 90 is the elevation position and 2 is the field of view of the camera)
# This will return a list of appropriate stars for the Ghana Dish, RA, DEC, time of transit, elevation angle at transit
#TO DO: 
#1) Add appropriate functionality for location (probably a text file with location names and lat/long/altitude).
#2) Also work out if the transit is North or South facing (pretty simple)

#Written by Charles Copley 28/7/2015
import pandas,numpy,ephem,sys
#define the column layout

elevation = float(sys.argv[1]) #first argument is the elevation position we are looking at
elevationRadians = numpy.deg2rad(float(elevation))
azimuth = 180
latitude = '5.75'
longitude = '-0.30'
fieldOfView = float(sys.argv[2]) #second argument is the field of view size



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




