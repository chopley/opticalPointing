#This is a script that reads in the starlist catalog
#It will eventually do the following:
#1) Start at a given time, latitude, longitude, elevation, azimuth, field of view
#2) Find the next star that will pass through the elevation,azimuth, as well as the time it will pass through
#3) Output the star, and the time it will pass through transit
#4) Grab images from the frame grabber for <TBD> time before and after zenith transit
#Written by Charles Copley,AVN Science

#Call as follows
# python findStars.py 90 2 2015/8/21 Kuntunse South (where 90 is the elevation position and 2 is the field of view of the camera)
# This will return a list of appropriate stars for the Ghana Dish, RA, DEC, time of transit, elevation angle at transit

#Rev 1.0 29/7/2015 
#TO DO: 
#1) Add appropriate functionality for location and time/date (probably a text file with location names and lat/long/altitude).
#2) Also work out if the transit is North or South facing (pretty simple)
#3) Add function to take images with the frame grabber during the transit?
#4) 

#Rev 2.0 29/7/2015 
#Done
#1) Added a text file (sites.txt) that is read to get appropriate details for the site location.
#2) Added a North/South option to allow you to select which side of zenith the star will be

#TO DO: 
#1) Also work out if the transit is North or South facing (pretty simple)
#2) Add function to take images with the frame grabber during the transit?


import pandas,numpy,ephem,sys
#define the column layout

elevation = float(sys.argv[1]) #first argument is the elevation position we are looking at
#elevation = 60. #first argument is the elevation position we are looking at
fieldOfView = float(sys.argv[2]) #second argument is the field of view size
#fieldOfView = 2. #second argument is the field of view size
date = sys.argv[3]
#date = '2015/08/22'
siteName = sys.argv[4]
#siteName = "Klerefontein"
direction = sys.argv[5]




#read in the list of sites
sites = pandas.read_csv('sites.txt')
siteDetail=sites[sites['Site']==siteName]
print siteDetail

colspecs=[(0,4),(6,15),(16,26),(28,40),(41,53),(54,65),(66,71)]
#read in the data frame from the starlist catalog
starList=pandas.read_fwf('starList.cat',colspecs=colspecs)

starList['DECd']=""
starList['SkyAngle']=""
starList['TransitTime']=""

#Set up the observing site using the information from the text file
Kuntunse=ephem.Observer()
Kuntunse.lat=(siteDetail['Lat'])*ephem.degree
Kuntunse.lon=(siteDetail['Lon'])*ephem.degree
#get the date from the input argument earlier 
Kuntunse.date = date 

print Kuntunse.lat,Kuntunse.lon,elevation


#and loop through and calculate everything
source = ephem.FixedBody()
fg = open('StarPositions.txt','w')
for i in range(0,len(starList)):
	source._ra = starList['RA'][i]
	source._dec = starList['DEC'][i]
	source._epoch=ephem.J2000
	source.compute(Kuntunse)
	Direction = (source._dec-Kuntunse.lat)
	Direction = (numpy.abs(Direction)/Direction)
	DirectionStr='North'
	if(Direction<0):	
		DirectionStr='South'
	#Need to check if the source will ever be up
	if(~source.neverup):
		try:
			zenithAngle = numpy.rad2deg(numpy.deg2rad(90) - source.transit_alt.norm	)
			outputString = '%s\t%s\t%s\t%s\t%s\t%s\t%2.2f' %(starList['Name'][i],starList['RA'][i],starList['DEC'][i],source.transit_time,source.transit_alt.norm,DirectionStr,zenithAngle)
			if(source.transit_alt.norm>numpy.radians(elevation-fieldOfView) and  (source.transit_alt.norm<numpy.radians(elevation+fieldOfView))):
				if(DirectionStr==direction):
					print outputString
					print aa
					fg.write('hello')
					fg.write('\n')
		except:
			outputString = '%s %s %s' %(starList['Name'][i],starList['RA'][i],starList['DEC'][i])

print 'closing the file'
fg.close()




