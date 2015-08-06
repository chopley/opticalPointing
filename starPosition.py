#Script that will do the following:
#1) Take as input star name from starList.cat
#2) Output positions of the star at a given time +-30 seconds
#Written by Charles Copley,AVN Science


#Rev 1.0 06/08/2015 


import pandas,numpy,ephem,sys
#define the column layout

#starName = sys.argv[1]
starName = 'KapCen'
#date = sys.argv[3]
date = '2015/08/02'
#siteName = sys.argv[2]
siteName = 'Klerefontein'

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

elevation=siteDetail['Altitude']
print Kuntunse.lat,Kuntunse.lon,elevation


#and loop through and calculate everything
source = ephem.FixedBody()

indexLoc = starList[starList['Name']==starName].index.tolist()

source._ra = starList['RA'][indexLoc[0]]
source._dec = starList['DEC'][indexLoc[0]]
source._epoch=ephem.J2000
now = datetime.datetime(2015,8,02,16,43,55) #set current time
Kuntunse.date=now


source.compute(Kuntunse)
print numpy.rad2deg(source.az),numpy.rad2deg(source.alt)







