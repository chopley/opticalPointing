#Script that will do the following:
#1) Take as input star name from starList.cat
#2) Output positions of the star at a given time +-30 seconds
#Written by Charles Copley,AVN Science


#Rev 1.0 06/08/2015 


import pandas,numpy,ephem,sys,cv2
#define the column layout

def extractImage(time):
	#function that reads an image at a given time, and outputs the offset of the brightest point in pixels
	filename=time+'.png'
	img = cv2.imread(filename, 0)
	rows,cols  = img.shape
	#remove some of the white stuff in the image
	img[20:rows-1,0:80]=0
	img[0:23]=0
	centerRow, centerCol = rows/2, cols/2
	#Find the brightest pixel position
	(sourceRow,sourceCol) = numpy.unravel_index(img.argmax(), img.shape)
	dCol = sourceCol-centerCol
	dRow = sourceRow-centerRow
	snr = numpy.max(img)/numpy.mean(img)
	return (dCol,dRow,snr)

def estimateStarPosition(timeObserved,starName,siteDetail,starList):
	Kuntunse=ephem.Observer()
	Kuntunse.lat=(siteDetail['Lat'])*ephem.degree
	Kuntunse.lon=(siteDetail['Lon'])*ephem.degree
	elevation=siteDetail['Altitude']
	print Kuntunse.lat,Kuntunse.lon,elevation
	indexLoc = starList[starList['Name']==starName].index.tolist()	
	source = ephem.FixedBody()
	source._ra = starList['RA'][indexLoc[0]]
	source._dec = starList['DEC'][indexLoc[0]]
	source._epoch=ephem.J2000
	Kuntunse.date=timeObserved
	source.compute(Kuntunse)
	a='%s %2.2f %2.2f' %(timeObserved,numpy.rad2deg(source.az),numpy.rad2deg(source.alt))
	print a
	return (timeObserved,source.az.norm/ephem.degree,source.alt/ephem.degree)

	
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
	

starName1 = 'KapCen'
time1='2015-08-02-164753'
time2='2015-08-02-164537'
time3='2015-08-02-164355'

starName2 = 'DelLup'
time4='2015-08-02-170957'
time5 = '2015-08-02-171031'
time6 = '2015-08-02-171139'

time1S= datetime.datetime.strptime(time1,'%Y-%m-%d-%H%M%S')
time2S= datetime.datetime.strptime(time2,'%Y-%m-%d-%H%M%S')
time3S= datetime.datetime.strptime(time3,'%Y-%m-%d-%H%M%S')
time4S= datetime.datetime.strptime(time4,'%Y-%m-%d-%H%M%S')
time5S= datetime.datetime.strptime(time5,'%Y-%m-%d-%H%M%S')
time6S= datetime.datetime.strptime(time6,'%Y-%m-%d-%H%M%S')

(timeout1,az1,el1)= estimateStarPosition(time1S,starName1,siteDetail,starList)
(timeout2,az2,el2)= estimateStarPosition(time2S,starName1,siteDetail,starList)
(timeout3,az3,el3)= estimateStarPosition(time3S,starName1,siteDetail,starList)
(pixcol1,pixrow1,snr1) = extractImage(time1)
(pixcol2,pixrow2,snr2) = extractImage(time2)
(pixcol3,pixrow3,snr3) = extractImage(time3)


(timeout4,az4,el4)= estimateStarPosition(time4S,starName2,siteDetail,starList)
(timeout5,az5,el5)= estimateStarPosition(time5S,starName2,siteDetail,starList)
(timeout6,az6,el6)= estimateStarPosition(time6S,starName2,siteDetail,starList)


(pixcol4,pixrow4,snr4) = extractImage(time4)
(pixcol5,pixrow5,snr5) = extractImage(time5)
(pixcol6,pixrow6,snr6) = extractImage(time6)




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
#ow = datetime.datetime(2015,8,02,16,43,55) #set current time
#ow = source.transit_time.datetime() #set current time
#ow = now-datetime.timedelta(seconds=150)

Kuntunse.date=timett
source.compute(Kuntunse)
a='%s %2.2f %2.2f' %(timett,numpy.rad2deg(source.az),numpy.rad2deg(source.alt))
print a







