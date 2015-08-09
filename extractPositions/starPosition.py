#Script that will do the following:
#1) Read in positions of stars from images in pixels
#2) Use a catalog to calculate the expected positions of the stars
#3) Calculate the az,el position of the centre of the image

#Written by Charles Copley,AVN Science


#Rev 1.0 06/08/2015 


import pandas,numpy,ephem,sys,cv2,datetime
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
	return (dCol,dRow,snr,sourceCol,sourceRow)

def estimateStarPosition(timeObserved,starName,siteDetail,starList):
	Kuntunse=ephem.Observer()
	Kuntunse.lat=(siteDetail['Lat'])*ephem.degree
	Kuntunse.lon=(siteDetail['Lon'])*ephem.degree
	elevation=siteDetail['Altitude']
	#print Kuntunse.lat,Kuntunse.lon,elevation
	indexLoc = starList[starList['Name']==starName].index.tolist()	
	source = ephem.FixedBody()
	source._ra = starList['RA'][indexLoc[0]]
	source._dec = starList['DEC'][indexLoc[0]]
	source._epoch=ephem.J2000
	Kuntunse.date=timeObserved
	source.compute(Kuntunse)
#	a='%s %s %2.2f %2.2f' %(starName,timeObserved,numpy.rad2deg(source.az),numpy.rad2deg(source.alt))
	return (timeObserved,source.az.norm/ephem.degree,source.alt/ephem.degree)

def cartRotate(x,y,angle):
	x2 = x*numpy.cos(angle) + y*numpy.sin(angle)
	y2 = -x*numpy.sin(angle) + y*numpy.cos(angle)
	return (x2,y2)	

def analyseImage(image1,image2,rotAngle,interpValues):
	az1 = image1[3]
	az2 = image2[3]
	dAz=az2-az1
	if(dAz>=180):
		dAz=dAz-360
	if(dAz<=-180):
		dAz=dAz+360
	el1 = image1[4]
	el2 = image2[4]
	dEl=el2-el1
	xpix1=image1[1]
	xpix2=image2[1]
	ypix1=image1[2]
	ypix2=image2[2]
	dXpix = float(xpix2-xpix1)
	dYpix = float(ypix2-ypix1)
	distanceAngle = numpy.sqrt(((dAz)*numpy.cos(numpy.radians(el2)))**2 + (dEl)**2)
	distancePixel = numpy.sqrt(((dXpix))**2 + (dYpix)**2)
	rotationAngle=numpy.radians(rotAngle)
	
	(xp1,yp1) =  cartRotate(xpix1,ypix1,rotationAngle)
	(xp2,yp2) =  cartRotate(xpix2,ypix2,rotationAngle)
	
	i1=(xp1,yp1)
	i2=(xp2,yp2)	
	distanceScale = distancePixel/distanceAngle

	zeroAz = interpValues[0]
	zeroEl = interpValues[1]
	return (i1,i2,distanceScale)

#values that require setting	
#Angle of the azimuth axis to the Image horizontal axis
rotationAngle=14.4
#Number of pixels to one degree
Scale=160.
	
#starName = sys.argv[1]
starName = 'KapCen'
#date = sys.argv[3]
date = '2015/08/02'
#siteName = sys.argv[2]
siteName = 'Klerefontein'

#read in the list of possible observing sites
sites = pandas.read_csv('sites.txt')
siteDetail=sites[sites['Site']==siteName]
print siteDetail

#file that contains the star name, and time of observations
starDetail = pandas.read_csv('starLog.txt')


colspecs=[(0,4),(6,15),(16,26),(28,40),(41,53),(54,65),(66,71)]
#read in the data frame from the starlist catalog
starList=pandas.read_fwf('starList.cat',colspecs=colspecs)
	

starDetail['Time2']='na'
starDetail['Az']='na'
starDetail['El']='na'
starDetail['dPixCol']='na'
starDetail['dPixRow']='na'
starDetail['SNR']='na'
starDetail['PixCol']='na'
starDetail['PixRow']='na'
starDetail['Image']='na'
starDetail['PixAz']='na'
starDetail['PixEl']='na'
starDetail['offAz']='na'
starDetail['offEl']='na'
starDetail['centreAz']='na'
starDetail['centreEl']='na'
for i in range(0,len(starDetail)):
	obsTime1=starDetail['Time'][i]
	obsTime2= datetime.datetime.strptime(obsTime1,'%Y-%m-%d-%H%M%S')
	starName=starDetail['StarName'][i]
	(timeout1,az1,el1)= estimateStarPosition(obsTime2,starName,siteDetail,starList)
	(pixcol1,pixrow1,snr1,sCol1,sRow1) = extractImage(obsTime1)
	starDetail['Time2'][i]=obsTime2	
	starDetail['Az'][i]=az1	
	starDetail['El'][i]=el1	
	starDetail['dPixCol'][i]=pixcol1	
	starDetail['dPixRow'][i]=pixrow1	
	starDetail['SNR'][i]=snr1	
	starDetail['PixCol'][i]=sCol1	
	starDetail['PixRow'][i]=sRow1	
	image1=(starDetail['StarName'][i],starDetail['dPixCol'][i],starDetail['dPixRow'][i],starDetail['Az'][i],starDetail['El'][i],starDetail['PixCol'][i],starDetail['PixCol'][i])
	starDetail['Image'][i]=image1
	(xp1,yp1) =  cartRotate(starDetail['dPixCol'][i],starDetail['dPixRow'][i],rotationAngle)	
	starDetail['PixAz'][i]=xp1
	starDetail['PixEl'][i]=yp1
	starDetail['offAz'][i]=xp1/Scale
	starDetail['offEl'][i]=yp1/Scale
	starDetail['centreEl'][i]=float(starDetail['El'][i])-float(starDetail['offEl'][i])
	starDetail['centreAz'][i]=float(starDetail['Az'][i])+float(starDetail['offAz'][i])/numpy.cos(numpy.radians(80))
	if(float(starDetail['centreAz'][i])>180.):
		starDetail['centreAz'][i]=starDetail['centreAz'][i]-360
	


	

rotAngle=14.4
print starDetail['StarName'][5], starDetail['StarName'][4]
#first just check the rotation angle is correct- 
# First make sure the two images are just azimuth movement
# we want i1a[0]-i2a[0] to be close to zero. This implies zero elevation motion
(i1a,i2a,distanceScalea)= analyseImage(starDetail['Image'][5],starDetail['Image'][4],rotAngle,(180,80))
print i1a[0]-i2a[0],i1a[0]-i1a[1]
print distanceScalea

(i1,i2,distanceScale)= analyseImage(starDetail['Image'][0],starDetail['Image'][5],rotAngle,(180,80))
print (i1[0]-i2[0])/distanceScale,(i1[0]-i1[1])/distanceScale
print distanceScale

(i1,i2,distanceScale)= analyseImage(starDetail['Image'][8],starDetail['Image'][6],rotAngle,(0,80))
print (i1[0]-i2[0])/distanceScale,(i1[0]-i1[1])/distanceScale
print distanceScale

(i1,i2,distanceScale)= analyseImage(starDetail['Image'][9],starDetail['Image'][6],rotAngle,(0,80))
print (i1[0]-i2[0])/distanceScale,(i1[0]-i1[1])/distanceScale
print distanceScale

(i1,i2,distanceScale)= analyseImage(starDetail['Image'][7],starDetail['Image'][6],rotAngle,(0,80))
print (i1[0]-i2[0])/distanceScale,(i1[0]-i1[1])/distanceScale
print distanceScale




print starDetail





