import ephem
import datetime
import numpy
##Script that will estimate the Sun altitude  in thirty minutes. If the Sun is within 0.1 radians of the horizon at that time, we should close the optical camera and shut things down.
##CJC 18/1/2015

now = datetime.datetime.utcnow() #get current time
now = datetime.datetime(2015,8,24) #get current time
##now calculate the time in 30 minutes
#ow = datetime.datetime.utcnow() #get current time
#ow = datetime.datetime(2015,01,18,17)
print('Time in UTC')
print(now)

Kuntunse=ephem.Observer()
Kuntunse.horizon='-5'
Kuntunse.lat='5.75'
Kuntunse.lon='-0.30'

galactic_center = ephem.Galactic(0, 0)
eq = ephem.Equatorial(galactic_center)
print 'RA:', eq.ra, 'dec:', eq.dec


body = ephem.FixedBody()
body._ra = eq.ra
body._dec = eq.dec
body._epoch = eq.epoch


g962 = ephem.FixedBody()
g962._ra = '18:03:15.98'
g962._dec = '-20:31:52.9'
g962._epoch = ephem.B1950
g962.compute(Kuntunse)


body.compute(Kuntunse)
print 'Az:', body.az, 'Alt:', body.alt


now = now+datetime.timedelta(0,0,2,0,5)
Kuntunse.date=now
tt= open("sunPositionData.txt","w")
ttMoon= open("moonPositionData.txt","w")
ttGalaxy= open("galaxyPositionData.txt","w")
ttg962= open("g962Data.txt","w")
##we use the time that it will be in thirty minutes to see if the Sun will be up?
for i in range(0,20000):
	now = now+datetime.timedelta(0,0,0,0,1)
 	Kuntunse.date=now
	body.compute(Kuntunse)
	g962.compute(Kuntunse)
	moontt = ephem.Moon(Kuntunse)
	suntt=ephem.Sun(Kuntunse)
	datePrint = now.strftime("%B %d %Y %H %M %S")
	aSun= "%s\t%2.2f\t%2.2f\t%2.2f\n" %(datePrint,numpy.rad2deg(suntt.az),numpy.rad2deg(suntt.alt),90-numpy.rad2deg(suntt.alt))
	aMoon= "%s\t%2.2f\t%2.2f\t%2.2f\n" %(datePrint,numpy.rad2deg(moontt.az),numpy.rad2deg(moontt.alt),90-numpy.rad2deg(moontt.alt))
	aGalaxy= "%s\t%2.2f\t%2.2f\t%2.2f\n" %(datePrint,numpy.rad2deg(body.az),numpy.rad2deg(body.alt),90-numpy.rad2deg(body.alt))
	aG962= "%s\t%2.2f\t%2.2f\t%2.2f\n" %(datePrint,numpy.rad2deg(g962.az),numpy.rad2deg(g962.alt),90-numpy.rad2deg(g962.alt))

	tt.write(aSun)
	ttMoon.write(aMoon)
	ttGalaxy.write(aGalaxy)
	ttg962.write(aG962)

tt.close()
ttMoon.close()
ttGalaxy.close()
ttg962.close()
sun=ephem.Sun()


