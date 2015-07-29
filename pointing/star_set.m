function star_set(lst,ellim,maglim,lat, minuteDiff)
% star_set(lst)
%
% Find the set of stars which are bright and available for pointing
% at a given LST
%
% To get the lst type "print $time(lst)" in szaViewer
%
% e.g.: star_set(19.1)
%
% Take the resulting list of numbers and paste it into
% /home/szadaq/carma/sza/array/sch/point_manual.sch
%
% Then run using:
% schedule /home/szadaq/carma/sza/array/sch/point_manual.sch
%
% Then center up each star using widget and hit the button
% Be sure to keep a note of which data file in
% /home/szadaq/arc contains the data from the star pointing run

if(~exist('ellim'))
  ellim=[];
end
if(~exist('maglim'))
  maglim=[];
end
if(~exist('lat'))
  lat=[];
end

if(isempty(ellim))
  ellim=30;
end
if(isempty(maglim))
  maglim=4.0;
end
if(isempty(lat))
  lat=37; % OVRO
end


% Read in star catalog data
[cat.name,ra,dec,cat.mag] = ...
  textread('stars.cat','J2000 %s %s %s %*s %*s # %f','commentstyle','matlab');

% Convert ra,dec to fractional degrees
[cat.ra,cat.dec]=ast2fracdeg(ra,dec);

% Calc hour angle for given lst
cat.ha=lst*15-cat.ra;

% Calc az,el for site at given latitude
d2r=pi/180; r2d=180/pi;
[cat.az,cat.el]=hdl2ae(cat.ha*d2r,cat.dec*d2r,lat*d2r);
cat.az=cat.az*r2d; cat.el=cat.el*r2d;

% Cut down to stars above elevation and magnitude limits
ind=cat.el>ellim&cat.mag<maglim;
cat=structcut(cat,ind);

% Cut out stars close to zenith as causes hangup
ind=cat.el<83;
cat=structcut(cat,ind);

% Move az<-90 sources to +ve
ind=cat.az<-90;
cat.az(ind)=cat.az(ind)+360;

% Sort to decreasing azimuth order to make for an efficent
% star pointing run
[x,ind]=sort(cat.az);
cat=structcut(cat,flipud(ind));

% Print list for inclusion into pointing script
junk = lst:minuteDiff/60:(lst + minuteDiff/60*length(cat.az));
for i=1:length(cat.az)
  disp(sprintf('{%8s, %3.2f}, # %7.2f %6.2f %5.2f',...
	       cat.name{i},junk(i),cat.az(i),cat.el(i),cat.mag(i)));
end

