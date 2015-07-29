function rad_set(lststart, lstspace, lstfinish, fluxLim, filename)
% function rad_set(lststart,lstspace, lstfinish)
%
% Plot the range of elevations and az of a list of sources in rad_sources.cat
% visible starting at lststart and spaced at lstspace intervals (the
% time it takes to do one cross)
%
% Use lst>24 here-program will take care of it

clf
lat=37;       % Latitude of the site
d2r=pi/180; r2d=180/pi;
ellim=25;   % elevation limit below which we don't look
  
% Initialize plot
figure(2)
clf
polar(10,91)
hold on
view(-90,90)

% Default lststart is 0
if (~exist('lststart'))
  lststart=0;
end

% default is to run for 12 hours
if (~exist('lstfinish'))
  lstfinish=lststart+12;
end

% Default lstspace is half hour
if (~exist('lstspace'))
  lstspace=0.25;
end

% Read sources in-num is how many sources are in catalog
[names,flux, ra ,dec]=textread(filename,'%s %f %s %s %*s\n','commentstyle','matlab');

% add a + for positive declination
for i=1:size(ra,1)
  if ~strcmp( dec{i}(1),'-')
    if ~strcmp( dec{i}(1), '+')
      dec{i}=['+',dec{i}];
    end
  end
end

ind = flux>fluxLim;
names = names(ind);
flux = flux(ind);
ra = ra(ind);
dec = dec(ind);


    
num=length(ra);
disp(sprintf('%1f sources read',num))


% Convert to fractional degrees
[ra,dec]=ast2fracdeg(ra,dec);
% Convert input params in degrees and hours to radians
lat=lat*d2r;
dec=dec*d2r;

% Loop through lst range
for lstp=lststart:lstspace:lstfinish
% keep lstp in normal range
  if lstp>24
    lst=lstp-24;
  else 
    lst=lstp;
  end
  
  % Hour angle is defined as lst minus ra, here given in degrees
  ha=lst*ones(num,1)-ra/15;
  ha=15*ha*d2r; 

  % Generate az, el for sources at this lst
  cat.name=names;
  cat.flux=flux;
  [cat.az,cat.el]=hdl2ae(ha,dec,lat);
  cat.az=cat.az*r2d; cat.el=cat.el*r2d;
  
  % Only interested in sources higher than 15 degrees
  ind=cat.el>ellim;
  cat=structcut(cat,ind);
  ind=cat.flux>fluxLim;
  cat=structcut(cat,ind);
  
  % If it's the first source, take the first object
  % Otherwise look for the object furthest from the already observed
  % sources
  if (lst==lststart)
    nextsource=1;
    obs.az=cat.az(1);
    obs.el=cat.el(1);
  else
    % Method=2 means next source is the one with the furthest nearest
    % neighbor (think about it)
    nextsource=neighbor(obs,cat);
    
    % Extend obs to hold the new point
    obs.az=[obs.az;cat.az(nextsource)];
    obs.el=[obs.el;cat.el(nextsource)];
  end
  
  %    disp(sprintf('%8s,',cat.name{nextsource}))
  disp(sprintf('{%8s, %2.2f, %2.2f},',cat.name{nextsource},lst, cat.flux(nextsource)))
  polar(d2r*cat.az(nextsource),90-cat.el(nextsource),'ob')
  
end

return


function [nextsource]=furthest(obs,cat)
% function [nextsource]=furthest(obs,cat)
%
% nextsource is the index of the source in cat(az,el) that maximizes the
% sum of space angles to already observed sources in obs(az,el)

% sumsa is the sum of the space angles to all observed sources
sumsa=zeros(length(cat.az),1);
for k=1:length(obs.az)
  % Space ANgle between kth already observed source and all cat
  % sources
  ob.az=obs.az(k)*ones(length(cat.az),1);
  ob.el=obs.el(k)*ones(length(cat.az),1);
  sa=spaceangle(ob.az,ob.el,cat.az,cat.el);
  sumsa=sumsa+sa;
end

[mx,nextsource]=max(sumsa);
return

function [nextsource]=neighbor(obs,cat)
% function [nextsource]=neighbor(obs,cat)
% 
% nextsource is the index of the source that has the largest value of
% the distance to its nearest neighbor among points already observed
nearest=zeros(length(cat.az),1);

for k=1:length(cat.az)
  % Find space angle between every observed source and this option for
  % the next source and take nearest neighbor
  el=cat.el(k)*ones(length(obs.el),1);
  az=cat.az(k)*ones(length(obs.az),1);
  sa=spaceangle(az,el,obs.az,obs.el);
  nearest(k,1)=min(sa);
end

% Choose the source with the largest nearest neighbor distance
[mx,nextsource]=max(nearest);

return
