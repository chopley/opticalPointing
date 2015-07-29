function [az, el] = calcAzElCs(d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [az, el] = calcAzElCs(d)
%
%  calculates the azimuth and elevation of the source to be tracked.
%
%  sjcm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definitions
long = d.antenna0.tracker.siteActual(1,1);
lat  = d.antenna0.tracker.siteActual(1,2);
elev = d.antenna0.tracker.siteActual(1,3);

% LST of data
lst    = d.antenna0.tracker.lst*15;
utc    = d.antenna0.receiver.utc;
indpos = lst > 180;
lst(indpos) = lst(indpos)-360;
lst    = interp1(d.array.frame.utc, lst, utc); % in degrees
indpos = interp1(d.array.frame.utc, single(indpos), utc, 'nearest');
indpos(isnan(indpos)) = 0;
indpos = logical(indpos);
lst(indpos) = lst(indpos)+360;

% start with the RA/DEC of the source
ra = interp1(d.array.frame.utc, d.antenna0.tracker.equat_geoc(:,1), utc);
ra = ra*15; % in degrees
dec = interp1(d.array.frame.utc, d.antenna0.tracker.equat_geoc(:,2), utc); % in degrees
ha = lst - ra;  % also in degrees

% convert all to radians:
HA = ha*pi/180;
de = dec*pi/180;
La = lat*pi/180;

% now we convert to az/el
[az, el] = hdl2ae(unwrap(HA),de,La);
% hdl2ae is equivalent to the below, which is in the control system.
%el1 = asin( sin(de).*sin(La) + cos(de).*cos(La).*cos(HA));
%az1 = atan2( -cos(de).*sin(HA) , (sin(de).*cos(La) - ...
%    cos(de).*sin(La).*cos(HA)));

% next we apply a horizontal parallax, if necessary
sources = unique(d.antenna0.tracker.source);
ephemSourceList = {'moon', 'sun', 'jupiter', 'mars', 'mercury', 'neptune', ...
    'pluto', 'saturn', 'uranus', 'venus'};
needCorr = 0;
for m=1:length(sources)
  if(~isempty(strmatch(sources{m}, ephemSourceList)))
    srcIndex  = strmatch(sources{m}, ephemSourceList);
    needCorr  = 1;
    indSource = zeros(size(d.antenna0.tracker.source,1), 1);
    f         = strmatch(ephemSourceList{srcIndex}, d.antenna0.tracker.source(:,1));
    indSource(f) = 1;
    srcName   = ephemSourceList{srcIndex};
  end
end

% if there is a source name, we need to read in the ephem and interpolate
% accordingly
if(needCorr)
  [dist err r_cent] = cbassMatEphem(srcName, d.antenna0.receiver.utc);
  
  % apply the correction
  el = el - atan2( cos(el), dist/r_cent - sin(el));
end

% next we correct for refraction
% need to re-interpolate the utc for the refraction corrections.
Avals = d.antenna0.tracker.refraction(:,1);  % degrees
Bvals = d.antenna0.tracker.refraction(:,2);  % degress
refraVal = d.antenna0.tracker.refraction(:,3); % degrees

A = interp1(d.array.frame.utc, Avals, utc) * pi/180;
B = interp1(d.array.frame.utc, Bvals, utc) * pi/180;

refracCorr = ( A.*cos(el).*(sin(el)).^3 + B.*sin(el).*(cos(el)).^3 ) ./ ( ...
    (sin(el)).^4 + A.*(sin(el)).^2 + 3*B.*(cos(el)).^2);
el = el + refracCorr;

% next the diurnal aberration
%[az, el] = applyDiurnalAberration(az, el);

% now this is the azapparent and el apparent.

% next we convert back to degrees and apply the pointing model corrections.
az = az*180/pi;
el = el*180/pi;
az(az<0) = az(az<0)+360;
% we need to add the online pointing model applied
model = [d.antenna0.tracker.flexure(1,:), d.antenna0.tracker.tilts(1,:), ...
      d.antenna0.tracker.fixedCollimation(1,:), ...
      d.antenna0.tracker.encoder_off(1,:)];

% apply the pointing model
[az, el] = pointing_model(model, az, el);

return;
