function [az, el] = calcAzEl(d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [az, el] = calcAzEl(d, secOff)
%
%  calculates the azimuth and elevation of the source to be tracked.
%
%  sjcm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %
utc = d.antenna0.receiver.utc;

% need to re-interpolate the utc for the refraction corrections.
Avals = d.antenna0.tracker.refraction(:,1);
Bvals = d.antenna0.tracker.refraction(:,2);
A = interp1(d.array.frame.utc, Avals, utc,'spline');
B = interp1(d.array.frame.utc, Bvals, utc,'spline');
% apply refraction corrections
Pvals = d.array.weather.pressure;
Tvals = d.array.weather.airTemperature;
P = interp1(d.array.frame.utc, Pvals, utc,'spline');
T = interp1(d.array.frame.utc, Tvals, utc,'spline');

long=-118.2822;
lat=37.2339;
el=1222.00;

% now we calculate the lst
%jd = utc + 2400000.5;
%jd0 = floor(utc) + 2400000.5;
%H = 24*(jd - jd0);
%D = jd - 2451545.0;
%D0 = jd0 - 2451545.0;

% now we get the LST:
%GMT = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H;
%GMT = rem(GMT, 24);
%eps = 23.4393 - 0.0000004*D;
%L = 280.47 + 0.98565*D;
%omega = 125.04 - 0.052954*D;
%deltaPsi = -0.000319*sin(omega) - 0.000024*sin(2*L);
%eqeq = deltaPsi.*cos(eps);
%GAST = GMT + eqeq;
%lst2 = GAST + long/15;
%lst2(lst2<0) = lst2(lst2<0)+24;
%lst = lst2*15;
lst = d.antenna0.tracker.lst*15;
indpos = lst > 180;
lst(indpos) = lst(indpos)-360;


lst = interp1(d.array.frame.utc, lst, utc); % in degrees
indpos = interp1(d.array.frame.utc, indpos, utc, 'nearest');
indpos(isnan(indpos)) = 0;
indpos = logical(indpos);
lst(indpos) = lst(indpos)+360;

% now for hour angle:
ra = interp1(d.array.frame.utc, d.antenna0.tracker.equat_geoc(:,1), utc);
ra = ra*15; % in degrees
dec = interp1(d.array.frame.utc, d.antenna0.tracker.equat_geoc(:,2), utc); % in degrees
ha = lst - ra;  % also in degrees

% convert all to radians:
HA = ha*pi/180;
de = dec*pi/180;
La = lat*pi/180;

[az, el] = hdl2ae(unwrap(HA),de,La);

% apply refraction correction.
R = 0.00452.*P./( (273.15+T ).*tan(el))*pi/180;
el = el+R;

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
