function [azErr, elErr, az, el] = trackTest(d, secOff)

  % d = read_arc('24-sep-2009:22:05:00', '25-sep-2009:02:00:00');
% d = read_arc('28-sep-2009:22:10:00', '29-sep-2009:06:24:27');
%d =read_arc('16-oct-2009:19:35:00', '16-oct-2009:20:00:00');


% the locations in servo.utc are actually from teh previous second
d.antenna0.servo.utc = d.antenna0.servo.utc+(secOff/(24*60*60));
azRec = d.antenna0.servo.fast_az_pos;
elRec = d.antenna0.servo.fast_el_pos;
azErr = d.antenna0.servo.fast_az_err;
elErr = d.antenna0.servo.fast_el_err;


% need to re-interpolate the utc for the refraction corrections.
Avals = d.antenna0.tracker.refraction(:,1);
Bvals = d.antenna0.tracker.refraction(:,2);
A = interp1(d.array.frame.utc, Avals, d.antenna0.servo.utc);
B = interp1(d.array.frame.utc, Bvals, d.antenna0.servo.utc);
% apply refraction corrections
Pvals = d.array.weather.pressure;
Tvals = d.array.weather.airTemperature;
P = interp1(d.array.frame.utc, Pvals, d.antenna0.servo.utc);
T = interp1(d.array.frame.utc, Tvals, d.antenna0.servo.utc);


% now we calculate the lst
jd = d.antenna0.servo.utc + 2400000.5;
jd0 = floor(d.antenna0.servo.utc) + 2400000.5;
H = 24*(jd - jd0);
D = jd - 2451545.0;
D0 = jd0 - 2451545.0;

long=-118.2822;
lat=37.2339;
el=1222.00;

% now we get the LST:
GMT = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H;
GMT = rem(GMT, 24);
eps = 23.4393 - 0.0000004*D;
L = 280.47 + 0.98565*D;
omega = 125.04 - 0.052954*D;
deltaPsi = -0.000319*sin(omega) - 0.000024*sin(2*L);
eqeq = deltaPsi.*cos(eps);
GAST = GMT + eqeq;
lst2 = GAST + long/15;
lst2(lst2<0) = lst2(lst2<0)+24;
lst = lst2*15;

% now for hour angle:
ra = median(d.antenna0.tracker.equat_geoc(:,1))*15; % in degrees
dec = median(d.antenna0.tracker.equat_geoc(:,2));  % in degrees
ha = lst - ra;  % also in degrees

% convert all to radians:
HA = ha*pi/180;
de = dec*pi/180;
La = lat*pi/180;

[az, el] = hdl2ae(HA,de,La);

% apply refraction correction.
R = 0.00452.*P./( (273.15+T ).*tan(el))*pi/180;
el = el+R;

az = az*180/pi;
el = el*180/pi;
az(az<0) = az(az<0)+360;
% we need to add the online pointing model applied
model = [-0.1155, -0.0976, -0.0245, -0.0277, -0.3234, 0.1846, 0, ...
	   -0.5755, -0.1497];
[az, el] = pointing_model(model, az, el);

% apply the refraction corrections to the elevation.
%el = el.*pi/180;
%sinel = sin(el);
%cosel = cos(el);

%el2 = (A.*cosel.*(sinel).^3 + B.*sinel.*cosel.^3)./( sinel.^4 + A.*sinel.^2 + 3*B.*cosel.^2 );
%el1 = el - (A.*cot(el) + B.*(cot(el)).^3);
%el2 = el2.*180/pi;

% naed to apply the encoder zeros:
%az = az + d.antenna0.tracker.encoder_off(1,1);
%el = el + d.antenna0.tracker.encoder_off(1,2);

% now for plots:

time = (lst-lst(1))/15;
clf;
figure(1)
plot( time,(elRec-el)*60, 'k');
xlabel('time');
ylabel('Error, arcminutes');
title('Elevation Error');

figure(2)
plot( time,(azRec-az)*60, 'r');
%legend('Elevation', 'Azimuth');
xlabel('time');
ylabel('Error, arcminutes');
title('Azimuth Error');


nanmean(elRec-el)*60*60
nanmean(azRec-az)*60*60


