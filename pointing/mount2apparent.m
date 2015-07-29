function d=mount2apparent(d)
%  function d=mount2apparent(d)
%
% Convert ideal iaz,iel to model maz,mel using pointing model
% parameters tm

d2r=pi/180;

% first we get the values of all our parameters, etc etc.
mount_az = d.antenna0.servo.az*d2r;
mount_el = d.antenna0.servo.el*d2r;

model(:,1)=d.antenna0.tracker.flexure(:,1);
model(:,2)=d.antenna0.tracker.flexure(:,2);
model(:,3)=d.antenna0.tracker.tilts(:,1);
model(:,4)=d.antenna0.tracker.tilts(:,2);
model(:,5)=d.antenna0.tracker.tilts(:,3);
model(:,6)=d.antenna0.tracker.fixedCollimation(:,1);
model(:,7)=d.antenna0.tracker.fixedCollimation(:,2);
model(:,8)=d.antenna0.tracker.encoder_off(:,1);
model(:,9)=d.antenna0.tracker.encoder_off(:,2);

model = model*d2r;

% now we work our way back.

% refraction correction.
% need to re-interpolate the utc for the refraction corrections.
Avals = d.antenna0.tracker.refraction(:,1);
Bvals = d.antenna0.tracker.refraction(:,2);
A = interp1(d.array.frame.utc, Avals, d.antenna0.receiver.utc);
B = interp1(d.array.frame.utc, Bvals, d.antenna0.receiver.utc);
% apply refraction corrections
Pvals = d.array.weather.pressure;
Tvals = d.array.weather.airTemperature;
P = interp1(d.array.frame.utc, Pvals, d.antenna0.receiver.utc);
T = interp1(d.array.frame.utc, Tvals, d.antenna0.receiver.utc);

R = 0.00452.*P./( (273.15+T ).*tan(mount_el))*pi/180;
el = mount_el - R;

% Encoder zero points
az = mount_az - model(8);
el = el - model(9);

% El collimation
% Exactly the same as el zero point shift
el = el - model(7);

% Cross-el Collimation
el = asin(sin(el).*cos(model(6)));
el(imag(el)~=0) = nan;
az = az + asin(-sin(model(6))./cos(el));
az(imag(az)~=0) = nan;


% El tilt
el = asin(sin(el).*cos(model(5)));
el(imag(el)~=0)=NaN;
az = az + asin(tan(model(5)).*sin(el)./cos(el));
az(imag(az)~=0)=NaN;

% Az tilt

% astronomy az is clockwise from north
% but want to work in frame with az anticlock from x as for matlab
% sph2cart function etc. Also x-east, y-north, z-up seems most
% natural.
az = -az + pi/2;
% Get normal vector to tilt plane
c=cross([1,0,tan(model(3))],[0,1,tan(-model(4))]);
% Find magnitude and dir of tilt
phi=atan2(c(2),c(1));
theta=atan(sqrt(c(1)^2+c(2)^2)/c(3));
% Apply rotation
[x,y,z]=sph2cart(az,el,ones(size(az)));
[x,y,z]=rotaboutz(x,y,z,+phi);  % rotate to x along tilt dir
[x,y,z]=rotabouty(x,y,z,-theta); % rotate by negative tilt angle
[x,y,z]=rotaboutz(x,y,z,-phi);   % rotate back
[az,el]=cart2sph(x,y,z);
% Convert back to az clock from north
az = -az + pi/2;


% Flexure
toc
tic 
for m=1:length(el)
  el(m) = fzero(@(x) (x-model(2)*cos(x) - el(m)), el(m));
  el(m) = fzero(@(x) (x-model(1)*sin(x) - el(m)), el(m));
end
toc

% and now, refraction corrections.


return
