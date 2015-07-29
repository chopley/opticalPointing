function [maz,mel]=pointing_model(model,iaz,iel)
% [maz,mel]=pointing_model(model,iaz,iel)
%
% Convert ideal iaz,iel to model maz,mel using pointing model
% parameters tm

if(length(model)==9)
  model(10) = 0;
  model(11) = 0;
end


% All input params in degrees - convert to rad
d2r=pi/180;
model=model*d2r;
az=iaz*d2r;
el=iel*d2r;

% cosine of azimuth
%az = az + model(10).*sin(az);
%el = el + model(11).*cos(az);


% Flexure
el=el-model(1)*sin(el);
el=el-model(2)*cos(el);

% Az tilt
%display('New Az Tilt');
%[az1, el1] = azTiltCheck(model, az, el);
% the two are equivalent.  they're off by about a tenth of an arcsecond.

% astronomy az is clockwise from north
% but want to work in frame with az anticlock from x as for matlab
% sph2cart function etc. Also x-east, y-north, z-up seems most
% natural.
az=-az+pi/2;

% Get normal vector to tilt plane
c=cross([1,0,tan(model(3))],[0,1,tan(-model(4))]);
% Find magnitude and dir of tilt
phi=atan2(c(2),c(1));
theta=atan(sqrt(c(1)^2+c(2)^2./c(3)));
% Apply rotation
[x,y,z]=sph2cart(az,el,ones(size(az)));
[x,y,z]=rotaboutz(x,y,z,phi);   % rotate to x along tilt dir
[x,y,z]=rotabouty(x,y,z,theta); % rotate by tilt angle
[x,y,z]=rotaboutz(x,y,z,-phi);  % rotate back
[az,el]=cart2sph(x,y,z);

% Convert back to az clock from north
az=-az+pi/2;



% El tilt
% There is no way to do this with vector rotation as axes not perp.
% Below is taken from Tim/Martin
el = asin(sin(el)./cos(model(5)));
el(imag(el)~=0)=NaN;
az=az-asin(tan(model(5)).*sin(el)./cos(el));
az(imag(az)~=0)=NaN;

% Cross-el Collimation
   
% There is no way to do cross-el collimation with vector rotation
% Below is taken from Tim/Martin for the case of deck=270
az=az-asin(-sin(model(6))./cos(el));
az(imag(az)~=0)=NaN;
el=asin(sin(el)./cos(model(6)));
el(imag(el)~=0)=NaN;

% El collimation
% Exactly the same as el zero point shift
el=el+model(7);


% term that's proportional to cos(az)
%el = el + model(10)*cos(az);
%el = el + model(11)*sin(az);
%az = az + model(12)*sin(az);


% Encoder zero points
az=az+model(8);
el=el+model(9);


% put phi into 0 to 2pi range
ind=az<0; az(ind)=az(ind)+2*pi;
ind=az>2*pi; az(ind)=az(ind)-2*pi;

% Convert back to deg
maz=az/d2r;
mel=el/d2r;

return
