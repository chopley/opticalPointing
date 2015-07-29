function d = apparentAzEl(d, includeError)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function d = apparentAzEl(d, includeError)
%
%  functions that retracts the pointing model and refraction corrections to
%  go from encoder az/el to sky az/el.
%
%  d - data structure
%  
%  creates new subfield:
%    d.antenna0.servo.apparent, which is [Nby4], where each column
%    corresponds to :
%    [apparentAzimuth, apparentElevation, error in fit for Az, error in fit
%    for El];
%
%
%aza is the apparent azimuth
%ela is the apparent elevation
%erraz is the error between the antenna coordinate recorder in the data and
%what we get when we use the apparent position to start with for refraction
%and pointing calcualtions
%similarly for erralt
%
%  CJC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ogk edit 1:
% Changed the condition in the while loop to
% Azimuth Error > threshold -> mod(Azimuth Error,360)>threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  sjcm Aug 10, 2011:  modified inputs so we don't report back the error in
%  az/el automatically (which i don't think actually gets used anywhere);
%  memory problem.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Paddy Apr 11 2012: used wrap180 consistently throughout, avoids
%  need for mod in threshold condition. Tidied up a bit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<2)
  includeError = 0;
end

%get the antenna coordinates that we start with
azprime=d.antenna0.servo.az;
elprime = d.antenna0.servo.el;
maxitr=20;
% if flag mask has been filled in, remove data with bits 0, 3 10 set
% (corresponds to frame not received, tracker lacking status, and possible
% servo timing issues).
if(issubfield(d, 'flags', 'mask','fast'))
    bitmask = [0 3 10];
    indgood  = ~bitand(d.flags.mask.fast, sum(2.^bitmask)) & ...
	d.antenna0.servo.el>5 & d.antenna0.servo.el<88;
    dc = framecut(d, indgood);
else
  indgood = true(size(azprime));
  dc = d;
end


%lets just see what sort of correction we get at this antenna angle (AZ')- 
%if we assume the difference between AZ (apparent) and AZ'(antenna) is
%small then we can use this number as a start approximation- here we have 3
%coordinates Az (apparent AZ) Az' (antenna Az) AZ'' (meaningless azimuth   
%position calculated by applying a pointing and refraction correction to ...
%    Az') and we
%make the assumptio that Az-Az' is approximately Az'-Az''
[azdoubleprime,eldoubleprime]=calcAzEl2(azprime(indgood),elprime(indgood),dc);

%now we calculate what sort of offset we have, daz and del-note that this
%offset is between AZ' and AZ''- ie daz=AZ'-AZ''- therefore to go back to
%AZ we use the equation daz~=AZ-AZ' (using approximation) and therefore ...
%    AZ=AZ'+daz

daz=azprime(indgood)-azdoubleprime;
del=elprime(indgood)-eldoubleprime;

%Now we assume that the difference between the apparent and antenna
%coordinates is small i.e that the same sort of correction would have been
%found if we'd done the calculation on the apparent coordinates- defining  
%apparent coordinates as AZ and antenna coordinates as AZ' we have
%daz = AZ-AZ' (not the difference and similarity (!!!) with the calculation
%above                                                          
% AZ=AZ'+daz

azapp = azprime(indgood)+daz;
elapp = elprime(indgood)+del;

%So these are our first guesses for the apparent azimuth and elevation
%values- they should be ok but lets just check by going from these values
%back to the antenna values and comparing

[azprime2,elprime2]=calcAzEl2(azapp,elapp,dc);

errazprime = wrap180(azprime(indgood)-azprime2);
errelprime = elprime(indgood)-elprime2;

maxazerr=max(abs(errazprime));
maxelerr=max(abs(errelprime));
vecazerr=[];
vecelerr=[];
a=isnan(daz);
daz(a)=0;
a=isnan(del);
del(a)=0;
%Now we iteratively adjust the Az (apparent coordinate) until it produces
%the 'correct' AZ

%alen=length(daz);
i=0;

while(((maxazerr >0.0001) || (maxelerr>0.0001)) && (i<maxitr))
  i=i+1;
  disp(['Iteration (max): ', num2str(i),'(',num2str(maxitr),')']);  
  [azprime2,elprime2]=calcAzEl2(azapp,elapp,dc);
  errazprime = wrap180(azprime(indgood) - azprime2);
  errelprime = elprime(indgood) - elprime2;
  a=isnan(errazprime);
  errazprime(a)=0;
  a=isnan(errelprime);
  errelprime(a)=0;
  
  maxazerr=max(abs(errazprime));
  maxelerr=max(abs(errelprime));
  daz = daz + errazprime;
  del = del + errelprime;
  azapp = wrap360(azprime(indgood)+daz);
  elapp = wrap360(elprime(indgood)+del);
  vecazerr=[vecazerr maxazerr];
  vecelerr=[vecelerr maxelerr];
%end
end

if i>=maxitr
	%Error handling in case the removal of the pointing correction isn't
	% working
	error('Apparent Az/El not converging')
end
%might need some error checking here but otherwise

%erraz=errazprime;
%errel=errelprime;
%aza=azapp;
%ela=elapp;
totalIndices = 1:1:length(azprime);

erraz = interp1(totalIndices(indgood), errazprime, totalIndices, ...
    'linear', 'extrap');
errel = interp1(totalIndices(indgood), errelprime, totalIndices, ...
    'linear', 'extrap');
aza   = interp1(totalIndices(indgood), azapp, totalIndices, ...
    'linear', 'extrap');
ela   = interp1(totalIndices(indgood), elapp, totalIndices, ...
    'linear', 'extrap');

if(includeError)
  d.antenna0.servo.apparent = [aza' ela' erraz' errel'];
else
  d.antenna0.servo.apparent = [aza' ela']; 
end

return;
