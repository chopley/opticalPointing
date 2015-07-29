function d = apparentAzElSouth(d, includeError)

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
 
if(nargin<2)
  includeError = 0;
end


%get the antenna coordinates that we start with
azprime=d.antenna0.servo.az;
elprime = d.antenna0.servo.el;
maxitr=20;

%lets just see what sort of correction we get at this antenna angle (AZ')- 
%if we assume the difference between AZ (apparent) and AZ'(antenna) is
%small then we can use this number as a start approximation- here we have 3
%coordinates Az (apparent AZ) Az' (antenna Az) AZ'' (meaningless azimuth   
%position calculated by applying a pointing and refraction correction to ...
%    Az') and we
%make the assumptio that Az-Az' is approximately Az'-Az''
[azdoubleprime,eldoubleprime]=calcAzEl2South(azprime,elprime,d);

%now we calculate what sort of offset we have, daz and del-note that this
%offset is between AZ' and AZ''- ie daz=AZ'-AZ''- therefore to go back to
%AZ we use the equation daz~=AZ-AZ' (using approximation) and therefore ...
%    AZ=AZ'+daz

daz=azprime-azdoubleprime;
del=elprime-eldoubleprime;

%Now we assume that the difference between the apparent and antenna
%coordinates is small i.e that the same sort of correction would have been
%found if we'd done the calculation on the apparent coordinates- defining  
%apparent coordinates as AZ and antenna coordinates as AZ' we have
%daz = AZ-AZ' (not the difference and similarity (!!!) with the calculation
%above                                                          
% AZ=AZ'+daz

azapp = azprime+daz;
elapp = elprime+del;

%So these are our first guesses for the apparent azimuth and elevation
%values- they should be ok but lets just check by going from these values
%back to the antenna values and comparing

[azprime2,elprime2]=calcAzEl2South(azapp,elapp,d);

errazprime = (azprime-azprime2);
errelprime = (elprime-elprime2);

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

alen=length(daz);
i=0;
%while(((maxazerr >0.0001) || (maxelerr>0.0001)) && (i<maxitr))
 % i=i+1;
  %disp(['Iteration (max): ', num2str(i),'(',num2str(maxitr),')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ogk edit 1: change while loop conditional statement to
% mod(maxazerr,360)<...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(((mod(maxazerr,360) >0.0001) || (maxelerr>0.0001)) && (i<maxitr))
  i=i+1;
  disp(['Iteration (max): ', num2str(i),'(',num2str(maxitr),')']);  
  [azprime2,elprime2]=calcAzEl2South(azapp,elapp,d);
  errazprime = wrap180(azprime - azprime2);
  errelprime = wrap180(elprime - elprime2);
  a=isnan(errazprime);
  errazprime(a)=0;
  a=isnan(errelprime);
  errelprime(a)=0;
  
  maxazerr=max(abs(errazprime));
  maxelerr=max(abs(errelprime));
  daz = daz + errazprime;
  del = del + errelprime;
  azapp= wrap360(azprime+daz);
  elapp = wrap360(elprime+del);
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

erraz=errazprime;
errel=errelprime;
aza=azapp;
ela=elapp;

if(includeError)
  d.antenna0.servo.apparent = [aza ela erraz errel];
else
  d.antenna0.servo.apparent = [aza ela]; 
end

return;
