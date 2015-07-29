function [az, el] = applyDiurnalAberration(az, el)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [az, el] = applyDiurnalAberration(az, el)
%
%
%  name is pretty self-explanatory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tiny = 1e-7;

theta = acos(cos(el).*sin(az));

% If the source is along the direction of the tangential velocity    
%  vector then the azimuth and elevation are unaffected.        
if(abs(theta)<tiny)
  return;
end

x = -cos(el).*cos(az);
if(x~=0 | sin(el) ~= 0)
  psi = atan2(x, sin(el));
else
  psi = 0;
end

% The aberration causes the source to appear to shift towards            
% the direction of the site's velocity. Since the tangential   
% velocity of the Earth is directed along the y axis, this     
% results in a reduction of theta.                               
	
theta -= actual_.velocity/cvel.*sin(theta);

% Compute the modified azimuth and elevation.                        
el = asin(sin(theta).*cos(psi));

% Compute the new azimuth, unless pointing at the zenith. The          
% azimuth is irrelevant at the zenith.                          
if(cos(theta) ~= 0 | (sin(theta).*sin(psi))~= 0)
  az = atan2(cos(theta), -sin(theta).*sin(psi));
else
  az = 0;
end

return;

