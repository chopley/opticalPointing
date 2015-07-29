function gof=gof_pt_mod(tm,modfunc,oaz,oel,iaz,iel,xerrs,yerrs)
% gof=gof_pt_mod(tm,modfunc,oaz,oel,iaz,iel,xerrs,yerrs)
%
% Goodness of fit function used to compare online oaz,oel points to
% those derived from ideal iaz,iel using modfunc(tm,iaz,iel)

[maz,mel]=feval(modfunc,tm,iaz,iel);

% take model values and observed values, find chi-square for az, el.
% weight by errors, and insert cos(el) depenedence to az
d2r=pi/180;
% Change to radians
%xerrs=xerrs/d2r;yerrs=yerrs/d2r;
%oaz=d2r*oaz;oel=oel*d2r;maz=maz*d2r;mel=mel*d2r;

% Need to change errors from XY to AzEl.  
if(nargin<8)
  az_errs = ones(size(iel));
  el_errs = ones(size(iel));
else
%  az_errs=xerrs./cos(d2r*iel);
%  el_errs=yerrs;
end
az_errs = xerrs;
el_errs = yerrs;

% Weight the offset by the error, weight the az differences by cos(el)
chi_az=(maz-oaz)./az_errs;
chi_el=(mel-oel)./el_errs;
chisq_az=sum(chi_az.^2);
chisq_el=sum(chi_el.^2);

gof=chisq_az+chisq_el;

