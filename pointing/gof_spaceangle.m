function gof=gof_spaceangle(tm,modfunc,oaz,oel,iaz,iel)
% gof=gof_spaceangle(tm,modfunc,oaz,oel,iaz,iel)
%
% Goodness of fit function used to compare online oaz,oel points to
% those derived from ideal iaz,iel using modfunc(tm,iaz,iel)

[maz,mel]=feval(modfunc,tm,iaz,iel);

sa=spaceangle(oaz,oel,maz,mel,'deg');

gof=sum(sa.^2);
