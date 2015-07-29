function [fm,sa,gof]=do_fit(pars,freepars,obs,ide, xerrs, yerrs)
% [fm,sa,gof]=do_fit(pars,freepars,obs,ide)
%
% Fit a set of points to the model
% Provide the space angle errors as extra out arg

%[fm,fme,gof,stat]=matmin('gof_spaceangle',pars,freepars,'pointing_model',obs.az,obs.el,ide.az,ide.el,xerrs,yerrs);

if(nargin==4)
  xerrs = ones(size(obs.az));
  yerrs = xerrs;
end

ind=isnan(obs.az)|isnan(obs.el);
obs.el=obs.el(~ind);
obs.az=obs.az(~ind);
ide.el=ide.el(~ind);
ide.az=ide.az(~ind);
xerrs=xerrs(~ind);
yerrs=yerrs(~ind);

[fm,fme,gof,stat]=matmin('gof_pt_mod',pars,freepars,'pointing_model',obs.az,obs.el,ide.az,ide.el,xerrs,yerrs);

% Find the model values
[mva.az,mva.el]=pointing_model(fm,ide.az,ide.el);

% Find the space angles between observed and model values
sa=spaceangle(obs.az,obs.el,mva.az,mva.el,'deg');

return
