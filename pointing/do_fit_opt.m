function [fm,sa,gof]=do_fit_opt(pars,freepars,obs,ide)
% [fm,sa,gof]=do_fit(pars,freepars,obs,ide)
%
% Fit a set of points to the model
% Provide the space angle errors as extra out arg
freepars = double(freepars);
ind=isnan(obs.az)|isnan(obs.el);
obs.el=obs.el(~ind);
obs.az=obs.az(~ind);
ide.el=ide.el(~ind);
ide.az=ide.az(~ind);

[fm,fme,gof,stat]=matmin('gof_spaceangle',pars,freepars,'pointing_model',obs.az,obs.el,ide.az,ide.el);

% Find the model values
[mva.az,mva.el]=pointing_model(fm,ide.az,ide.el);


% Find the space angles between observed and model values
sa=spaceangle(obs.az,obs.el,mva.az,mva.el,'deg');

return
