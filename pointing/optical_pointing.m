function [fm,ide,obs, sa, om]=optical_pointing(d, filename)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     % [fm,ide,obs]=optical_pointing(d, filename)
%  function should give us an optical pointing model.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  
% create some easier to parse variables:

%obs.az=d.antenna0.tracker.horiz_topo(:,1);
%obs.el=d.antenna0.tracker.horiz_topo(:,2);

%ide.az=d.antenna0.tracker.horiz_mount(:,1);
%ide.el=d.antenna0.tracker.horiz_mount(:,2);

ide.az=d.antenna0.tracker.horiz_topo(:,1);
ide.el=d.antenna0.tracker.horiz_topo(:,2);

obs.az=d.antenna0.tracker.horiz_mount(:,1);
obs.el=d.antenna0.tracker.horiz_mount(:,2);

tilt.x=zeros(size(ide.az));
tilt.y=zeros(size(ide.az));

om(:,1)=d.antenna0.tracker.flexure(:,1);
om(:,2)=d.antenna0.tracker.flexure(:,2);
om(:,3)=d.antenna0.tracker.tilts(:,1);
om(:,4)=d.antenna0.tracker.tilts(:,2);
om(:,5)=d.antenna0.tracker.tilts(:,3);
om(:,6)=d.antenna0.tracker.fixedCollimation(:,1);
om(:,7)=d.antenna0.tracker.fixedCollimation(:,2);
om(:,8)=d.antenna0.tracker.encoder_off(:,1);
om(:,9)=d.antenna0.tracker.encoder_off(:,2);

off.az=d.antenna0.tracker.horiz_off(:,1);
off.el=d.antenna0.tracker.horiz_off(:,2);

date=utc2date(d.array.frame.utc(1));
%sources = unique(d.antenna0.tracker.source);


% The actual online values are not recorded - reconstruct them by
% adding in encoder zero points
obs.az=obs.az+om(:,8);
obs.el=obs.el+om(:,9);

% From now on just use online model as it was at start of run -
% shouldn't matter if it changes during a run although there is no
% reason why it ever should.
for i=1:size(om,2)
  if(any(om(:,i)~=om(1,i)))
    warning(sprintf('Online model parameter %d changed during run',i));
  end
end
om=om(1,:);

% Put az values into 0-360 range
ind=obs.az<0; obs.az(ind)=obs.az(ind)+360;
ind=obs.az>360; obs.az(ind)=obs.az(ind)-360;

% Convert tilt meter values to degrees
tilt.x=tilt.x/60; tilt.y=tilt.y/60;

disp(sprintf('Read in %d data points',length(obs.az)));
%disp(sprintf('  ... %d unique sources', length(sources)));


% Fit the model
%fm=singlestep_fit(om,ide,obs,filename);
[fm sa]=stepbystep_fit(om,ide,obs,filename);

% what to add to the files
display(sprintf('encoder_zeros  %1.4f, %1.4f, 0', fm(8), fm(9)));
display(sprintf('tilts  %4.4f, %4.4f, %4.4f', fm(3), fm(4), fm(5)));
display(sprintf('collimate_fixed optical, %4.4f, %4.4f, ptel=0+1+2', fm(6), ...
    fm(7)));
display(sprintf('flexure optical, %4.4f, %4.4f, ptel=0+1+2', fm(1), fm(2)));



return;


