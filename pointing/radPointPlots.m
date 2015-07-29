function [obs, off, ide, out] = radioPointPlots(d, plotparams, type, field, ...
    maindir, useMel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [obs, off, ide] = radioPointPlots(d, plotparams, type, field, maindir)
%
%  pipeline radio pointing analysis.
%
%   d is the data structure
%   type is either 
%   'raster' -- schedule did a full raster on the source (az/el)
%   'scan'  -- schedule did an az/el scan
%   'cross'   -- schedule did discrete steps and integrated
%   gen is whether to generate the plots or not.
%
%   off are the offsets.
%  sjcm
%
%  modified:  3/23/2011 - sjcm:  made it compatible with data from before
%  October 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<5)
  maindir = [];
  useMel  = 0;
end
if(nargin<6)
  useMel = 0;
end



% we must determine the az/el offset between the source and its "thought-of"
% location 
% find data where the tracker utc is not 
if(length(d.antenna0.tracker.utc)~=length(unique(d.antenna0.tracker.utc)))
  d.antenna0.tracker.utc = d.array.frame.utc + 1/60/60/24;
end

numCols = size(d.antenna0.receiver.data,2);

if( (mean(d.antenna0.thermal.coldLoad) > 20) & (useMel == 0) )
  % isn't needed with Mel's stuff.
  % I channels are negative
  I1min = min(d.antenna0.receiver.data(:,1));
  I2min = min(d.antenna0.receiver.data(:,numCols));
  d.antenna0.receiver.data(:,1) =   d.antenna0.receiver.data(:,1) - I1min;
  d.antenna0.receiver.data(:,numCols) =   d.antenna0.receiver.data(:,numCols) - I2min;  
end



switch type
  case 'raster'
%    [obs ide] = calcOffRaster(d, gen, field);
 
  case 'scan'
    if(useMel)
      [obs ide off out] = calcOffScanMel(d, plotparams, field, maindir);          
    else
      [obs ide off out] = calcOffScan(d, plotparams, field, maindir);
    end

    

  case 'cross'
    [obs ide off] = calcOffCross(d, plotparams, field);    
end

% Put az values into 0-360 range
ind = obs.az<0;   obs.az(ind) = obs.az(ind)+360;
ind = obs.az>360; obs.az(ind) = obs.az(ind)-360;

% now we just return the obs, ide, off



return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [obs ide off] = calcOffScan(d, plotparams, field, maindir)
%
%setPlotDisplay(plotparams.plot);
%
%% load up the current model
%om(:,1)=d.antenna0.tracker.flexure(:,1);
%om(:,2)=d.antenna0.tracker.flexure(:,2);
%om(:,3)=d.antenna0.tracker.tilts(:,1);
%om(:,4)=d.antenna0.tracker.tilts(:,2);
%om(:,5)=d.antenna0.tracker.tilts(:,3);
%om(:,6)=d.antenna0.tracker.fixedCollimation(:,1);
%om(:,7)=d.antenna0.tracker.fixedCollimation(:,2);
%om(:,8)=d.antenna0.tracker.encoder_off(:,1);
%om(:,9)=d.antenna0.tracker.encoder_off(:,2);
%om = mean(om);
%
%% I. split them up into respective scans
%% Ia.  find the start and stop index of each scan.
%[s e] = findStartStop(d.index.radio_point_scan.slow);
%
%if(isempty(e) | isempty(s))
%  obs.az = nan;
%  obs.el = nan;
%  ide.az = nan;
%  ide.el = nan;
%  off.az = nan;
%  off.el = nan;
%  goodPoint = [0 0];
%  display('No crosses');
%  s = [];
%  return;
%  
%end
%
%if(e(1)<s(1))
%  obs.az = nan;
%  obs.el = nan;
%  ide.az = nan;
%  ide.el = nan;
%  off.az = nan;
%  off.el = nan;
%  goodPoint = [0 0];
%  display('No crosses');
%  s = [];
%  return;
%end
%
%betaAzAll = nan(length(s), 2, 3);
%betaElAll = nan(length(s), 2, 3);
%
%plotIndex = 0;
%for m=1:length(s)
%  % cut the data for each scan
%  ind = zeros(size(d.array.frame.features));
%  ind(s(m):e(m)) = 1;
%  ind = logical(ind);
%  dc  = framecut(d, ind);
%
%  idevals = dc.antenna0.tracker.horiz_topo(1,:);
%  [obsvals(1) obsvals(2)] = pointing_model(om, idevals(1), idevals(2));
%  ide.az(m) = idevals(1);
%  ide.el(m) = idevals(2);
%  obs.az(m) = obsvals(1);
%  obs.el(m) = obsvals(2);
%  
%  % next we split it into the azimuth/elevation scan
%  % for the elevation scan, we fit a line and remove the slope (effectively
%  % a sky dip), and fit a gaussian to the remainder
%  % for the azimuth scan, we just fit a gaussian to it.
%
%  % get the offsets from the control system.
%  azApp = interp1(dc.antenna0.tracker.utc, ...
%      dc.antenna0.tracker.horiz_topo(:,1), dc.antenna0.receiver.utc);
%  azOffSave = azApp - dc.antenna0.servo.apparent(:,1);
%  azOffSave = -azOffSave;
%  
%  elApp = interp1(dc.antenna0.tracker.utc, ...
%      dc.antenna0.tracker.horiz_topo(:,2), dc.antenna0.receiver.utc);
%  elOffSave = elApp - dc.antenna0.servo.apparent(:,2);
%  elOffSave = -elOffSave;
%
%  % assuming the slowest we would ever scan is at 0.2 degrees per second, a
%  % az scan happens when the az speed is mostly constant, increasing, and
%  % greater than 0.2 degrees/s, which is 
%  if(dc.array.frame.utc(1)<date2mjd(2010,10,1,0,0,0))
%    azInd = (deriv(dc.antenna0.servo.az))<-0.2/100;
%  else
%    azInd = (deriv(dc.antenna0.servo.az))>0.1/100;
%  end
%  azIs  = dc.antenna0.receiver.data(azInd,[1 6]);
%  azIs(:,1) = smooth(azIs(:,1));
%  azIs(:,2) = smooth(azIs(:,2));
%
%  % if amplitude is negative, no peak in data
%  
%  % same goes for elevation
%  if(dc.array.frame.utc(1)<date2mjd(2010,10,1,0,0,0))
%    elInd = deriv(dc.antenna0.servo.el)<-0.2/100;
%  else
%    elInd = deriv(dc.antenna0.servo.el)>0.1/100;
%  end
%  elIs  = dc.antenna0.receiver.data(elInd,[1 6]);
%  elIs(:,1) = smooth(elIs(:,1));
%  elIs(:,2) = smooth(elIs(:,2));
%  
%  
%  testVals = [120:30:800];
%  for mm=1:2
%    
%    if(~isempty(find(azInd)))
%      % fit a gaussian about the max in both cases
%      azOff = azOffSave;
%      azI = azIs(:,mm);
%      f = find(azI == max(azI));
%      azOff = azOff(azInd);
%      
%      stopAz = 0;
%      
%      azIndex = 1;
%      while(stopAz==0)
%	azFit = find(azInd);
%	if(f<=testVals(azIndex) | (f+testVals(azIndex))>=length(azI))
%	  badAzPoint = 1;
%	  azFit = [];
%	  xoff = nan;
%	  azIndex = azIndex+1;
%	  betaAz(1:3) = nan;
%	else
%	  azFit = azOff((f-testVals(azIndex):f+testVals(azIndex)));
%	  azIFit = azI(f-testVals(azIndex):f+testVals(azIndex));
%	  
%%	  display(sprintf('testing value: %d', testVals(azIndex)));
%	  
%	  OPTIONS = optimset('Display', 'off', 'MaxIter', 200000, ...
%	      'MaxFunEvals', 20000, 'TolX', 0.00005, 'TolFun', 0.00005);
%	  LB = [0 azFit(testVals(azIndex))-0.2  -10 -inf -inf];
%	  UB = [max(azIFit) azFit(testVals(azIndex))+0.2 10 inf inf];
%	  x0 = [1 azFit(testVals(azIndex)-1) 1 1 1];
%	  [betaAz resnorm residual exitflag] = lsqnonlin(@(x) mchisq2(azIFit, gaussfit(x, azFit)), x0, LB, UB, OPTIONS);
%	  
%	  amplitude = betaAz(1);
%	  xoff      = betaAz(2);
%	  sigma     = betaAz(3);
%	  
%	  % other terms remove a baseline
%	  if (amplitude<0 | abs(xoff)>1 | resnorm>500)
%	    badAzPoint = 1;
%	  else
%	    badAzPoint = 0;
%	  end 
%	  azFit = gaussfit(betaAz, azOff);
%	  
%	  if(abs(xoff)<2)
%	    stopAz = 1;
%	  else
%	    azIndex = azIndex+1;
%	  end
%	end
%	
%	if(azIndex>length(testVals))
%	  stopAz = 1;
%	end
%      end
%      
%      betaAzAll(m,mm,:) = betaAz(1:3);
%    end
%    
%    
%    if(~isempty(elInd))
%      elOff = elOffSave;
%      elI = elIs(:,mm);
%      % first we remove the offset from slewing the atmosphere
%      % from the sky dip code
%      x = 1./(sind(dc.antenna0.servo.el(elInd)));
%      y = abs(elI);
%      [tau Tground] = linfit(x,y);
%      % remove the effect
%      elI = elI - Tground(1) - tau(1)./sind(dc.antenna0.servo.el(elInd));
%      
%      % fit a gaussian about the max
%      f = find(elI == max(elI));
%      elOff2 = elOff(elInd);
%      
%      stopEl = 0;
%      elIndex = 1;
%      while(stopEl==0)
%	elFit = find(elInd);
%	if(f<=testVals(elIndex) | (f+testVals(elIndex))>=length(elI))
%	  badElPoint = 1;
%	  elFit = [];
%	  yoff = nan;
%	  elIndex = elIndex+1;
%	  betaEl = [nan nan nan];
%	else
%	  elFit = elOff2((f-testVals(elIndex):f+testVals(elIndex)));
%	  elIFit = elI(f-testVals(elIndex):f+testVals(elIndex));
%	  
%%	  display(sprintf('el testing value: %d', testVals(elIndex)));
%	  [betaEl] = nlinfit(elFit, elIFit, @gaussfit, [1 elFit(testVals(elIndex)) 1 1 1]);
%	  amplitude = betaEl(1);
%	  yoff      = betaEl(2);
%	  sigma     = betaEl(3);
%	  % other terms remove a baseline
%	  if (amplitude<0 | abs(yoff)>1)
%	    badElPoint = 1;
%	  else
%	    badElPoint = 0;
%	  end 
%	  elFit = gaussfit(betaEl, elOff2);
%	  
%	  if(abs(yoff)<1)
%	    stopEl = 1;
%	  else
%	    elIndex = elIndex+1;
%	  end
%	end
%	if(elIndex>length(testVals))
%	  stopEl = 1;
%	end
%	
%      end
%      betaElAll(m,mm,:) = betaEl(1:3);      
%    end
%    
%%    display('here');
%    badPoint(m,mm) = badAzPoint | badElPoint;
%    
%    if(~badPoint(m,mm))
%      % then we plot it
%
%      if(~any(isnan(betaAzAll(m,mm,:))))
%	
%	subplot(2,1,1)
%	plot((azOff)*cos(dc.antenna0.servo.el(1)*pi/180), azI');
%	xlabel('Az Offset (degrees)');
%	ylabel('Intensity');
%	if(~isempty(azFit))
%	  hold on; plot(azOff*cos(dc.antenna0.servo.el(1)*pi/180), azFit, 'r'); hold off
%	  axis([min(azOff) max(azOff) min(azI) max(azI)]);
%	end
%	
%      end
%      
%      if(~any(isnan(betaElAll(m,mm,:))))
%	subplot(2,1,2)
%	plot(elOff(elInd), elI);
%	xlabel('El Offset (degrees)');
%	ylabel('Intensity');
%	if(~isempty(elFit))
%	  hold on; plot(elOff(elInd), elFit, 'r'); hold off
%	  axis([min(elOff) max(elOff) min(elI) max(elI)]);
%	end
%      end
%      eval(sprintf('gtitle(''Pointing scan on %s - Chan %d'');', ...
%	  dc.antenna0.tracker.source{1}, mm));
%      
%      display(sprintf('Source Name: %s', dc.antenna0.tracker.source{1}));
%	
%      if(plotparams.save==1)
%	plotIndex = plotIndex+1;
%	if(isempty(maindir))
%	  maindir = getMainDir(d, field);
%	end
%	dbclear if error
%	set(gcf,'paperposition',[0 0 6.0 9.0])
%	filename = sprintf('%s/pointing/fig%d.png', maindir, plotIndex);
%	eval(sprintf('print -dpng -r70 %s', filename));
%	dbstop if error
%      end
%    end
%  end
%  
%  clear azIs;
%  clear elIs;
%  pause(0.1);
%  sources{m} = unique(dc.antenna0.tracker.source);
%  timeVal(m) = dc.antenna0.receiver.utc(1);
%end
%
%badPoint = logical(badPoint);
%elOff = betaElAll(:,:,2);
%azOff = betaAzAll(:,:,2);
%elOff(badPoint) = nan;
%azOff(badPoint) = nan;
%
%
%% if hte values for both channels don't agree, toss them.
%dEl = abs(diff(elOff, [], 2));
%dAz = abs(diff(azOff, [], 2).*cos(obs.el*pi/180)');
%
%% if they differ by more than 20 arcminutes, throw them out
%f = find(dEl>20/60 | dAz>20/60);
%elOff(f,:) = nan;
%azOff(f,:) = nan;
%
%off.az = nanmean(azOff, 2);
%off.el = nanmean(elOff, 2);
%off.azErr = nanstd(azOff, [], 2);
%off.elErr = nanstd(elOff, [], 2);
%% if there's only one sample, we need to set a weight equal to the other
%% ones.
%indSingle = off.azErr == 0;
%off.azErr(indSingle) = median(off.azErr(~indSingle));
%off.elErr(indSingle) = median(off.elErr(~indSingle));
%
%display('Done with all crosses');
%
%% add final offests to the observed az/el
%obs.az = obs.az + off.az';
%obs.el = obs.el + off.el';
%
%
%% throw out the data that is bad
%ind = isnan(obs.az);
%obs.az(ind) = [];
%obs.el(ind) = [];
%ide.az(ind) = [];
%ide.el(ind) = [];
%off.az(ind) = [];
%off.el(ind) = [];
%off.azErr(ind) = [];
%off.elErr(ind) = [];
%off.sources = sources(~ind);
%off.timeVal = timeVal(~ind);
%
%
%return;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obs ide] = calcOffRaster(d, gen)

% I. split them up into respective scans
% Ia.  find the start and stop index of each scan.
[s e] = findStartStop(d.array.frame.features);


for m=1:length(s)
  % cut the data for each scan
  ind = zeros(size(d.array.frame.features));
  ind(s(m):e(m)) = 1;
  ind = logical(ind);
  dc  = framecut(d, ind);

  % calculate the offset from nominal for this scan
  [obs.az(m) obs.el(m)] = sourceScanMap(dc, 0.5, 0);
  
  % calculate the az/el for the source
  ide.az(m) = mean(dc.antenna0.tracker.horiz_mount(:,1));
  ide.el(m) = mean(dc.antenna0.tracker.horiz_mount(:,2));
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [s e] = findStartStop(features)
%
%% start indices are when the difference in the feature value increases, stop
%% is when it decreases
%s = find(diff(double(features))>0);
%
%% end is when it decreases
%e = find(diff(double(features))<0);
%
%if(length(e)<length(s))
%  % the dimensions don't match
%  if(s(1)<e(1))
%    % if start is less than end, that means we're usually missing hte last
%    % endpoint
%    e(length(s)) = length(features);
%  else
%    s(1) = [];
%  end
%elseif(length(s)<length(e))
%  s(1) = [];
%end
%s = s+1;
%
%goodPoint = ones(size(s));
%for m=1:length(s)
%  thisFeat = unique(features(s(m):e(m)));
%  if(thisFeat ~= 1)
%    goodPoint(m) = 0;
%  end
%end  
%
%s = s(logical(goodPoint));
%e = e(logical(goodPoint));
%
%return;


function x = mchisq2(dataVals, modelVals)

x = dataVals - modelVals;

f = find(isnan(x));
x(f) = [];
f = find(isinf(x));
x(f) = [];

if(isempty(x))
  display('damn it all to hell')
  keyboard;
end

return;

	    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obs ide off] = calcOffCross(d, plotparams, field);    

setPlotDisplay(plotparams.plot);

% load up the current model
om(:,1)=d.antenna0.tracker.flexure(:,1);
om(:,2)=d.antenna0.tracker.flexure(:,2);
om(:,3)=d.antenna0.tracker.tilts(:,1);
om(:,4)=d.antenna0.tracker.tilts(:,2);
om(:,5)=d.antenna0.tracker.tilts(:,3);
om(:,6)=d.antenna0.tracker.fixedCollimation(:,1);
om(:,7)=d.antenna0.tracker.fixedCollimation(:,2);
om(:,8)=d.antenna0.tracker.encoder_off(:,1);
om(:,9)=d.antenna0.tracker.encoder_off(:,2);
om = mean(om);

% calculates the offsets when doing a pointing cross:  fits for a gaussian
% in x and y offsets, converts to offsets in az/el, and calculates the
% observed location of the source
% should be very similar to the sza

[s e] = findCrossStartStop(d);

% now for each of these fields, we split them up into the x and y offsets,
% fit a gaussian, and find the offsets.  convert offsets to az/el offsets,
% and be done.
for m=1:length(s)
  % cut the data for each scan
  ind = zeros(size(d.array.frame.features));
  ind(s(m):e(m)) = 1;
  ind = logical(ind);
  dc  = framecut(d, ind);
  
  [intStart intStop] = findCrossStepStartStop(dc);
  if(length(intStart)~=26)
    display('You may have some mixed up data.  not the right number of points');
    obs = [];
    ide = [];
    off = [];
    return;
  end
  
  for n=1:length(intStart)
    ind = zeros(size(dc.array.frame.features));
    ind(intStart(n):intStop(n)) = 1;
    ind = logical(ind);
    dcc = framecut(dc, ind);
    
    offsets(n,:) = median(dcc.antenna0.tracker.sky_xy_off);
    meanIs(n,:) = mean(dcc.antenna0.receiver.data(:,[1 6]));
    stdIs(n,:)  = std(dcc.antenna0.receiver.data(:,[1 6]));
    
    ideal(n,:) = dcc.antenna0.tracker.horiz_topo(2,:);
    observed(n,:) = dcc.antenna0.tracker.horiz_mount(2,:) + ...
	dcc.antenna0.tracker.encoder_off(2,:);
  end
  fy = find(offsets(:,2)~=0);
  yIndex = fy(1):last(fy);
  fx = find(offsets(:,1)~=0);
  xIndex = fx(1):last(fx);
  
  % fit a gaussian to both data sets.
  xOffVals = offsets(xIndex,1);
  yOffVals = offsets(yIndex,2);
  xIntVals = abs(meanIs(xIndex,:));
  yIntVals = abs(meanIs(yIndex,:));
  xErrVals = stdIs(xIndex,:);
  yErrVals = stdIs(yIndex,:);

  % find the x offset
  x0 = [max(xIntVals(:,1)) 0 0.2];
  ind = abs(xOffVals)<=0.5;
  [px(1,:) pex(1,:)] = matmin('chisq', x0, [], 'gauss', ...
      xIntVals(ind,1), xErrVals(ind,1), xOffVals(ind));
  [px(2,:) pex(2,:)] = matmin('chisq', x0, [], 'gauss', ...
      xIntVals(ind,2), xErrVals(ind,2), xOffVals(ind));
  fullX = min(xOffVals):0.01:max(xOffVals);
  aa(:,1) = gauss(px(1,:), fullX);
  aa(:,2) = gauss(px(2,:), fullX);
  offset.x(m) = vwstat(px(:,2), pex(:,2).^2);

  elFit = xOffVals(ind);
  elIFit = xIntVals(ind,2);

  [betaEl] = nlinfit(elFit, elIFit, @gaussfit, [1 0 1 1 1]);
  amplitude = betaEl(1);
  yoff      = betaEl(2);
  sigma     = betaEl(3);
  xoffvals  = min(xOffVals):0.01:max(xOffVals);
  elFitFin = gaussfit(betaEl, xoffvals);
  
  y0 = [max(yIntVals(:,1)) 0 0.2];
  ind = abs(yOffVals)<=0.5;
  [py(1,:) pey(1,:)] = matmin('chisq', y0, [], 'gauss', ...
      yIntVals(ind,1), yErrVals(ind,1), yOffVals(ind));
  [py(2,:) pey(2,:)] = matmin('chisq', y0, [], 'gauss', ...
      yIntVals(ind,2), yErrVals(ind,2), yOffVals(ind));  
  fullY = min(yOffVals):0.01:max(yOffVals);
  bb(:,1) = gauss(py(1,:), fullY);
  bb(:,2) = gauss(py(2,:), fullY);  
  offset.y(m) = vwstat(py(:,2), pey(:,2).^2);

  if(any(abs(px(:,2)) > max(xOffVals)))
    display('Suspect cross');
    badXPoint = 1;
  else
    badXPoint = 0;
  end
  if(abs(diff(px(:,2)))>0.05)
    display('Best Fit X offset does not agree between channels');
    badXPoint = 1;
  end
  

  if(any(abs(py(:,2)) > max(xOffVals)))
    display('Suspect cross');
    badYPoint = 1;
  else
    badYPoint = 0;
  end
  if(abs(diff(py(:,2)))>0.05)
    display('Best Fit Y offset does not agree between channels');
    badYPoint = 1;
  end
  
  
  
  % plot the results to see the fit.
  meanX = mean(xIntVals);
  xIntVals(:,1) = xIntVals(:,1) - meanX(1);
  xIntVals(:,2) = xIntVals(:,2) - meanX(2);
  aa(:,1) = aa(:,1) - meanX(1);
  aa(:,2) = aa(:,2) - meanX(2);
  
  meanY = mean(yIntVals);
  yIntVals(:,1) = yIntVals(:,1) - meanY(1);
  yIntVals(:,2) = yIntVals(:,2) - meanY(2);
  bb(:,1) = bb(:,1) - meanY(1);
  bb(:,2) = bb(:,2) - meanY(2);

  subplot(2,1,1)
  plot(xOffVals, xIntVals, '*-');
  hold on
  plot(fullX, aa, 'r')
  axis([min(xOffVals) max(xOffVals) min(xIntVals(:)) max(xIntVals(:))]);
  legend('Channel 1', 'Channel 2', 'Fit Gaussian');
  xlabel('X offset');
  ylabel('Counts');
  hold off;
  
  subplot(2,1,2)
  plot(yOffVals, yIntVals, '*-');
  hold on
  plot(fullY, bb, 'r')
  axis([min(yOffVals) max(yOffVals) min(yIntVals(:)) max(yIntVals(:))]);
  legend('Channel 1', 'Channel 2', 'Fit Gaussian');
  xlabel('Y offset');
  ylabel('Counts');
  hold off;
  
  if(plotflag.save==1)
    plotIndex = plotIndex+1;
    maindir = getMainDir(d, field);
    dbclear if error
    set(gcf,'paperposition',[0 0 6.0 9.0])
    filename = sprintf('%s/pointing/fig%d.png', maindir, plotIndex);
    eval(sprintf('print -dpng -r70 %s', filename));
    dbstop if error
  end
    
  display(sprintf('Source Name: %s', dc.antenna0.tracker.source{1}));
  eval(sprintf('display(''X offset: %3.4f, Y offset: %3.4f'');', ...
      offset.x(m), offset.y(m)));
  if(badXPoint | badYPoint)
    badPoint = 1;
    display('Cross Suspect');
  end
  
  
  % check if we want to keep this source
%  r = input('Keep This point? [y/n]: ', 's');
%  
%  if(strcmp(r, 'y') || strcmp(r, 'Y'))
%    badPoint = 0;
%  elseif(strcmp(r, 'n') ||  strcmp(r, 'N'))
%    badPoint = 1;
%  elseif(strcmp(r, 'k'))
%    keyboard;
%  end  
  
  if(~badPoint)
    % now for the ideal and observed
    f = find(offsets(:,1)==0 & offsets(:,2)==0);
    
    ide.az(m) = ideal(f(1),1);
    ide.el(m) = ideal(f(1),2);
    obs.az(m) = observed(f(1),1);
    obs.el(m) = observed(f(1),2);
    
    % convert great circle offsets into az/el offsets
    
    a.az = obs.az(m);
    a.el = obs.el(m);
    b = xyoff_to_azel(a, offset.x(m), offset.y(m));
    off.x(m) = offset.x(m);
    off.y(m) = offset.y(m);
    obs.az(m) = b.az;
    obs.el(m) = b.el;
  
  else
    ide.az(m) = nan;
    ide.el(m) = nan;
    obs.az(m) = nan;
    obs.el(m) = nan;
    off.x(m) = offset.x(m);
    off.y(m) = offset.y(m);
  end
end  


return;

function [s e] = findCrossStartStop(d)
% check for when the radio pointing flag switches.

index = d.index.radio_point_cross.slow;

s = find(diff(index)>0)+1;
e  = find(diff(index)<0);

% that's it.
if(length(e)<length(s))
  % the dimensions don't match
  if(s(1)<e(1))
    % if start is less than end, that means we're usually missing hte last
    % endpoint
    %    e(length(s)) = length(d.array.frame.features);
    % cut that last start out.
    s(length(s)) = [];
  else
    s(1) = [];
  end
elseif(length(s)<length(e))
  s(1) = [];
end

return;


function [s e] = findCrossStepStartStop(d)
% check for when the radio pointing flag switches.

index = d.index.source.slow;

s = find(diff(index)>0)+1;
e  = find(diff(index)<0);

% that's it.
if(length(e)<length(s))
  % the dimensions don't match
  if(s(1)<e(1))
    % if start is less than end, that means we're usually missing hte last
    % endpoint
    e(length(s)) = length(d.array.frame.features);
  else
    s(1) = [];
  end
elseif(length(s)<length(e))
  s(1) = [];
end

return;


%%%%%%%%%%%%%%%
function azel=xyoff_to_azel(azel,x,y)
% Apply xy sky offset to az,el coord
% following Martin Shepherd in tracker.c
% Convert all to radians
d2r=pi/180;
azel.az=azel.az*d2r; azel.el=azel.el*d2r;
x=x*d2r; y=y*d2r;

% Convert x,y to r,theta on surface of sphere
t=atan2(sin(x),tan(y));
r=acos(cos(x).*cos(y));

% Apply offset
e=azel.el;
top=sin(t).*sin(r);
bot=cos(e).*cos(r)-sin(e).*sin(r).*cos(t);
azel.az=azel.az+atan2(top,bot);
azel.el=asin(sin(e).*cos(r)+cos(e).*sin(r).*cos(t));

% Back to deg
azel.az=azel.az/d2r; azel.el=azel.el/d2r;

return
