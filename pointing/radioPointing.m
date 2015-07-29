function par = radioPointing(d, type, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function par = radioPointing(d, type, filename)
%
%  radio pointing analysis schedule.
%
%   d is the data structure
%   type is either 
%   'raster' -- schedule did a full raster on the source (az/el)
%   'scan'  -- schedule did an az/el scan
%   'cross'   -- schedule did discrete steps and integrated
%   filename -- string that will be used to name your plots.
%
%  sjcm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we must determine the az/el offset between the source and its "thought-of"
% location 

% find data where the tracker utc is not 
if(length(d.antenna0.tracker.utc)~=length(unique(d.antenna0.tracker.utc)))
  d.antenna0.tracker.utc = d.array.frame.utc + 1/60/60/24;
end


if(d.array.frame.utc(1) < date2mjd(2010, 6, 15))
  d.antenna0.receiver.data = abs(d.antenna0.receiver.data);
end


switch type
  case 'raster'
    [obs ide] = calcOffRaster(d);
 
  case 'scan'
      
    [obs ide off] = calcOffScan(d);
%   [obs1 ide1 off1] = calcOffScanMel(d);
    
  case 'cross'
    [obs ide off] = calcOffCross(d);    
end

% Put az values into 0-360 range
ind = obs.az<0;   obs.az(ind) = obs.az(ind)+360;
ind = obs.az>360; obs.az(ind) = obs.az(ind)-360;

sources = unique(d.antenna0.tracker.source);
disp(sprintf('Read in %d data points',length(obs.az)));
disp(sprintf('  ... %d unique sources', length(sources)));
    
% Fit the model
% we only want to fit the collimation and flexure
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

off.azErr(isnan(off.azErr)) = max(off.azErr);
off.elErr(isnan(off.elErr)) = max(off.elErr);

om = om(1,:);

figure(1)
[rmssao ind sa0] = plot_res(om, ide, obs, 'No Fit', 0);
sa0 = mean(sa0);


if(rmssao*60 < 6)
  display('**********ATTENTION*************');
  display('The RMS of the current pointing is very low');
  display('The model might be good with no update');
else
  display('**********ATTENTION*************');
  display('The RMS of the current pointing is not very low');
  display('Will try to fit for the new solution');
end

old_m = om;
%opt_m = [ -0.0500 -0.0438 -0.0243 -0.0291 -0.1331 0.1142 0 -0.6572 -0.1666];
%opt_m = [ -0.0479 -0.0377 -0.0303 -0.0302 -0.0716 0.1943 0 -0.7072  -0.1594];
%opt_m = [-0.0905   -0.0627   -0.0238   -0.0292   -0.0736    0.1910 0   -0.7956   -0.2366];

%om = opt_m;


figure(2)
% only fit for collimation
fitPars = ones(size(om));
fitPars([1:5 8:9]) = 0;
[par sa, gof] = do_fit(om, fitPars, obs, ide, off.azErr', off.elErr');
plot_res(par, ide, obs, 'Collimation Fit Only', 0)

figure(3)
% only fit for flexure
fitPars = ones(size(om));
fitPars(3:9) = 0;
[par sa, gof] = do_fit(om, fitPars, obs, ide, off.azErr', off.elErr');
plot_res(par, ide, obs, 'Flexure Fit Only', 0)

figure(4)
% fit both collimation and flexure
fitPars = ones(size(om));
fitPars([3:5 8:9]) = 0;
[par sa, gof] = do_fit(om, fitPars, obs, ide, off.azErr', off.elErr');
[par sa, gof] = do_fit(om, fitPars, obs, ide, ones(size(off.azErr')), ones(size(off.elErr')));
[rmssaf ind saf] = plot_res(par, ide, obs, 'Collimation and Flexure', 0);
%saf = mean(saf);

par
figure(5)
[rmssaf ind saf] = plot_res(par, ide, obs, 'Feb Full Solution', 0);

fm = par;

if(rmssaf*60 > 10)
  display('**********WARNING*************');
  display('Fit to data not so great');
  display('Suggest scheduling an optical pointing run');
end

if(rmssaf > rmssao)
  display('**********WARNING*************');
  display('Old model better than new fit');
end


if(any(abs(par)>7))
  display('**********WARNING*************');
  display('Some pointing model terms fit are quite large');
  display('Be suspect of the resulting model');
  display('i.e., CHECK the solution before committing to the control system');
end


if(any(abs(om-fm) > 1.5))
  display('**********WARNING*************');
  display('Some pointing model terms changed by over 1.5 degree');
  display('Be suspect of the resulting model');
  display('i.e., CHECK the solution before committing to the control system');  
end



display('Pointing solution derived');
display('Add these lines to bottom of pointing.init');
display(sprintf('encoder_zeros %3.4f, %3.4f, 0', fm(8:9)));
display(sprintf('tilts %3.4f, %3.4f, %3.4f', fm(3:5)));
display(sprintf('collimate_fixed radio, %3.4f, %3.4f', ...
    fm(6:7)));
display(sprintf('flexure radio, %3.4f, %3.4f', fm(1:2)));



  
% that does it.
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obs ide off] = calcOffScan2(d)

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

% I. split them up into respective scans
% Ia.  find the start and stop index of each scan.
[s e] = findStartStop(d.array.frame.features);
betaAzAll = nan(length(s), 2, 3);
betaElAll = nan(length(s), 2, 3);
  
for m=1:length(s)
  % cut the data for each scan
  ind = zeros(size(d.array.frame.features));
  ind(s(m):e(m)) = 1;
  ind = logical(ind);
  dc  = framecut(d, ind);

  % find the observation az/el
%  indOnSource = dc.antenna0.tracker.offSource==0;
%  indAzOff    = dc.antenna0.tracker.horiz_off(:,1) == 0 & ...
%      dc.antenna0.tracker.scan_off(:,1) ==0;
%  indElOff    = dc.antenna0.tracker.horiz_off(:,2) == 0 & ...
%      dc.antenna0.tracker.scan_off(:,2) ==0;
%  f = find(indOnSource & indAzOff & indElOff);
%  if(length(f)==0)
%    f = find(indOnSource);
%  end
%  if(length(f)>1)
%    obsvals = (dc.antenna0.tracker.horiz_mount(f(1),:)) + dc.antenna0.tracker.encoder_off(f(1),:);      
%    idevals = (dc.antenna0.tracker.horiz_topo(f(1),:));
%    
%  elseif(length(f) == 1)
%    obsvals = dc.antenna0.tracker.horiz_mount(f,:) + dc.antenna0.tracker.encoder_off(f,:);    
%    idevals = (dc.antenna0.tracker.horiz_topo(f,:));
%  end


  idevals = dc.antenna0.tracker.horiz_topo(1,:);
  [obsvals(1) obsvals(2)] = pointing_model(om, idevals(1), idevals(2));
  ide.az(m) = idevals(1);
  ide.el(m) = idevals(2);
  obs.az(m) = obsvals(1);
  obs.el(m) = obsvals(2);
  
  
  
  
  
  % next we split it into the azimuth/elevation scan
  % for the elevation scan, we fit a line and remove the slope (effectively
  % a sky dip), and fit a gaussian to the remainder
  % for the azimuth scan, we just fit a gaussian to it.

  % get the offsets from the control system.
  azApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,1), dc.antenna0.receiver.utc);
  azOffSave = azApp - dc.antenna0.servo.apparent(:,1);
%  [aa, ee ] =calcAzElCs(dc, 1);
%  ao = aa - dc.antenna0.servo.az;
%  azOffSave = ao;
%  azOffSave(azOffSave<-180) = azOffSave(azOffSave<-180) + 360;
%  azOffSave(azOffSave>180) = azOffSave(azOffSave>180) - 360;
   azOffSave = -azOffSave;

  elApp = interp1(dc.antenna0.tracker.utc, ...
      dc.antenna0.tracker.horiz_topo(:,2), dc.antenna0.receiver.utc);
  elOffSave = elApp - dc.antenna0.servo.apparent(:,2);
  %elOffSave = ee - dc.antenna0.servo.el;
  elOffSave = -elOffSave;
  
  % assuming the slowest we would ever scan is at 0.2 degrees per second, a
  % az scan happens when the az speed is mostly constant, increasing, and
  % greater than 0.2 degrees/s, which is 
  azInd = (deriv(dc.antenna0.servo.az))>0.1/100;
  if(isempty(find(azInd)))
    % fudge it so it doesn't crash
    azInd(1:2) = 1;
    azInd = logical(azInd);
    azIs  = dc.antenna0.receiver.data(azInd,[1 6]);
    azIs(:) = nan;
    azFlags = ones(size(azIs));
  else
    [startaz endaz] = findStartStop(azInd);
    f = find( (endaz - startaz) == max(endaz - startaz));
    azInd = zeros(size(azInd));
    azInd(startaz(f):endaz(f)) = 1;
    azInd = logical(azInd);
    azIs  = dc.antenna0.receiver.data(azInd,[1 6]);
    azIs(:,1) = smooth(azIs(:,1));
    azIs(:,2) = smooth(azIs(:,2));
    
    azFlags       = dc.flags.fast(azInd,[1 3]);
    azIs(azFlags) = nan;
  end
  
  % if amplitude is negative, no peak in data
  
  % same goes for elevation
  elInd = deriv(dc.antenna0.servo.el)>0.1/100;
  if(isempty(find(elInd)))
    elInd(1:2) = 1;
    elInd = logical(elInd);
    elIs  = dc.antenna0.receiver.data(elInd,[1 6]);
    elIs(:) = nan;
    elFlags = ones(size(elIs));    
  else
    [startel endel] = findStartStop(elInd);
    f = find( (endel-startel) == max(endel - startel));
    elInd = zeros(size(elInd));
    elInd(startel(f):endel(f)) = 1;
    elInd = logical(elInd);
    elIs  = dc.antenna0.receiver.data(elInd,[1 6]);
    elIs(:,1) = smooth(elIs(:,1));
    elIs(:,2) = smooth(elIs(:,2));
    
    elFlags       = dc.flags.fast(elInd, [1 3]);
    elIs(elFlags) = nan;
  end
    
  % if most of the scan is flagged, don't bother fitting.
  azFlagPer = length(find(azFlags))./length(azIs(:));
  elFlagPer = length(find(elFlags))./length(elIs(:));
  if(azFlagPer < 0.15 & elFlagPer < 0.15)
    doFit = 1;
  else
    doFit = 0;
  end


  if(doFit)
    testVals = [120:30:800];
    for mm=1:2
      
      if(~isempty(find(azInd)))
	% fit a gaussian about the max in both cases
	azOff = azOffSave;
	azI = azIs(:,mm);
	f = find(azI == max(azI));
	azOff = azOff(azInd);
	ampVal = nan;
	stopAz = 0;
	
	azIndex = 1;
	while(stopAz==0)
	  azFit = find(azInd);
	  if(f<=testVals(azIndex) | (f+testVals(azIndex))>=length(azI))
	    badAzPoint = 1;
	    azFit = [];
	    xoff = nan;
	    azIndex = azIndex+1;
	    betaAz(1:3) = nan;
	  else
	    azFit = azOff((f-testVals(azIndex):f+testVals(azIndex)));
	    azIFit = azI(f-testVals(azIndex):f+testVals(azIndex));
	    indgood= zeros(size(azI));
	    indgood(f-testVals(azIndex):f+testVals(azIndex)) = 1;
	    
	    
	    display(sprintf('testing value: %d', testVals(azIndex)));
	    
	    OPTIONS = optimset('Display', 'off', 'MaxIter', 200000, ...
		'MaxFunEvals', 20000, 'TolX', 0.00005, 'TolFun', 0.00005);
	    LB = [0 azFit(testVals(azIndex))-0.2  -10 -inf -inf];
	    UB = [max(azIFit) azFit(testVals(azIndex))+0.2 10 inf inf];
	    x0 = [1 azFit(testVals(azIndex)-1) 1 1 1];
	    [betaAz resnorm residual exitflag] = lsqnonlin(@(x) mchisq2(azIFit, gaussfit(x, azFit)), x0, LB, UB, OPTIONS);
	    
	    amplitude = betaAz(1);
	    xoff      = betaAz(2);
	    sigma     = betaAz(3);
	    xoffIndex = find(abs(azFit - xoff) == min(abs(azFit - xoff)));
	    ampVal    = azIFit(xoffIndex);
	    
	    
	    % other terms remove a baseline
	    if (amplitude<0 | abs(xoff)>8 | resnorm>500)
	      badAzPoint = 1;
	    else
	      badAzPoint = 0;
	    end 
	    azFit = gaussfit(betaAz, azOff);
	    
	    if(abs(xoff)<4)
	      stopAz = 1;
	    else
	      azIndex = azIndex+1;
	    end
	  end
	  
	  if(azIndex>length(testVals))
	    stopAz = 1;
	  end
	end
	% save for weighting
	if(azIndex >= length(testVals))
	  azIndex = azIndex-1;
	end
	azWidth(m,mm) = testVals(azIndex);
	azMax(m,mm)   = ampVal;
	if(~badAzPoint)
	  azBase(m,mm)  = median(azI(~indgood));
	  azRms(m,mm)   = std(azI(~indgood));
	else
	  azBase(m,mm)  = nan;
	  azRms(m,mm)   = nan;
	end
	betaAzAll(m,mm,:) = betaAz(1:3);
      end
      
      
      if(~isempty(elInd))
	elOff = elOffSave;
	elI = elIs(:,mm);
	% first we remove the offset from slewing the atmosphere
	% from the sky dip code
	x = 1./(sind(dc.antenna0.servo.el(elInd)));
	y = abs(elI);
	x(isnan(y)) = [];
	y(isnan(y)) = [];
	[tau Tground] = linfit(x,y);
	% remove the effect
	elI = elI - Tground(1) - tau(1)./sind(dc.antenna0.servo.el(elInd));
	
	% fit a gaussian about the max 
	frange = elOff(elInd);
	indrange = abs(frange)<4;
	f = find(elI(indrange) == max(elI(indrange))) + find(indrange,1);
	elOff2 = elOff(elInd);
	
	stopEl = 0;
	elIndex = 1;
	while(stopEl==0)
	  elFit = find(elInd);
	  if(f<=testVals(elIndex) | (f+testVals(elIndex))>=length(elI))
	    badElPoint = 1;
	    elFit = [];
	    yoff = nan;
	    elIndex = elIndex+1;
	    betaEl = [nan nan nan];
	  else
	    elFit = elOff2((f-testVals(elIndex):f+testVals(elIndex)));
	    elIFit = elI(f-testVals(elIndex):f+testVals(elIndex));
	    indgood= zeros(size(elI));
	    indgood(f-testVals(elIndex):f+testVals(elIndex)) = 1;
	    
	    display(sprintf('el testing value: %d', testVals(elIndex)));
	    [betaEl] = nlinfit(elFit, elIFit, @gaussfit, [1 elFit(testVals(elIndex)) 1 1 1]);
	    amplitude = betaEl(1);
	    yoff      = betaEl(2);
	    sigma     = betaEl(3);
	    yoffIndex = find(abs(elFit - yoff) == min(abs(elFit - yoff)));
	    ampVal    = elIFit(yoffIndex);	    
	    
	    % other terms remove a baseline
	    if (amplitude<0 | abs(yoff)>8)
	      badElPoint = 1;
	    else
	      badElPoint = 0;
	    end 
	    elFit = gaussfit(betaEl, elOff2);
	    
	    if(abs(yoff)<1)
	      stopEl = 1;
	    else
	      elIndex = elIndex+1;
	    end
	  end
	  if(elIndex>length(testVals))
	    stopEl = 1;
	  end
	  
	end
	% save for weighting
	if(elIndex >= length(testVals))
	  elIndex = elIndex - 1;
	end
	elWidth(m,mm) = testVals(elIndex);
	elMax(m,mm)   = ampVal;
	if(~badElPoint)
	  elBase(m,mm)   = median(elI(~indgood));
	  elRms(m,mm)    = std(elI(~indgood));
	else
	  elBase(m,mm)   = nan;
	  elRms(m,mm)    = nan;
	end
	betaElAll(m,mm,:) = betaEl(1:3);      
      end
      
      display('here');
      
      display(sprintf('Point %d of %d', m, length(s)));
      
      
      clf
      %if(~any(isnan(betaAzAll(m,mm,:))))
      
      subplot(2,1,1)
      plot((azOff)*cos(dc.antenna0.servo.el(1)*pi/180), azI');
      xlabel('Az Offset (degrees)');
      ylabel('Intensity');
      if(~isempty(azFit))
	hold on; plot(azOff*cos(dc.antenna0.servo.el(1)*pi/180), azFit, 'r'); hold off
	axis([min(azOff) max(azOff) min(azI) max(azI)]);
      end
      
      %end
      
      %if(~any(isnan(betaElAll(m,mm,:))))
      subplot(2,1,2)
      plot(elOff(elInd), elI);
      xlabel('El Offset (degrees)');
      ylabel('Intensity');
      if(~isempty(elFit))
	hold on; plot(elOff(elInd), elFit, 'r'); hold off
	axis([min(elOff) max(elOff) min(elI) max(elI)]);
      end
      %end
      eval(sprintf('gtitle(''Pointing scan on source %s'');', dc.antenna0.tracker.source{1}));
      
      
      
      display(sprintf('Source Name: %s', dc.antenna0.tracker.source{1}));
      r = input('Keep This point? [y/n]: ', 's');
      
      if(strcmp(r, 'y') || strcmp(r, 'Y'))
	badPoint(m,mm) = 0;
      elseif(strcmp(r, 'n') ||  strcmp(r, 'N'))
	badPoint(m,mm) = 1;
      elseif(strcmp(r, 'k'))
	keyboard;
      end
    end
    display('here 3');
    
    clear azIs;
    clear elIs;
    
    pause(0.1);
  else
    display('All data in this cross have been flagged');
    azWidth(m,1:2) = nan;
    azMax(m,1:2)   = nan;
    azBase(m,1:2)  = nan;
    azRms(m,1:2)   = nan;
    betaAzAll(m,1:2,1:3) = nan;
    elWidth(m,1:2) = nan;
    elMax(m,1:2)   = nan;
    elBase(m,1:2)  = nan;
    elRms(m,1:2)   = nan;
    betaElAll(m,1:2,1:3) = nan;
  end
end

badPoint = logical(badPoint);
elOff = betaElAll(:,:,2);
azOff = betaAzAll(:,:,2);
elOff(badPoint) = nan;
azOff(badPoint) = nan;

% first we assign a weight with the source.
azSig = (azMax - azBase)./azRms;
azWidth = abs(betaAzAll(:,:,3)/2);
azErr = azWidth./azSig;


elSig = (elMax - elBase)./elRms;
elWidth = abs(betaElAll(:,:,3)/2);
elErr   = elWidth./elSig;

% throw out values where the error is large or negative.
elBad = elErr < 0 | elErr > 1;
azBad = azErr < 0 | azErr > 1;

elOff(elBad) = nan;
azOff(azBad) = nan;
elErr(elBad) = nan;
azErr(elBad) = nan;

% if hte values for both channels aren't within 10 arcmin, toss them.
dEl = abs(diff(elOff, [], 2));
dAz = abs(diff(azOff, [], 2).*cos(obs.el*pi/180)');

f = find(dEl>10/60 | dAz>10/60);
elOff(f,:) = nan;
azOff(f,:) = nan;
elErr(f,:) = nan;
azErr(f,:) = nan;


display(sprintf('Throwing out %d points', length(f)));

[off.az  off.azErr] = vwstat(azOff, azErr, 2);
[off.el  off.elErr] = vwstat(elOff, elErr, 2);

display('Done with all crosses');

% add final offests to the observed az/el
obs.az = obs.az + off.az';
obs.el = obs.el + off.el';

% throw out the data that is bad
ind = isnan(obs.az);
obs.az(ind) = [];
obs.el(ind) = [];
ide.az(ind) = [];
ide.el(ind) = [];
off.az(ind) = [];
off.el(ind) = [];
off.azErr(ind) = [];
off.elErr(ind) = [];


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obs ide] = calcOffRaster(d)

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
function [s e] = findStartStop(features)

% start indices are when the difference in the feature value increases, stop
% is when it decreases
s = find(diff(double(features))>0);

% end is when it decreases
e = find(diff(double(features))<0);

% last one
lastOneBad = 0;
if(length(e)<length(s))
  % the dimensions don't match
  if(s(1)<e(1))
    % if start is less than end, that means we're usually missing hte last
    % endpoint
    e(length(s)) = length(features);
    lastOneBad = 1;
  else
    s(1) = [];
  end
elseif(length(s)<length(e))
  s(1) = [];
end
s = s+1;

goodPoint = zeros(size(s));
for m=1:length(s)
  thisFeat = unique(features(s(m):e(m)));
  if(thisFeat == 129 | thisFeat ==128)
    goodPoint(m) = 1;
  end
end  

if(lastOneBad==1)
  goodPoint(length(s)) = 0;
end


s = s(logical(goodPoint));
e = e(logical(goodPoint));

return;


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
function [obs ide off] = calcOffCross(d)

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
    keyboard;
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

  figure(1)
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
