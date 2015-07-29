[s w] = unix('ls --color=never -lh /home/cbass');

entVals = strfind(w, '.fit');

for m=1:length(entVals)
  entry{m} = w(entVals(m)-15:entVals(m)+3);
end

eloff = -2:0.15:2.01;
azoff = -2:0.15:2.01;

keyboard;
index = 1;
fullMap = [];
for m=1:length(eloff)
  thisRow = [];
  for n=1:length(azoff)
    eval(sprintf('thisFits = fitsread2(''/home/cbass/%s'');', entry{index}));
thisFits(thisFits<0) = nan;
thisFits(thisFits>50) = nan;
    thisRow = [thisRow thisFits];
    tf(:,:,n) = thisFits;
    index  = index+1;
  end
  imagesc(thisRow);
  eval(sprintf('title(''Elevation Row %d: %2.2f'');', m, eloff(m)));
  meanRowAz(m) = mean(tf(:));
colorbar;

pause(1);
end

keyboard;
