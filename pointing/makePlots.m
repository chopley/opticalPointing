[s w] = unix('ls --color=never -lh /mnt/data/sjcm/jup_sec');

entVals = strfind(w, '.fit');

for m=1:length(entVals)
  entry{m} = w(entVals(m)-15:entVals(m)+3);
end

%eloff = 0.75:0.25:4.76;
%azoff = 0.15:0.25:4.16;
eloff = -2:0.15:2;
azoff = -2:0.15:2;

keyboard;
index = 1;
fullMap = [];
for m=1:length(eloff)
  thisRow = [];
  for n=1:length(azoff)
    eval(sprintf('thisFits = fitsread2(''/mnt/data/sjcm/jup_sec/%s'');', entry{index}));
%thisFits(thisFits<0) = nan;
%thisFits(thisFits>50) = nan;
    tf(:,:,n) = thisFits;
    index  = index+1;
  end
    if m==16
    display('at right row');
keyboard;
end

    
  badPix = mean(tf,3)>100 | std(tf,0,3)==0;
flat = ones(size(thisFits));
flat(badPix) = nan;
flat = repmat(flat, [1 1 n]);
tf = tf./flat;
thisRow = [];
for n=1:length(azoff)
  thisRow = [thisRow tf(:,:,n)];
end
imagesc(thisRow);
eval(sprintf('title(''Elevation Row %d: %2.2f'');', m, eloff(m)));
meanRowAz(m) = mean(tf(:));
colorbar;

pause(1);
end

keyboard;
