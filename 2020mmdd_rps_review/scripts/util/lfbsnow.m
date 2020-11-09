function snow=lfbsnow(tags)
%snow=lfbsbow(tags) looks for snow accumulated on window by processing LFB
%IR camera images for tags given
%   Looks at all LFB images between the start of the first tag and the end
%   of the second tag.  Tags must be cell.
%   Applies patented snow finding algorithm and outputs result in arbitrary
%   snow units (ASU).  Snow will be very recognizable but as a rule of
%   thumb, anything over 1700 ASU's is probably snow.
%
%   Output:
%       snow.acc  = snow accumulation for set of tags
%       snow.jpgs = cell structure containing the LFB images
%       snow.t    = time axis
%
%       example code for plotting LFB images:
%           figure; imagesc(snow.jpgs{xx}); colormap('Gray'); axis off; caxis([0 200]);

if(~iscell(tags))
    tags={tags};
end

ri=get_run_info(tags);
tstart=min(datenum(ri.tstart,'dd-mmm-yyyy:HH:MM:SS'));
tend=max(datenum(ri.tend,'dd-mmm-yyyy:HH:MM:SS'));

% Get list of LFB images between dates
[ljpgs tjpgs]=get_jpgs(tstart,tend);

if(~isempty(ljpgs))
  snow.acc=findsnow(ljpgs);
  snow.jpgs=ljpgs;
  snow.t=tjpgs;
else
  snow.acc=[];
  snow.jpgs=[];
  snow.t=[];
end

end

function [obj dateobj]=get_jpgs(ts,te)
% For the timespan of the tags, load the jpgs from the LFB camera

% Get list of jpgs using ls for each day that the tags span
f=[];
for dn=floor(ts):ceil(te)
  fnames=sprintf('snow/jpg/%s*.jpg',datestr(dn,'yyyymmdd'));
  f=[f;dir(fnames)];
end
fnames=char(f.name);
fnames=cellstr(fnames(:,1:15));
dnf=datenum(fnames,'yyyymmdd_HHMMSS');

% Pick out only the images between tstart and tend
pickind=dnf>=ts & dnf<=te;
list=fnames(pickind);
dnf=dnf(pickind);

% load jpgs
obj=[];
dateobj=[];
if(~isempty(list))
    xxi=1;
    for xx=1:length(list)
        iserr=0;
        try
            obj{xxi}=imread(['snow/jpg/' list{xx} '.jpg']);
            dateobj(xxi)=dnf(xx);
        catch err
            disp(['Problem loading ' list{xx}]);
            iserr=1;
            continue
        end
    
        if(~iserr)
            xxi=xxi+1;
        end
    end
else
    disp(['ERROR: no LFB images found between ' datestr(ts) ' and ' datestr(te)]);
end

end

function avgfft=findsnow(obj)
% Apply snow finding algorithm for all the LFB images in range

for jj=1:length(obj)
    % NaN out stripes    
    obj{jj}(240:241,:)=NaN;
    obj{jj}(:,256:257)=NaN;
    
    % Determine if the image is crap
%     avgbw(jj)=nanmean(nanmean(im2bw(obj{jj}(221:end,:),.08)));
    
    % take fft in region below the flare from the Al membrane
    fftobj{jj}=fftshift(fft2(obj{jj}(221:end,:)));
    
    % remove stripes in fft
    fftobj{jj}(129:132,:)=NaN;
    fftobj{jj}(:,255:260)=NaN;

    % take avg over small section of fft around center--this should
    % probably be done over a circle but a rectangle is trivial to code
    % and works fine
    avgfft(jj)=nanmean(nanmean(abs(fftobj{jj}(130-30:130+30,256-55:256+55))));
end
end
