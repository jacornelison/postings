function [mask_100 mask_150] = ffbm_maskMaker(mirrorpos)
% function [mask_100 mask_150] = ffbm_maskMaker(mirrorpos)
% 
% Determines which detectors are on/off the mirror
% Mask detailed here:
% http://bicep.caltech.edu/~spuder/analysis_logbook/
%        analysis/20131217_compositebm/index_position.html
%
% Mask value 0 allows that pixel to be included in measurements
% Mask value 1 are detectors with absolutely no beams
% Mask value 2 are detectors with beams, but are falling off the mirror
%
% New for 2014: position flag
%               2012/2013/2014a: 'back'
%               2014b: 'front'

% 10 positions, detailed in posting, or in dk_rx_pos_ffbm
% Position 3 is at base of drum.
mask = zeros(10,528);

if ~exist('mirrorpos','var')
  mirrorpos = 'back';
end

% Sometimes we're stupid and use 'forward' instead of 'front'
if strcmp(mirrorpos,'forward')
  mirrorpos = 'front';
end

% Load up fp_data file and get pixel locations
ind_100 = 1:528;
ind_150 = 529:1056;
p = get_array_info(20140301);
x_100 = p.pix_phys_x(ind_100); x_150 = p.pix_phys_x(ind_150);
y_100 = p.pix_phys_y(ind_100); y_150 = p.pix_phys_y(ind_150);
ti_100 = p.tile(ind_100);      ti_150 = p.tile(ind_150);

% 150 is hilariously counterintuitive, so mirror it to match 100
for i = 1:528
  switch ti_150(i)
    case 1
      x_150(i) = -x_150(i);
    case 2
      y_150(i) = -y_150(i);
    case 3
      x_150(i) = -x_150(i);
    case 4
      y_150(i) = -y_150(i);
  end
end

if strcmp(mirrorpos,'back')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 1
  
  % Tile 1 = 2
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  
  % Tile 2 = 0
  
  % Tile 3 row 1 = 2
  for ro=1:1
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  % Tile 3 row 2 col 5:8 = 2
  for ro=2:2
    for co=5:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  % Tile 3 row 3 col 7:8 = 2
  for ro=3:3
    for co=7:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  % Tile 3 col 8 row 4:5 = 2
  for ro=4:5
    for co=8:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  
  % Tile 4 = 2
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Position 2

  % Tile 1 = 0
  
  % Tile 2 = 0
  
  % Tile 3 = 0
  
  % Tile 4 row 1:4 col 8 = 2
  for ro=1:4
    for co=8:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  % Tile 4 row 1:3 col 7 = 2
  for ro=1:3
    for co=7:7
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  % Tile 4 row 1:2 col 6 = 2
  for ro=1:2
    for co=6:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  % Tile 4 row 1 col 4:5 = 2
  for ro=1:1
    for co=4:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Position 3 -- all good
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Position 4 
  
  % Tile 1 = 0
  
  % Tile 2 = 0
  
  % Tile 3 row 1:4 col 1 = 2
  for ro=1:4
    for co=1:1
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  % Tile 3 row 1:3 col 2 = 2
  for ro=1:3
    for co=2:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  % Tile 3 row 1:2 col 3 = 2
  for ro=1:2
    for co=3:3
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  % Tile 3 row 1 col 4:5 = 2
  for ro=1:1
    for co=4:5
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  
  % Tile 4 = 0
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 5

  % Tile 1 row 1:4 col 1:2 = 2
  for ro=1:4
    for co=1:2
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  % Tile 1 row 1:6 col 3 = 2
  for ro=1:6
    for co=3:3
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  % Tile 1 row 1:7 col 4 = 2
  for ro=1:7
    for co=4:4
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  % Tile 1 row 1:8 col 5:8 = 2
  for ro=1:8
    for co=5:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end

  % Tile 2 = 2
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 3 = 2
  for ro=1:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 4 = 0

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 6

  % Tile 1 = 2
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end

  % Tile 2 col 1 = 2
  for ro=1:8
    for co=1:1
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  % Tile 2, col 2:8 = 1
  for ro=1:8
    for co=2:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1;      
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end
  
  % Tile 3, col 1:4 = 1
  for ro=1:8
    for co=1:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1;   
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end
  % Tile 3, col 5 row 5:8 = 1
  for ro=5:8
    for co=5:5
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1;
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end

  % Tile 3,4 = 2 (except 1s from above)
  for ro=1:8
    for co=1:8
      for ti=3:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end    
	if mask(6,B)==0 mask(6,B)=2; end
	B=det2gcp(ro,co,'B',ti);
	if mask(6,B)==0 mask(6,B)=2; end
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Positions 7:9 = 1
  mask(7:9,:)=1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 10

  % Tile 1 = 1
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end
  
    % Tile 2 = 2
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  
  % Tile 4 col 5:8 = 1
  for ro=1:8
    for co=5:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end
  % Tile 4 col 4 row 5:8 = 1
  for ro=5:8
    for co=4:4
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end

  % Tile 3,4 = 2 (except 1s from above) 
  for ro=1:8
    for co=1:8
      for ti=3:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	if mask(10,B)==0 mask(10,B)=2; end 
	B=det2gcp(ro,co,'B',ti);
	if mask(10,B)==0 mask(10,B)=2; end 
      end
    end
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % OLD STUFF
  % Plot the mask
  % Now do a tile map with sigmas 
  %tileMapA=NaN(17,17);
  %tileMapB=NaN(17,17);
  %for ii=1:528
  %  if ii==0, ii=528; end
  %  if ii==528, index=0; else index=ii; end
  %  [ro co pol ti]=gcp2det(index); 
  %  if (ti==1 && strcmp(pol,'A'))
  %    tileMapA(ro,co)=mask(6,ii);
  %  elseif (ti==1 && strcmp(pol,'B'))
  %    tileMapB(ro,co)=mask(6,ii);
  %  elseif (ti==2 && strcmp(pol,'A'))
  %    tileMapA(ro,co+9)=mask(6,ii);
  %  elseif (ti==2 && strcmp(pol,'B')) 
  %    tileMapB(ro,co+9)=mask(6,ii);
  %  elseif (ti==3 && strcmp(pol,'A'))
  %    tileMapA(18-ro,18-co)=mask(6,ii);
  %  elseif (ti==3 && strcmp(pol,'B')) 
  %    tileMapB(18-ro,18-co)=mask(6,ii);
  %  elseif (ti==4 && strcmp(pol,'A'))
  %    tileMapA(18-ro,9-co)=mask(6,ii);
  %  elseif (ti==4 && strcmp(pol,'B')) 
  %    tileMapB(18-ro,9-co)=mask(6,ii);
  %  end
  %end
  % END OLD STUFF
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

elseif strcmp(mirrorpos,'front')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 1
  
  % Tile 1 = 2 (all except left column)
  for ro=1:8
    for co=2:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end

  % Tile 2 = 2 (all except right column)
  for ro=1:8
    for co=1:6
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=7:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=1; 
      end
    end
  end
  
  % Tile 3 = 2
  for ro=1:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  
  % Tile 4 = 2 (all except left columns)
  for ro=1:8
    for co=1:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 2
  
  % Tile 1 = 1 (upper triangle)
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=1; 
      end
    end
  end
  % Tile 1 = 2 (lower triangle)
  for ro=1:8
    for co=1:1
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=2:8
    for co=2:2
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=3:8
    for co=3:3
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=4:8
    for co=4:4
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=5:8
    for co=5:5
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=6:8
    for co=6:6
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=7:8
    for co=7:7
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=8:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  
  % Tile 2 = 1
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=1; 
      end
    end
  end
  
  % Tile 3 = 1 (upper triangle)
  for ro=1:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=1; 
      end
    end
  end
  % Tile 3 = 2 (lower triangle)
  for ro=1:8
    for co=6:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:7
    for co=5:5
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=3:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:5
    for co=2:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=1:1
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  
  % Tile 4 = 2
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 3
  
  % Tile 1,2 = 1
  for ro=1:8
    for co=1:8
      for ti=1:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(3,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(3,B)=1; 
      end
    end
  end
  
  % Tile 3,4 = 1 (top rows) 
  for ro=8:8
    for co=1:8
      for ti=3:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(3,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(3,B)=1; 
      end
    end
  end
  for ro=1:7
    for co=1:8
      for ti=3:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(3,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(3,B)=2; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 4
  
  % Tile 1 = 1
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=1; 
      end
    end
  end
  
  % Tile 2 = 1 (upper left triangle)
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=1; 
      end
    end
  end
  % Tile 2 = 2 (lower right triangle)
  for ro=2:8
    for co=7:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=3:8
    for co=6:6
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=4:8
    for co=5:5
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=5:8
    for co=4:4
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=6:8
    for co=3:3
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=7:8
    for co=2:2
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=1:1
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
 
  % Tile 3 = 2
  for ro=1:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  
  % Tile 4 = 1 (upper left triangle)
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=1; 
      end
    end
  end
  % Tile 4 = 2 (lower right triangle)
  for ro=1:8
    for co=1:3
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:7
    for co=4:4
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=5:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=5:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:5
    for co=6:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=7:7
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:3
    for co=8:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 5
  
  % Tile 1 = 1 (left columns)
  for ro=1:8
    for co=1:2
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  for ro=1:8
    for co=3:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 2 = 2 (except right column)
  for ro=1:8
    for co=1:7
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 3 = 2 (except right columns)
  for ro=1:8
    for co=3:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 4 = 2
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 6
  
  % Tile 1 = 2 (upper 5 rows) 
  for ro=1:5
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end

  % Tile 2 = 2 (most of it)
  for ro=1:5
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end

  for ro=6:6
    for co=3:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=7:7
    for co=4:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=5:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  
  % Tile 3 = 2 (right half)
  for ro=1:8
    for co=1:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end

  % Tile 4 = 0
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 7
  
  % Tile 1 = 2 (upper half, bottom half upper triangle)
  for ro=1:4
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=5:5
    for co=5:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=6:6
    for co=6:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=7:7
    for co=7:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=8:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  
  % Tile 2 = 2
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  
  % Tile 3 = 2 (upper right triangle)
  for ro=8:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=7:7
    for co=1:7
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=6:6
    for co=1:6
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=5:5
    for co=1:5
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=4:4
    for co=1:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=3:3
    for co=1:3
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=2:2
    for co=1:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=1:1
    for co=1:1
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end

    
  % Tile 4 = 0

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 8
  
  % Tile 1 = 2 
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(8,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(8,B)=2; 
      end
    end
  end

  % Tile 2 = 2 
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(8,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(8,B)=2; 
      end
    end
  end
  
  % Tile 3 = 2 (top half)
  for ro=5:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(8,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(8,B)=2; 
      end
    end
  end
  
  % Tile 4 = 2 (top half)
  for ro=5:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(8,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(8,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 9
    
  % Tile 1 = 2, 1 (upper left corner)
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=1:2
    for co=1:2
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=1; 
      end
    end
  end
  
  % Tile 2 = 2 (upper left triangle)
  for ro=1:1
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=2:2
    for co=1:7
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=3:3
    for co=1:6
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=4:4
    for co=1:5
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=5:5
    for co=1:4
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=6:6
    for co=1:3
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=7:7
    for co=1:2
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=1:1
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end

  % Tile 3 = 0
  
  % Tile 4 = 2 (upper left triangle)
  for ro=8:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=7:7
    for co=2:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=6:6
    for co=3:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=5:5
    for co=4:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=4:4
    for co=5:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=3:3
    for co=6:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=2:2
    for co=7:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end
  for ro=1:1
    for co=8:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(9,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(9,B)=2; 
      end
    end
  end


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 10
  
  % Tile 1 = 2 (left wedge)
  for ro=1:8
    for co=1:6
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=7:7
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:3
    for co=8:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end

  % Tile 2 = 0
  
  % Tile 3 = 0
  
  % Tile 4 = 2 (left wedge)
  for ro=1:8
    for co=6:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=4:8
    for co=5:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=5:8
    for co=4:4
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  

elseif strcmp(mirrorpos,'forwardest')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 1
  
  % Tile 1 = 0,1, or 2...
  for ro=1:4
    for co=3:3
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=4:7
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=8:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=1; 
      end
    end
  end
  for ro=5:8
    for co=8:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end

  % Tile 2 = 1 
  for ro=1:8
    for co=1:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=1; 
      end
    end
  end
  
  % Tile 3 = 1
  for ro=1:8
    for co=1:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=1; 
      end
    end
  end
  
  % Tile 4 = 0 or 2
  for ro=1:8
    for co=1:2
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  for ro=6:8
    for co=3:3
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end
  for ro=8:8
    for co=4:4
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(1,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(1,B)=2; 
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 2  
  
  % Tile 1,2,3 = 1
  for ro=1:8
    for co=1:8
      for ti=1:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=1; 
      end
    end
  end
  
  % Tile 4 = 1 or 2
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=1; 
      end
    end
  end
  for ro=1:1
    for co=2:2
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:2
    for co=3:3
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:3
    for co=4:4
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=5:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:5
    for co=6:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=7:7
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  for ro=1:7
    for co=8:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(2,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(2,B)=2; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 3
  
  % All tiles = 1
  for ro=1:8
    for co=1:8
      for ti=1:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(3,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(3,B)=1; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 4
  
  % All tiles = 1, except for lower triangle of tile 3
  for ro=1:8
    for co=1:8
      for ti=1:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=1; 
      end
    end
  end
  
  % Tile 3 = 2 in lower triangle
  for ro=1:6
    for co=1:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:5
    for co=3:3
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=4:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:3
    for co=5:6
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  for ro=1:2
    for co=7:7
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(4,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(4,B)=2; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 5
  
  % Tile 1 = 1
  for ro=1:8
    for co=1:8
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  
  % Tile 2 = 0,1, or 2...
  for ro=1:8
    for co=1:1
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  for ro=1:4
    for co=2:2
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  for ro=5:8
    for co=2:2
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  for ro=1:1
    for co=3:3
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  for ro=2:8
    for co=3:3
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=4:5
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=6:6
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=7:7
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 3 = 0 or 2
  for ro=6:8
    for co=4:4
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=5:8
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=2; 
      end
    end
  end
  
  % Tile 4 = 1
  for ro=1:8
    for co=1:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(5,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(5,B)=1; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 6
  
  % Tile 1 = 0 or 2
  for ro=1:8
    for co=1:2
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=3:8
    for co=3:4
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=5:8
    for co=5:5
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=6:8
    for co=6:6
      for ti=1:1
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  
  % Tiles 2,3 = 0
  
  % Tile 4 = 0,1, or 2...
  for ro=1:5
    for co=1:1
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=2:5
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=1:2
    for co=6:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end
  for ro=3:8
    for co=6:6
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=1:5
    for co=7:7
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end
  for ro=6:8
    for co=7:7
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=8:8
      for ti=4:4
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(6,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(6,B)=1; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 7
  
  % All Tiles = 0 except for chunk of Tile 2
  
  % Tile 2 = 2 in all but one corner
  for ro=1:5
    for co=1:1
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=1:6
    for co=2:2
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=1:7
    for co=3:3
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=4:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(7,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(7,B)=2; 
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 8
  
  % All Tiles = 0
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 9
  
  % All Tiles = 0
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Position 10
  
  % Tile 1 = 0
  
  % Tile 2 = 2 in right half
  for ro=7:8
    for co=4:4
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=5:8
    for co=5:6
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=4:8
    for co=7:7
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=8:8
      for ti=2:2
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  
  % Tile 3 = 0,1, or 2...
  for ro=1:8
    for co=1:1
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end
  for ro=1:5
    for co=2:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end
  for ro=6:8
    for co=2:2
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:2
    for co=3:3
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=1; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=1; 
      end
    end
  end
  for ro=3:8
    for co=3:3
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:8
    for co=4:6
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  for ro=1:4
    for co=7:7
      for ti=3:3
	B=det2gcp(ro,co,'A',ti);
	if B==0 B=528; end
	mask(10,B)=2; 
	B=det2gcp(ro,co,'B',ti);
	mask(10,B)=2; 
      end
    end
  end
  
  % Tile 4 = 0
  
end





% All the above have been saved in gcp index with gcp0 in element 528
% Shift indices for matching rgls
mask_150 = circshift(mask,[0,1]); 

% Make the 100 GHz mask
mask_100 = zeros(10,528);
for i = 1:10
  for j = 1:528
    % For ease of plotting, NaN out 100 GHz pixels which don't exist
    if isnan(x_100(j))
      mask_100(i,j) = NaN;
    end
    % Do the same for 150, I guess
    if isnan(x_150(j))
      mask_150(i,j) = NaN;
    end
    
    % Check if we're a real pixel
    if ~isnan(x_150(j))
      % Check if the 150 GHz mask is 1 or 2
      if mask_150(i,j) == 1 || mask_150(i,j) == 2
	% Find distance of all 100 GHz pixels from it
	dist = sqrt( (x_100 - x_150(j)).^2 + (y_100 - y_150(j)).^2 );
	ind = find( dist == min(dist) );
	% If it's zero in the 100 GHz mask, change it
	if mask_100(i,ind) == 0
	  mask_100(i,ind) = mask_150(i,j);
	end
      end
    end
  end
end



return

