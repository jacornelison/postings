function nfbm_plotter(x,y,px,sq,ramin,ra,xra,gauss,rotate,theta,orto)
% x = x coord of scan
% y = y coord of scan
% px(row, col, tile) = structure containing element px.a_map and px.b_map
% ra = z range for plots
% xra = x and y range, used for xlim and ylim
% rotate - set keyword to rotate images
% theta = angle by which image is rotated
% orto - set to 'small' to generate only thumbnails, set to 'large' to only
% generate large maps.

figure(1)
cf = gcf;
%set(cf,'Visible','off')

if nargin < 10
  orto = [];
end

if nargin < 7
  gauss = [];
end

if nargin < 8
  rotate = [];
end

% px = structure into which all beam files and beam parameters are stored
if nargin == 4
  ra = 1;
end

if(1)
  % NaN out pixels with values above the threshold
  threshold = 1000;
  for rr = 1:8
    for cc = 1:8
      for ti = 1:4
	for ii = 1:length(px(rr,cc,ti).a_map(:,1))
	  for jj = 1:length(px(rr,cc,ti).a_map(1,:))
	    if px(rr,cc,ti).a_map(ii,jj) > threshold
	      px(rr,cc,ti).a_map(ii,jj) = NaN;
	    end
	    if px(rr,cc,ti).b_map(ii,jj) > threshold
	      px(rr,cc,ti).b_map(ii,jj) = NaN;
	    end
	  end
	end
	if length(px(rr,cc,ti).a_map) > 1	
	  for ii = 1:2:length(px(rr,cc,ti).a_map(:,1))
	    for jj = 1:2:length(px(rr,cc,ti).b_map(1,:))
	      AVec = [px(rr,cc,ti).a_map(ii,jj);px(rr,cc,ti).a_map(ii+1,jj);px(rr,cc,ti).a_map(ii,jj+1);px(rr,cc,ti).a_map(ii+1,jj+1)];
	      BVec = [px(rr,cc,ti).b_map(ii,jj);px(rr,cc,ti).b_map(ii+1,jj);px(rr,cc,ti).b_map(ii,jj+1);px(rr,cc,ti).b_map(ii+1,jj+1)];
	      A((ii+1)/2,(jj+1)/2) = nanmean(AVec);
	      B((ii+1)/2,(jj+1)/2) = nanmean(BVec);
	    end
	  end
	  px(rr,cc,ti).a_map = A;
	  px(rr,cc,ti).b_map = B;
	  
	  for ii = 1:length(px(rr,cc,ti).a_map(:,1))
	    for jj = 1:length(px(rr,cc,ti).a_map(1,:))
	      if isnan(px(rr,cc,ti).a_map(ii,jj))
		px(rr,cc,ti).a_map(ii,jj) = threshold;
	      end
	      if isnan(px(rr,cc,ti).b_map(ii,jj))
		px(rr,cc,ti).b_map(ii,jj) = threshold;
	      end
	    end
	  end
	end
      end
    end
  end
end

% Individual plots
% Thumbnails

if strcmp(orto,'small') || isempty(orto);
  if strcmp(orto,'small')
    disp('Skipping large')
  end
  
  disp('Small maps')
  for rr=1:8
    for cc=1:8
      for ti=1:4
	if strcmp(rotate,'rotate')
	  A=imrotate(px(rr,cc,ti).a_map,theta,'bicubic');
	  B=imrotate(px(rr,cc,ti).b_map,theta,'bicubic');
	else
	  A=px(rr,cc,ti).a_map;
	  B=px(rr,cc,ti).b_map;
	end

	imagesc(x,y,A)
	set(gca,'YDir','normal')
	%xlim([-xra,xra])
	%ylim([-xra,xra])
	axis square
	%caxis([ramin,ra])
	colormap('JET')
	set(gca,'xtick',[],'ytick',[])
	img_str=['small/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_A_sm'];
	print('-dpng', img_str);
	clf
	imagesc(x,y,B)
	set(gca,'YDir','normal')
	%xlim([-xra,xra])
	%ylim([-xra,xra])
	axis square
	caxis([ramin,ra])
	colormap('JET')
	set(gca,'xtick',[],'ytick',[])
	img_str=['small/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_B_sm'];
	print('-dpng', img_str);
	clf
	if ~isempty(A) && ~isempty(B)
	  %%%%% Add gaussfit here (JET) 20101231
	       
	  if strcmp(gauss,'gauss')
	    disp('doing gauss fit')
	    isnana=find(isnan(A));                 
	    A(isnana)=0 ;
	    norma=normfit2d([1:size(A,1)], [1:size(A,2)], A);
            
	    isnanb=find(isnan(B));                 
	    B(isnanb)=0 ;
	    normb=normfit2d([1:size(B,1)], [1:size(B,2)], B);
	    
	    imagesc(x,y,A/norma(1)-B/normb(1))
	    %%%%%%%%%%%
	    
	  else
	    imagesc(x,y,A/max(max(A))-B/max(max(B)))
	  end
	  set(gca,'YDir','normal')
	  %xlim([-xra,xra])
	  %ylim([-xra,xra])
	  axis square
	  caxis([-.1,.1])
	  colormap('JET')
	  set(gca,'xtick',[],'ytick',[])
	  img_str=['small/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_dif_sm'];
	  print('-dpng', img_str);
	end
	clf
      end
    end
  end
  %     disp('Cropping...')
  %     system('mogrify -crop ''720x720+253+67'' small/*')

  disp('Shrinking')
  system('mogrify -resize ''7%x7%'' small/*')
  
end
	
if strcmp(orto,'large') || isempty(orto)
  if strcmp(orto,'large')
    disp('skipping small')
  end

  disp('Large maps')

  %Individual plots

  for rr=1:8
    for cc=1:8
      for ti=1:4
	if strcmp(rotate,'rotate')
	  A=imrotate(px(rr,cc,ti).a_map,theta,'bicubic');
	  B=imrotate(px(rr,cc,ti).b_map,theta,'bicubic');
	else
	  A=px(rr,cc,ti).a_map;
	  B=px(rr,cc,ti).b_map;
	end
	imagesc(x,y,A)
	%xlim([-xra,xra])
	%ylim([-xra,xra])
	set(gca,'YDir','normal')
	axis square
	xlabel('Inches')
	ylabel('Inches')
	colorbar
	caxis([ramin,ra])
	colormap('JET')
	title(['Tile ' int2str(ti) ' Row ' int2str(rr) ' Col ' int2str(cc) ' Pol A'])
	img_str=['full/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_A'];
	print('-dpng', img_str);
	clf
	imagesc(x,y,B)
	set(gca,'YDir','normal')
	%xlim([-xra,xra])
	%ylim([-xra,xra])
	axis square
	xlabel('Inches')
	ylabel('Inches')
	colorbar
	caxis([ramin,ra])
	colormap('JET')
	title(['Tile ' int2str(ti) ' Row ' int2str(rr) ' Col ' int2str(cc) ' Pol B'])
	img_str=['full/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_B'];
	print('-dpng', img_str);
	clf
	if ~isempty(A) && ~isempty(B)
	 
	  %%%%Add gaussfit here (JET) 20101231
              
	  if strcmp(gauss,'gauss')
	    isnana=find(isnan(A));                 
	    A(isnana)=0; 
	    norma=normfit2d([1:size(A,1)], [1:size(A,2)], A);
	    
	    isnanb=find(isnan(B));                 
	    B(isnanb)=0; 
	    normb=normfit2d([1:size(B,1)], [1:size(B,2)], B);
  
	    imagesc(x,y,A/norma(1)-B/normb(1))
	    %%%%%%%%%%%

	  else
	    
	    imagesc(x,y,A/max(max(A))-B/max(max(B)))
	  end 
	  set(gca,'YDir','normal')
	  %xlim([-xra,xra])
	  %ylim([-xra,xra])
	  axis square
	  xlabel('Inches')
	  ylabel('Inches')
	  colorbar
	  caxis([-.1,.1])
	  colormap('JET')
	  title(['Tile ' int2str(ti) ' Row ' int2str(rr) ' Col ' int2str(cc) ', A-B'])
	  img_str=['full/det_row_' int2str(rr) '_col_' int2str(cc) '_tile_' int2str(ti) '_dif'];
	  print('-dpng', img_str);
	end
	clf
      end
    end

    %----------------------------------------
    %Plot dark Squids (JET) 20110102
    dsq=[1,34,67,100,133,166,199,232,265,298,331,364,397,430,463,496];
    for i=1:16
      A=sq(i).map;
      imagesc(x,y,A)
      %xlim([-xra,xra])
      %ylim([-xra,xra])
      set(gca,'YDir','normal')
      axis square
      xlabel('Inches')
      ylabel('Inches')
      colorbar
      title(['Dark Squid, GCP=' int2str(dsq(i))])
      img_str=['full/dark_sq_gcp' int2str(dsq(i))];
      print('-dpng', img_str);
    end
    %----------------------------------------
    
  end
end
    
disp('Complete.')
