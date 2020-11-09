function [map ac]=make_map(ac,m,coaddopt)
% map = make_map(ac,m,coaddopt)
%
% This function used to be the last step of reduc_coaddpairmaps.  Now you
% load the output file, which contains a coadded ac, and feed it into this
% function independently.
%
% If ac is a 2-element cell array, then the first element is assumed to be
% a noise sim, and the second is assumed to be a signal sim.  The two are
% combined to make a signal+noise map.i
%
% coaddopt only needed for coaddtype. If not supplied, assumes coaddtype=0.
%
% coaddopt can be cell arrays if maps from different obsevations/experiments
% were coadded. In this case require assume the coaddoption of the first cell

if ~exist('coaddopt','var')
  coaddtype=0;
elseif iscell(coaddopt)
  coaddtype=coaddopt{1}.coaddtype;
else
  coaddtype=coaddopt.coaddtype;
end


if ~iscell(ac)
  % usual map
  map=make_map_iteration(ac,m,coaddtype);
elseif length(ac)==1
  % usual map mistakenly is a cell
  map=make_map_iteration(ac{1},m,coaddtype);
elseif length(ac)>1
  % noise sim
  map=make_map_iteration(ac{1},m,coaddtype);
  for jj=2:length(ac)
    % signal sim
    signal_map=make_map_iteration(ac{jj},m,coaddtype);
    % sum the maps - we just sum signal maps and keep var maps as the noise
    % only sim ones
    for j=1:numel(map)
      map(j).T=map(j).T+signal_map(j).T;
      if(isfield(map,'Q'))
        map(j).Q=map(j).Q+signal_map(j).Q;
        map(j).U=map(j).U+signal_map(j).U;
      end
      if(isfield(map,'Qpsub') & isfield(signal_map,'Qpsub'))
        map(j).Qpsub=map(j).Qpsub+signal_map(j).Qpsub;
        map(j).Upsub=map(j).Upsub+signal_map(j).Upsub;
        map(j).Qgsub=map(j).Qgsub+signal_map(j).Qgsub;
        map(j).Ugsub=map(j).Ugsub+signal_map(j).Ugsub;
      end      
      if(isfield(map,'Qd') & isfield(signal_map,'Qd'))
        for k=1:numel(map(j).Qd)
          map(j).Qd{k}=map(j).Qd{k}+signal_map(j).Qd{k};
          map(j).Ud{k}=map(j).Ud{k}+signal_map(j).Ud{k};
        end
        map(j).Qdsub=map(j).Qdsub+signal_map(j).Qdsub;
        map(j).Udsub=map(j).Udsub+signal_map(j).Udsub;
      end

    end
  end
else
  error('The input ac must be a struct array or N-element cell array.')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map=make_map_iteration(ac,m,coaddtype)

sizes=size(ac);

for i=numel(ac):-1:1 %done in descending order so we can clear ac as we go

  % copy axis ticks
  map(i).x_tic=m.x_tic;
  map(i).y_tic=m.y_tic;

  % Calculate T map
  map(i).T=ac(i).wz./ac(i).wsum;
  map(i).Tvar=ac(i).wwv./ac(i).wsum.^2;
  
  % Calculate T int time maps
  map(i).Titime=ac(i).sitime;
  if isfield(ac(i),'swmax')
    map(i).Twitime=ac(i).switime./ac(i).swmax;
  end

  % calc T poly-sub and gnd-sub maps
  if isfield(ac(i),'wz_psub')
    map(i).Tpsub=ac(i).wz_psub./ac(i).wsum;
    map(i).Tgsub=ac(i).wz_gsub./ac(i).wsum;
  end

  % only make polarization maps if requested
  switch coaddtype
    
    case {0,1,2,5}
      
      % normalize such that sum weight is 1
      ac(i).wcz=ac(i).wcz./ac(i).w;
      ac(i).wsz=ac(i).wsz./ac(i).w;
      
      ac(i).wcc=ac(i).wcc./ac(i).w;
      ac(i).wss=ac(i).wss./ac(i).w;
      ac(i).wcs=ac(i).wcs./ac(i).w;
      
      ac(i).wwccv=ac(i).wwccv./ac(i).w.^2;
      ac(i).wwssv=ac(i).wwssv./ac(i).w.^2;
      ac(i).wwcsv=ac(i).wwcsv./ac(i).w.^2;

      % calc T 2x2 matrix per pixel
      x=ac(i).wcc; y=ac(i).wss; z=ac(i).wcs;
      ac(i).n=1./(x.*y-z.^2);
      % remove regions where only precision limits have prevented matrix
      % from being singular
      ac(i).n(abs(ac(i).n)>1e3)=NaN;   %used to be 10^12
      % simple matrix inverse of [x,z;z,y] = [a,b;b,c] - Ken T
      ac(i).e=ac(i).n.*y; ac(i).f=ac(i).n.*-z; ac(i).g=ac(i).n.*x;
      
      % long form 2x2 * 2x1 product
      map(i).Q=ac(i).e.*ac(i).wcz+ac(i).f.*ac(i).wsz;
      map(i).U=ac(i).f.*ac(i).wcz+ac(i).g.*ac(i).wsz;

      % calc Q/U var maps according to Ken recipe
      map(i).Qvar=ac(i).e.^2.*ac(i).wwccv+ac(i).f.^2.*ac(i).wwssv+...
	  2*ac(i).e.*ac(i).f.*ac(i).wwcsv;
      map(i).Uvar=ac(i).f.^2.*ac(i).wwccv+ac(i).g.^2.*ac(i).wwssv+...
	  2*ac(i).f.*ac(i).g.*ac(i).wwcsv;
      map(i).QUcovar=ac(i).e.*ac(i).f.*ac(i).wwccv+ac(i).e.*ac(i).g.*ac(i).wwcsv+...
	  ac(i).f.*ac(i).f.*ac(i).wwcsv+ac(i).f.*ac(i).g.*ac(i).wwssv;

      % calc Q/U poly-sub and gnd-sub maps
      if isfield(ac(i),'wcz_psub')
        % normalize
        ac(i).wcz_psub=ac(i).wcz_psub./ac(i).w;
        ac(i).wsz_psub=ac(i).wsz_psub./ac(i).w;
        ac(i).wcz_gsub=ac(i).wcz_gsub./ac(i).w;
        ac(i).wsz_gsub=ac(i).wsz_gsub./ac(i).w;
        % calc the maps
        map(i).Qpsub=ac(i).e.*ac(i).wcz_psub+ac(i).f.*ac(i).wsz_psub;
        map(i).Upsub=ac(i).f.*ac(i).wcz_psub+ac(i).g.*ac(i).wsz_psub;
        map(i).Qgsub=ac(i).e.*ac(i).wcz_gsub+ac(i).f.*ac(i).wsz_gsub;
        map(i).Ugsub=ac(i).f.*ac(i).wcz_gsub+ac(i).g.*ac(i).wsz_gsub;
      end
      
      % calc Q/U removed by deprojection
      if(isfield(ac,'wcd'))
        % for each deproj mode
        for k=1:numel(ac(i).wcd)
          % normalize
          ac(i).wcd{k}=ac(i).wcd{k}./ac(i).w;
          ac(i).wsd{k}=ac(i).wsd{k}./ac(i).w;
          % calc the ind maps
          map(i).Qd{k}=ac(i).e.*ac(i).wcd{k}+ac(i).f.*ac(i).wsd{k};
          map(i).Ud{k}=ac(i).f.*ac(i).wcd{k}+ac(i).g.*ac(i).wsd{k};
          % accumulate
          if k==1
            sum_wcd=ac(i).wcd{k};
            sum_wsd=ac(i).wsd{k};
          else
            sum_wcd=ac(i).wcd{k}+sum_wcd;
            sum_wsd=ac(i).wsd{k}+sum_wsd;
          end            
        end
        % make the total maps
        map(i).Qdsub=ac(i).e.*sum_wcd+ac(i).f.*sum_wsd;
        map(i).Udsub=ac(i).f.*sum_wcd+ac(i).g.*sum_wsd;
      end
      
    case 3

      % make ind pair diff maps
      map(i).D=ac(i).wzdiff./ac(i).w;
      map(i).Dvar=ac(i).wwvd./ac(i).w.^2;
      % calc T poly-sub and gnd-sub maps
      if isfield(ac(i),'wz_psub')
        map(i).Dpsub=ac(i).wzd_psub./ac(i).w;
        map(i).Dgsub=ac(i).wzd_gsub./ac(i).w;
      end
      if(isfield(ac,'wd'))
        % for each deproj mode
        for k=1:numel(ac(i).wd)
          % calc the ind maps
          map(i).Dd{k}=ac(i).wd{k}./ac(i).w;
          % accumulate
          if k==1
            sum_wd=ac(i).wd{k};
          else
            sum_wd=ac(i).wd{k}+sum_wd;
          end
        end
        % make the total maps
        map(i).Ddsub=sum_wd./ac(i).w;
      end
      
  end
  
  % Calculate P int time maps
  if isfield(ac,'ditime')
    map(i).Pitime=ac(i).ditime; % (identical to Titime)
  end
  if isfield(ac,'dwmax')
    map(i).Pwitime=ac(i).dwitime./ac(i).dwmax;
  end

  % Clear for memory's sake
  ac(i)=[];
end

map=reshape(map,sizes);

return
