function subplot_grid2(n,m,x,dox,doy)
% subplot_grid2(n,m,x)
%
% call after a stacked subplot setup with subplot_stack to finalize
% axis labels
%
% e.g. subplot_grid(3,2,1)
%      plot(randn(10))
%      subplot_grid2(3,2,1)
%
% Somehow, splitting the call to subplot_grid and subplot_grid2 gives different results
% than calling in the same line. That is,
%
% subplot_grid(5,5,1);plot(1:10);
% subplot_grid2(5,5,1);
% 
% ...does not equal
%
% subplot_grid(5,5,1);plot(1:10);subplot_grid2(5,5,1);
%
% It's best to put a drawnow command after the call to plot, or to do all the plots
% like 
%
% h=subplot_grid(5,5,1);plot(1:10);subplot_grid2(5,5,1);
%
% then go back and do set(gcf,'CurrentAxes',h) before the call to subplot_grid2. 


if ~exist('dox','var') || isempty(dox)
  dox=false;
end

if ~exist('doy','var') || isempty(doy)
  doy=true;
end

% for plots below top row remove uppermost y tick label if it is too
% close to top of panel to avoid it conflicting with bottom y tick
% label on the panel above
ypos=ceil(x/m);
if doy
  if(ypos~=1)
    ytv=get(gca,'YTick');
    yl=ylim;
    
    if strcmp(get(gca,'YDir'),'normal')
      ind1=2;
      ind2=numel(ytv);
    else
      ind1=1;
      ind2=1;
    end
    
    if(abs((yl(ind1)-ytv(ind2))/diff(yl))<0.1)
      % This is messed up, but on a log plot, doing
      % set(gca,'YTickLabel',get(gca,'YTickLabel')) does not work. We have no choice
      % but to remove the tickval entirely for log plots.
      if strcmp(get(gca,'YScale'),'linear')
        ytl=get(gca,'YTickLabel');
        % Can be either cell or character array, depending on matlab version
        if ~iscell(ytl)
          ytl(ind2,:)=' ';
        else
          ytl{ind2} = '';
        end
        set(gca,'YTickLabel',ytl);
      else
        set(gca,'YTick',ytv(1:end-1));
      end
    end
  end
  set(gca,'YTickMode','manual');
end



% for plots to the right of the leftmost column, remove the last x tick label if it is
% too close to the bordering tickmark to the right
xpos = mod(x-1,m)+1;
if dox  
  if(xpos~=m)
    xtv=get(gca,'XTick');
    xl=xlim;

    if strcmp(get(gca,'XDir'),'normal')
      ind1=2;
      ind2=numel(xtv);
    else
      ind1=1;
      ind2=1;
    end
    
    if(abs((xl(ind1)-xtv(ind2))/diff(xl))<0.1)
      if strcmp(get(gca,'XScale'),'linear')
        xtl=get(gca,'XTickLabel');
        % Can be either cell or character array, depending on matlab version
        if ~iscell(xtl)
          xtl(ind2,:)=' ';
        else
          xtl{ind2} = '';
        end
        set(gca,'XTickLabel',xtl);
      else
        set(gca,'XTick',xtv(1:end-1));
      end
    end
  end
  set(gca,'XTickMode','manual');
end


if dox
  if xpos~=1
    set(gca,'YtickLabel',[])
    ylabel([]);
  end
end

if doy
  if ypos~=n
    set(gca,'XtickLabel',[])
    xlabel([]);
  end
end

return
