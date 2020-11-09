function gof=gof_spaceangle(tm,modfunc,obs,ide)
% gof=gof_spaceangle(tm,modfunc,obs,ide)
%
% Goodness of fit function used to compare online oaz,oel points to
% those derived from ideal iaz,iel using modfunc(tm,iaz,iel)

mva=feval(modfunc,tm,ide);

%sa=spaceangle(obs.az,obs.el,mva.az,mva.el,'deg');
%
%
% JMK: the above is WRONG when the encoder offsets are not small
%  Fixed this by subtracting encoder offsets from both obs and mva coords
%  coordinates, to make them close to ide coords.
sa=spaceangle(obs.az-tm(10),obs.el-tm(11),mva.az-tm(10),mva.el-tm(11),'deg');


gof=sum(sa.^2)/0.1;
