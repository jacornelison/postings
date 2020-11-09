function r=get_run_info(str,csvfile)
% r=get_run_info(str)
%
% Get a structure containing info on all runs
%
% Optional input str will cause cutdown list of only those runs
% containing string str
%
% e.g: r=get_run_info('20100103a')

if (nargin<2)||isempty(csvfile)
  csvfile='aux_data/tag_list.csv';
end

r=[];

% Any runs to be added by hand go here
% r=ae(r,'04-jan-2010:01:20:16','04-jan-2010:03:15:30','dummy','20100101a','dummy');
% r=ae(r,'03-jan-2010:01:20:16','03-jan-2010:03:15:30','9_CMB_30_1B_006.sch','20100103a','cmb');
% r=ae(r,'03-jan-2010:01:20:16','03-jan-2010:05:08:21','9_CMB_30_1B_006.sch','20100103a','cmb');
% r=ae(r,'06-feb-2010:08:46:52','06-feb-2010:21:33:51','1_ffflat_moon_raster_07_001.sch','20100206a','moon'); % moon      
% r=ae(r,'06-feb-2010:08:46:52','06-feb-2010:09:46:52','1_ffflat_moon_raster_07_001.sch','20100206dk180','moon'); % moon dk=-180


if exist(csvfile,'file')
  [p k]=ParameterRead(csvfile);
  if exist('str','var')
    ind=false(size(p.tag));
    if iscell(str)
      for i=1:length(str)
        nchar=length(str{i});
        indx=strncmpi(p.tag,str{i},nchar);
        ind=indx | ind;
      end 
    else
      nchar = length(str);
      ind=strncmpi(p.tag,str,nchar);
    end
    p=structcut(p,ind);
  end
  if isfield(p,'tag')
    r.tstart=p.tstart;
    r.tend=p.tend;
    r.schedule=p.sch;
    r.tag=p.tag;
    r.type=p.type;
  end
else
  disp (['No tag list file ' csvfile ' found!']);
end

return