function mapopt=strip_opt(mapopt)
% mapopt=strip_opt(mapopt)
%
% remove redundant tags in mapopt to prevent it from growing too
% big. mapopt.p, mapopt.ind, and mapopt.simopt gets set to empty apart
% from the first tag of each season.
%
%

n_tags=length(mapopt.tags);
year0=mapopt.tags{1}(1:4);
for i=2:n_tags
  year=mapopt.tags{i}(1:4);
  if strcmp(year,year0)
    mapopt.ind{i}=[];
    mapopt.p{i}=[];
    if(isfield(mapopt,'simopt'))
      mapopt.simopt{i}=[];
    end
    %sprintf('deleting ind,p, simopt for %i %s',i, year)
  else
    year0=year;     
  end
end

return
  
  
  