function gifsicle(tagz,mname)
% gifsicle(tagz,mname)
%
% compose gifsicle command line
%
% gifsicle(get_tags,'maps_h5.25d50_dk-3.gif')

cmd='gifsicle -l -d100';
for i=1:length(tagz)
  cmd=[cmd,' ',tagz{i},'/',mname];
end
system([cmd,' > ',mname]);
  
return
