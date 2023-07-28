function dir2html()

%%
figdir = 'c:/Users/James/Documents/GitHub/postings/';
d = dir(figdir);

dates = [];
for dirind = 1:length(d)
di = d(dirind).date;
dates(dirind) = datenum(di,'dd-mmm-yyyy HH:MM:SS');
end


[s, si] = sort(dates);
si = si(end:-1:1);
%
fileID = fopen(fullfile(figdir,'dirlist.html'),'w');

for dirind = 1:length(d)
    di = d(si(dirind)).name;
    dt = d(si(dirind)).date(1:11);
    if length(di)>4 & ~isempty(str2num(di(1:4))) %~strcmp(di(1),'.')
        idxname = fullfile(figdir,di,'index.html');
        if ~exist(idxname,'file')
            continue
        end
        idxID = fopen(idxname);
        tline = fgetl(idxID);
        while ischar(tline)
            if strfind(tline,'<h1>')
                ttl = strtrim(tline);
                ttl = strrep(ttl,'<h1>','');
                ttl = strrep(ttl,'</h1>','');
                break
            else
                tline = fgetl(idxID);
            end
        end

        txt = sprintf('<p><a href=''%s/index.html''>%s: %s</a></p>\n',di,dt,ttl);
        fprintf(fileID,txt,d(dirind).name);
    end
end

fclose(fileID);