function [data,x,y,phi,az,el,dk,quad] = rps_concat_data(tods,ch)

[data,x,y,phi,az,el,dk,quad] = deal([]);

for i=1:length(tods)
    ch_ind = find(tods{i}.ch == ch);
    data = [data; tods{i}.todcos(:,ch_ind)];
    quad = [quad; tods{i}.todsin(:,ch_ind)];
    x = [x; tods{i}.x];
    y = [y; tods{i}.y];
    phi = [phi; tods{i}.phi];
    az = [az; tods{i}.az];
    el = [el; tods{i}.el];
    dk = [dk; tods{i}.dk];
    
end
