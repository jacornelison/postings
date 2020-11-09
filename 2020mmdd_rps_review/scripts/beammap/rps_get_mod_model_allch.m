function model = rps_get_mod_model_allch(param,eta,cdata,source,rs)
% function model = rps_get_mod_model_allch(param,eta,cdata,source,rs)
% 

if ~exist('rs','var')
    rs = true;
end

chans = unique(cdata.ch);

% Figure out where all of the parameters are
lp = length(chans);
lx = lp;
ln = 2;
la = 2;
lg = size(cdata.data,1);

model = NaN(length(cdata),length(eta));
for i = 1:lg
    % find the parameters specific to this mod curve:
    ind = find(ismember(chans,cdata.ch(i)));
    p = param(ind);
    x = param(lp + ind);
    n1 = param(lp + lx + 1);
    n2 = param(lp + lx + 2);
    a1 = param(lp + lx + ln + 1);
    a2 = param(lp + lx + ln + 2);
    g = param(lp + lx + ln + la + i);
    
    aparam = [p x n1 n2 a1 a2 g];
    
    A = cdata.A(i,:);
    B3 = cdata.B3(i,:);
    B1 = cdata.B1(i,:);
    
    model(i,:) = rps_get_mod_model_vectors(aparam,eta,A,B3,B1,source);
    
end

if rs
    model = reshape(model,1,[]);
end
