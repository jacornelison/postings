function model = rps_get_ort_model(param,rs)
% function model = rps_get_ort_model(param,rs)

% reshape to single column by default
if ~exist('rs','var')
    rs = true;
end

dk = [0 45 90 135]+1.25;
ang = [];
az = param(1);
el = param(2);

ang.a0 = param(3);
ang.b0 = param(4);
ang.a1 = param(5);
ang.b1 = param(6);
ang.a2 = param(7);
ang.b2 = param(8);
ang.a3 = param(9);
ang.b3 = param(10);

th = -180:30:180;
rot = th;
freepar.free = [1 1 1 1 1];
freepar.lb = [-360 -1 -1e4 -1e4 -1e4];
freepar.ub = [360 1 1e4 1e4 1e4];

model = NaN(4,4,2);
for i = 1:4 % dks
    dks = dk(i);
    for j = 0:3 % mce's 1thru3 or just 0
        for k = 1:2 % pol A/B
            if k==1
                pol = 'a';
            else
                pol = 'b';
            end
            
            x = ang.([pol sprintf('%i',j)]);

            psi = dks+x;
            
            A = ((cosd(psi)*cosd(el)+sind(psi)*sind(az)*sind(el))*cosd(th)+sind(psi)*cosd(az)*sind(th)).^2;
            N = abs(cosd(el)*cosd(th)).^2+abs(cosd(th)*sind(az)*sind(el)+cosd(az)*sind(th)).^2;
            guess = [psi 0.000 0 0 1/2];
            
            %             % Estimate parameters
            %             [aparam] = matmin('chisq',...
            %                 guess, freepar,'rps_get_mod_model',A,1,th);
            %
            fun = @(p)nansum((A./N-...
                p(5)*(cosd(2*(rot-p(1)))-(p(2)+1)/(p(2)-1)).*...
                (p(3)*cosd(rot)+p(4)*sind(rot)+1)).^2);
            aparam = fminsearch(fun,guess);
            model(i,j+1,k) = aparam(1)-dks;

        end

    end

end

if rs
    model = reshape(model,1,[]);
end
