function [model] = rps_get_model(param,tods,ch,rpsopt)
% rps_get_model(param,tods,ch,rpsopt)
% ---
% Construct modeled 2D Gaussian maps from input params.
% then reshape and concat into a long Nx1 array.

model = [];

% Loops through source angles and
for i=1:length(tods)
    % Start rotation at zero
    rot = tods{i}.scan.abs_rot-rpsopt.grid.horizontal;
    ch_ind = find(tods{i}.ch == ch);
    
    % Offset beam center.
    dx = tods{i}.x - param(1);
    dy = tods{i}.y - param(2);
    
    % Covariance matrix for beam ellipticity.
    cov = [param(3)^2, param(5) * param(3) * param(4); ...
        param(5) * param(3) * param(4), param(4)^2];
    
    % The inverse function really hates if there's NaN's in the covariance matrix.
    % Instead of doing anything productive about this, I've decided to simply suppress
    % the output for now. Apologies to the reader.
    [T,icov] = evalc('inv(cov)');
    
    % Elliptical Gaussian model (Vector)
    g = exp(-0.5 * ...
        (dx.^2 * icov(1,1) + dy.^2 * icov(2,2) + 2 * dx .* dy * icov(1,2)));
    
    % Amplitude model (Scalar)
    %A = param(11)*(cosd(2*(rot+param(7)))-(param(8)+1)/(param(8)-1))*(param(9)*cosd(rot+param(10))+1);
    A = param(5+i); %Fitting beam first, modulation second.
    % Complete model
    m = A*g;% + param(6); param(6) used to be background.
    
    model = [model; m];
end
