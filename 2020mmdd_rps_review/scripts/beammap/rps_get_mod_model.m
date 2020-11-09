function model = rps_get_mod_model(param,rot)

model = param(5)*(cosd(2*(rot-param(1)))-(param(2)+1)/(param(2)-1))...
    .*(param(3)*cosd(rot)+param(4)*sind(rot)+1);