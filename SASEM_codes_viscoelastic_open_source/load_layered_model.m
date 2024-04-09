function gmodel = load_layered_model(modelfile)
%load_layered_model: This function import model parameters of multilayered models 

% Input parameters
% filename: a csv file that contains four columns 
%       Column 1-thickness (m)
%       Column 2-vp (m/s) of different layers
%       Column 3-vs (m/s) of different layers
%       Column 4-density (kg/m^3) of different layers
%       Column 5-Qp of different layers
%       Column 6-Qs of different layers

% We always assume a homogeneous half-space. At the last row, the depth
% should be Inf, and the corrsponding parameters belong to the half-space
% substrate.

% Output parameter
% structure gmodel:
%     vp - P-wave velocities of different layers
%     vs -  S-wave velocities of different layers
%     dns - Densities of different layers
%     thk - thickness
%     Qp - P-wave Q values of different layers
%     Qs - S-wave Q values  of different layers

% Copyright 2022 by Caiwang Shi.

Para_temp=csvread(modelfile,1,0);%'layered_model.csv'
% check file
if ~isinf(Para_temp(end,1))
    error('Thickness of layered model must ends with inf.');
end

gmodel.thk=Para_temp(1:end-1,1).';
gmodel.vp=Para_temp(:,2).';
gmodel.vs=Para_temp(:,3).';
gmodel.dns=Para_temp(:,4).';
gmodel.Qp=Para_temp(:,5).';
gmodel.Qs=Para_temp(:,6).';

% check fluid layers
if std(gmodel.vp(find(gmodel.vs==0)))>1e-6
    error('Only homogeneous fluid layers are supported in this version.');
end




