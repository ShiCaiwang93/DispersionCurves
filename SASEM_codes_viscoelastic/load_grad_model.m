function gmodel = load_grad_model(modelfile,layer_info)
%load_grad_model: This function import model parameters of gradient models 

% Input parameters
% filename: a csv file that contains four columns 
%       Column 1-depth (m)
%       Column 2-vp (m/s) at different depths
%       Column 3-vs (m/s) at different depths
%       Column 4-density (kg/m^3) at different depths
%       Column 5-Qp at different depths
%       Column 6-Qs at different depths

% layer_info: optional input filename that contains prior information of sublayer interfaces

% We always assume a homogeneous half-space. At the last row, the depth
% should be Inf, and the corrsponding parameters belong to the half-space
% substrate.

% Output parameter
% structure gmodel:
%     vp - P-wave velocities at different depths
%     vs -  S-wave velocities at different depths
%     dns - Densities at different depths
%     Qp - P-wave Q values  at different depths
%     Qs -  S-wave Q values  at different depths
%     depth - Vertical coordinates
%     depth_layer - depths of the interfaces between sublayers
%     vp_sub - P-wave velocity of halfspace
%     vs_sub - S-wave velocity of halfspace
%     dns_sub - Density of halfspace
%     Qp_sub - P-wave Q value of halfspace
%     Qs_sub - S-wave Q value of halfspace

% Copyright 2022 by Caiwang Shi.

Para_temp=csvread(modelfile,1,0);%'gradient_model.csv'
% check file
if Para_temp(1,1)~=0
    error('Depth of gradient model must starts from zero.');
elseif ~isinf(Para_temp(end,1))
    error('Depth of gradient model must ends with inf.');
end

depth=Para_temp(1:end-1,1).';
vp=Para_temp(1:end-1,2).';
vs=Para_temp(1:end-1,3).';
dns=Para_temp(1:end-1,4).';
Qp=Para_temp(1:end-1,5).';
Qs=Para_temp(1:end-1,6).';
%%
if nargin<2
    
    % divide solid model into several sublayers
    range_sub=[5,10];% min and max number of sublayers
    
    vs_solid=vs(vs>0);depth_solid=depth(vs>0);
    gradient=diff(vs_solid)./diff(depth_solid);% gradients of the model
    [B,I] = sort(abs(gradient),'descend');% sort the gradient
    [pks1,locs1] = findpeaks(vs_solid);% local maxima
    [pks2,locs2] = findpeaks(-vs_solid);% local minima
    extremum=[1,sort(union(locs1,locs2),'ascend'),length(vs_solid)];
    if length(extremum)<range_sub(1)+1
        while length(extremum)<range_sub(1)+1
           idx=find(diff(extremum)==max(diff(extremum)),1);
           extremum=[extremum(1:idx),round((extremum(idx)+extremum(idx+1))/2),extremum(idx+1:end)];
        end
    elseif length(extremum)>11
        [B,I]=sort(abs(diff(vs_solid(extremum))),'descend');
        idx=I(range_sub(2):end);idx=idx(idx>1);
        extremum(idx)=[];
        % check extremum
        if extremum(1)~=1
            extremum=[1,extremum];
        elseif extremum(end)~=length(vs_solid)
            extremum(end+1)=length(vs_solid);
        end    
    end
    gmodel.depth_layer=depth_solid(extremum);
    
    % top fluid layer requires no division
    if vs(1)==0
        gmodel.depth_layer=[0,gmodel.depth_layer];
    end
else
    gmodel.depth_layer=load(layer_info);
end


gmodel.depth=depth;
gmodel.vp=vp;
gmodel.vs=vs;
gmodel.Qp=Qp;
gmodel.Qs=Qs;
gmodel.dns=dns;
gmodel.vp_sub=Para_temp(end,2);
gmodel.vs_sub=Para_temp(end,3);
gmodel.dns_sub=Para_temp(end,4);
gmodel.Qp_sub=Para_temp(end,5);
gmodel.Qs_sub=Para_temp(end,6);

% check fluid layers
if std(gmodel.vp(find(gmodel.vs==0)))>1e-6
    error('Only homogeneous fluid layers are supported in this version.');
end
