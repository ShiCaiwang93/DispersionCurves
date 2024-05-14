function [fluid,solid] = psv_discrete_homo(gmodel,freq,zgll,zgrl)
% This function discretize the multi-layered homogeneous model into
% numerous element controlled by Gauss-Lobatto-Legendre and Gauss-Radau-Laguerre nodes.

% Input parameters
% structure gmodel:
%     vp - P-wave velocities in all layers (m/s)
%     vs - S-wave velocities in all layers (m/s);for fluid layer, vs=0
%     dns - Densities in all layers (kg/m^3)
%     thk - thicknesses of all finite layers (m)
%     Qp - Q values for P waves
%     Qs - Q value for S waves
% freq: frequency (Hz)
% mode_type : 1 for fundamental mode only and 2 for multi modes

% zgll coordinates of standard GLL nodes in the range of [-1,1]
% zgrl coordinates of standard GRL nodes in the range of [0,inf]

% Output parameters: two structure variables: fluid and solid.
% fluid: element information for fluid layers ,including
%        efluid - number of fluid elements in fluid layers
%        idx - indexes of fluid nodes
%        nfluid - number of fluid nodes
%        coor - coordinates of fluid nodes
%        Vp - P-wave velocities of fluid nodes
%        Dns - densities of fluid nodes
%        Qp - P-wave Q values velocities of fluid nodes

% solid: element information for solid layers ,including
%        esolid - number of solid elements in finite layers
%        idx - indexes of solid nodes in finite layers (including the buffer layer)
%        nsolid - number of solid nodes in finite layers
%        coor - coordinates of solid nodes in finite layers
%        vp - P-wave velocities of solid nodes in finite layers
%        vs - S-wave velocities of solid nodes in finite layers
%        Dns - Density of solid nodes in finite layers
%        Qp - P-wave Q values of solid nodes in finite layers
%        Qs - S-wave Q values of solid nodes in finite layers
%        grl_idx - indexes of grl nodes
%        vp_ife - P-wave velocities of solid nodes in infinite layers
%        vs_ife - S-wave velocities of solid nodes in infinite layers
%        Dns_ife - Dns of solid nodes in infinite layers
%        Qp_ife - Qp of solid nodes in infinite layers
%        Qs_ife - Qs of solid nodes in infinite layers
%        coor_p_grl - P-wave grl-node coordinates of solid nodes in infinite layers
%        coor_s_grl - S-wave grl-node coordinates of solid nodes in infinite layers
%        scale factor from stand grl scale to physical scale

% Copyright 2022 by Caiwang Shi.

vp=gmodel.vp;
vs=gmodel.vs;
dns=gmodel.dns;
thk=gmodel.thk;
qp=gmodel.Qp;
qs=gmodel.Qs;

% check fluid layers
fluid_idx=find(vs==0);
nlayer_fluid=length(fluid_idx);
if (nlayer_fluid>0 && min(fluid_idx) >1) || (nlayer_fluid>1 && max(diff(fluid_idx)) >1)
   error('Fluid layers should appear on the top of the model.');
elseif nlayer_fluid>=length(thk)
   error('Pure fluid model is not supported.');
end

solid_idx=(nlayer_fluid+1):length(vp);
vp_solid=vp(solid_idx);
vs_solid=vs(solid_idx);
dns_solid=dns(solid_idx);
qp_solid=qp(solid_idx);
qs_solid=qs(solid_idx);
thk_solid=thk(solid_idx(1):end);

% global parameters
global PPW NGLL NGRL;

% maximum and minimum wavelength
wlength_finite_min=vs_solid(1:end-1)/freq;
wlength_infinite_max=vp_solid(end)/freq;
wlength_infinite_min=vs_solid(end)/freq;

% Number of GLL elements used to discretize finite solid layers 
NEL_FEM=ceil(thk_solid./wlength_finite_min*(PPW-1)/(NGLL-1));
NEL_cumsum=cumsum(NEL_FEM);

% Coordinates of the top/bottom interfaces in GLL elements
Z_FEM=zeros(sum(NEL_FEM)+1,1);
depth=cumsum(thk_solid);
for layer_idx=1:length(solid_idx)-1
    if layer_idx==1
        Z_FEM(1:(NEL_cumsum(1)+1))=linspace(0,depth(layer_idx),NEL_FEM(layer_idx)+1);
    else
        Z_FEM((NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)+1))=...
            linspace(depth(layer_idx-1),depth(layer_idx),NEL_FEM(layer_idx)+1);
    end
end

% Nodes for wavefields
NEL=sum(NEL_FEM);
nfield = NEL*(NGLL-1)+1; % number of nodes	
ifield = zeros(NGLL,NEL+1);% global indexes of these 
coor_solid = zeros(nfield,1);% coordinates of GLL nodes		
for  e=1:NEL
    dze = Z_FEM(e+1)-Z_FEM(e);
    ifield(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';
    coor_solid(ifield(:,e)) = 0.5*(Z_FEM(e)+Z_FEM(e+1)) + 0.5*dze*zgll;
end

% parameters at each nodes in finite solid layers
Dns = zeros(NGLL,NEL);
Vs = zeros(NGLL,NEL);
Vp = zeros(NGLL,NEL);
Qs = zeros(NGLL,NEL);
Qp = zeros(NGLL,NEL);
for layer_idx=1:length(solid_idx)-1
    if layer_idx==1
        Dns(:,1:NEL_cumsum(1))=dns_solid(1); 
        Vs(:,1:NEL_cumsum(1))=vs_solid(1); 
        Vp(:,1:NEL_cumsum(1))=vp_solid(1);
        Qs(:,1:NEL_cumsum(1))=qs_solid(1); 
        Qp(:,1:NEL_cumsum(1))=qp_solid(1);
    else
        Dns(:,(NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)))=dns_solid(layer_idx); 
        Vs(:,(NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)))=vs_solid(layer_idx); 
        Vp(:,(NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)))=vp_solid(layer_idx);
        Qs(:,(NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)))=qs_solid(layer_idx); 
        Qp(:,(NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)))=qp_solid(layer_idx);
    end
end

% Discretization of the half-space
    
% a buffer layer using GLL element (to model displacements wavefields)
e=NEL+1;
dze=Z_FEM(end)-Z_FEM(end-1);%0.1*wlength_infinite_min;
ifield(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';
coor_solid(ifield(:,e)) = (Z_FEM(e)+0.5*dze) + 0.5*dze*zgll;

% semi-infinite element
e=NEL+2;
scale=(10*wlength_infinite_min/max(zgrl));
ifield_ife = (e-1)*(NGLL-1)+(1:NGRL)';
coor_p_grl = Z_FEM(NEL+1) + dze + scale*zgrl;
coor_s_grl = Z_FEM(NEL+1) + dze + scale*zgrl;

Dns_ife=dns_solid(end);
Vs_ife=vs_solid(end);
Vp_ife=vp_solid(end);
Qs_ife=qs_solid(end);
Qp_ife=qp_solid(end);

% add fluid layers
fluid=[];
if nlayer_fluid>0
    %wavelength of P waves in water
    wlength_water=vp(fluid_idx)/freq;
    % Number of GLL elements used to discretize finite solid layers 
    NEL_FEM_fluid=2*ceil(thk(fluid_idx)./wlength_water*(PPW-1)/(NGLL-1));
    NEL_fluid_cumsum=cumsum(NEL_FEM_fluid);

    % Coordinates of the top/bottom interfaces in GLL elements
    Z_FEM_fluid=zeros(sum(NEL_FEM_fluid)+1,1);
    depth_fluid=cumsum(thk(fluid_idx));
    for layer_idx=1:length(fluid_idx)
        if layer_idx==1
            Z_FEM_fluid(1:(NEL_fluid_cumsum(1)+1))=linspace(0,depth_fluid(layer_idx),NEL_FEM_fluid(layer_idx)+1);
        else
            Z_FEM_fluid((NEL_fluid_cumsum(layer_idx-1)+1):(NEL_fluid_cumsum(layer_idx)+1))=...
                linspace(depth_fluid(layer_idx-1),depth_fluid(layer_idx),NEL_FEM_fluid(layer_idx)+1);
        end
    end

    % Nodes for wavefields
    NEL_fluid=sum(NEL_FEM_fluid);
    nfluid = NEL_fluid*(NGLL-1)+1; % number of nodes	
    ifluid = zeros(NGLL,NEL_fluid);% global indexes of these 
    coor_fluid = zeros(nfluid,1);% coordinates of GLL nodes		
    for  e=1:NEL_fluid
        dze = Z_FEM_fluid(e+1)-Z_FEM_fluid(e);
        ifluid(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';
        coor_fluid(ifluid(:,e)) = 0.5*(Z_FEM_fluid(e)+Z_FEM_fluid(e+1)) + 0.5*dze*zgll;
    end

    % parameters at each nodes in finite solid layers
    Dns_fluid = zeros(NGLL,NEL_fluid);
    Vp_fluid = zeros(NGLL,NEL_fluid);
    Qp_fluid = zeros(NGLL,NEL_fluid);
    for layer_idx=1:length(fluid_idx)
        if layer_idx==1
            Dns_fluid(:,1:NEL_fluid_cumsum(1))=dns(1); 
            Vp_fluid(:,1:NEL_fluid_cumsum(1))=vp(1);
            Qp_fluid(:,1:NEL_fluid_cumsum(1))=qp(1);
        else
            Dns_fluid(:,(NEL_fluid_cumsum(layer_idx-1)+1):(NEL_fluid_cumsum(layer_idx)))=dns(layer_idx); 
            Vp_fluid(:,(NEL_fluid_cumsum(layer_idx-1)+1):(NEL_fluid_cumsum(layer_idx)))=vp(layer_idx);
            Qp_fluid(:,(NEL_fluid_cumsum(layer_idx-1)+1):(NEL_fluid_cumsum(layer_idx)))=qp(layer_idx);
        end
    end   
    
    % generate fluid-related structural variable
    fluid.e=NEL_FEM_fluid;
    fluid.idx=ifluid;
    fluid.n=nfluid;
    fluid.coor=coor_fluid;
    fluid.Vp=Vp_fluid;
    fluid.Dns=Dns_fluid;
    fluid.Qp=Qp_fluid;
    
    % modify the solid informations
    coor_solid=coor_solid+coor_fluid(end);
    coor_p_grl=coor_p_grl+coor_fluid(end);
    coor_s_grl=coor_s_grl+coor_fluid(end);
end

% generate solid-related structural variable
solid.e=NEL_FEM;
solid.idx=ifield;
solid.n=nfield;
solid.coor=coor_solid;
solid.Vp=Vp;
solid.Vs=Vs;
solid.Dns=Dns;
solid.Qp=Qp;
solid.Qs=Qs;
solid.grl_idx=ifield_ife;
solid.Vp_ife=Vp_ife;
solid.Vs_ife=Vs_ife;
solid.Dns_ife=Dns_ife;
solid.Qp_ife=Qp_ife;
solid.Qs_ife=Qs_ife;
solid.coor_p_grl=coor_p_grl;
solid.coor_s_grl=coor_s_grl;
solid.scale=scale;
return;