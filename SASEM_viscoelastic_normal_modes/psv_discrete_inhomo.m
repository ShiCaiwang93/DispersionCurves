function [fluid,solid] = psv_discrete_inhomo(gmodel,freq,zgll,zgrl)
% This function discretize the multi-layered homogeneous model into
% numerous element controlled by Gauss-Lobatto-Legendre and Gauss-Radau-Laguerre nodes.

% Input parameters
% structure gmodel:
%     vp - P-wave velocities according to depth (km/s)
%     vs - S-wave velocities according to depth (km/s);for fluid layer, vs=0
%     dns - Densities according to depth (g/cm^3)
%     Qp - P-wave Q values  at different depths
%     Qs -  S-wave Q values  at different depths
%     depth - Vertical coordinate(km)
%     depth_layer - number of nodes in sublayers
%     vp_sub - P-wave velocity of halfspace
%     vs_sub - S-wave velocity of halfspace
%     dns_sub - Density of halfspace
%     Qp_sub - P-wave Q value of halfspace
%     Qs_sub - S-wave Q value of halfspace
% freq: frequency (Hz)
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

% Copyright 2022 by Caiwang Shi .

vp=gmodel.vp;
vs=gmodel.vs;
dns=gmodel.dns;
depth=gmodel.depth;
qp=gmodel.Qp;
qs=gmodel.Qs;
depth_layer=gmodel.depth_layer;
vp_sub=gmodel.vp_sub;
vs_sub=gmodel.vs_sub;
dns_sub=gmodel.dns_sub;
qp_sub=gmodel.Qp_sub;
qs_sub=gmodel.Qs_sub;

% check fluid layers
fluid_idx=find(vs==0);
n_fluid=length(fluid_idx);
if (n_fluid>0 && min(fluid_idx) >1) || (n_fluid>1 && max(diff(fluid_idx)) >1)
   error('Fluid layers should appear on the top of the model.');
elseif n_fluid>=length(depth)
   error('Pure fluid model is not supported.');
end
fluid_idx=1;

if vs(1)==0
    solid_idx=2:length(depth_layer);
    sea_depth=depth_layer(2);
else
    solid_idx=1:length(depth_layer);
    sea_depth=0;
end
thk=diff(depth_layer);
thk_solid=thk(solid_idx(1):end);
% global parameters
global PPW NGLL NGRL;

% maximum and minimum velocity in solid sublayers
vs_solid_min=zeros(1,length(depth_layer));
vp_solid_max=zeros(1,length(depth_layer));

for no_layer=1:length(depth_layer)
   if no_layer==length(depth_layer)
       vs_solid_min(no_layer)=vs_sub;
       vp_solid_max(no_layer)=vp_sub;
   else
       vs_solid_min(no_layer)=min(vs(find(depth>=depth_layer(no_layer) & depth<depth_layer(no_layer+1))));
       vp_solid_max(no_layer)=max(vp(find(depth>=depth_layer(no_layer) & depth<depth_layer(no_layer+1))));
   end
end
if vs(1)==0
    vs_solid_min(1)=[];vp_solid_max(1)=[];
end
% maximum and minimum wavelength
wlength_finite_min=vs_solid_min(1:end-1)/real(freq);
wlength_infinite_max=vp_sub/real(freq);

% Number of GLL elements used to discretize finite solid layers 
NEL_FEM=ceil(thk_solid./wlength_finite_min*(PPW-1)/(NGLL-1));
NEL_cumsum=cumsum(NEL_FEM);

% Coordinates of the top/bottom interfaces in GLL elements
Z_FEM=zeros(sum(NEL_FEM)+1,1);
interface_solid=cumsum(thk_solid);
for layer_idx=1:length(solid_idx)-1
    if layer_idx==1
        Z_FEM(1:(NEL_cumsum(1)+1))=linspace(0,interface_solid(layer_idx),NEL_FEM(layer_idx)+1);
    else
        Z_FEM((NEL_cumsum(layer_idx-1)+1):(NEL_cumsum(layer_idx)+1))=...
            linspace(interface_solid(layer_idx-1),interface_solid(layer_idx),NEL_FEM(layer_idx)+1);
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
for e=1:NEL
    Vp(:,e)=interp1(depth-1e-30,vp,sea_depth+coor_solid(ifield(:,e)),'linear');
    Vs(:,e)=interp1(depth-1e-30,vs,sea_depth+coor_solid(ifield(:,e)),'linear');
    Dns(:,e)=interp1(depth-1e-30,dns,sea_depth+coor_solid(ifield(:,e)),'linear');
    Qp(:,e)=interp1(depth-1e-30,qp,sea_depth+coor_solid(ifield(:,e)),'linear');
    Qs(:,e)=interp1(depth-1e-30,qs,sea_depth+coor_solid(ifield(:,e)),'linear');
end
index=find(isnan(Vp(:,NEL)));
if ~isempty(index)
    Vp(index,NEL)=interp1(depth-1e-30,vp,sea_depth+coor_solid(ifield(:,e)),'previous');
    Vs(index,NEL)=interp1(depth-1e-30,vs,sea_depth+coor_solid(ifield(:,e)),'previous');
    Dns(index,NEL)=interp1(depth-1e-30,dns,sea_depth+coor_solid(ifield(:,e)),'previous');
    Qp(index,NEL)=interp1(depth-1e-30,qp,sea_depth+coor_solid(ifield(:,e)),'previous');
    Qs(index,NEL)=interp1(depth-1e-30,qs,sea_depth+coor_solid(ifield(:,e)),'previous');
end
% Discretization of the half-space
    
% a buffer layer using GLL element (to model displacements wavefields)
e=NEL+1;
dze=Z_FEM(end)-Z_FEM(end-1);%0.2*vs_sub/real(freq);
ifield(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';
coor_solid(ifield(:,e)) = (Z_FEM(e)+0.5*dze) + 0.5*dze*zgll;

% semi-infinite element
e=NEL+2;
scale=(5*wlength_infinite_max/max(zgrl));
ifield_ife = (e-1)*(NGLL-1)+(1:NGRL)';
coor_p_grl = Z_FEM(NEL+1) + dze + scale*zgrl;
coor_s_grl = Z_FEM(NEL+1) + dze + scale*zgrl;

Dns_ife=dns_sub;
Vs_ife=vs_sub;
Vp_ife=vp_sub;
Qs_ife=qs_sub;
Qp_ife=qp_sub;

% add fluid layers
fluid=[];
if n_fluid>0
    %wavelength of P waves in water
    wlength_water=vp(fluid_idx(1))/real(freq);
    % Number of GLL elements used to discretize finite solid layers 
    NEL_FEM_fluid=2*ceil(thk(1)./wlength_water*(PPW-1)/(NGLL-1));
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