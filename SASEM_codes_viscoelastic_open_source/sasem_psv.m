function [vc,hw,wavefields] = sasem_psv(gmodel,freqs,model_type,mode_type,output_v)
% This function calculates the PSV normal modes of the stratified model

% Input parameters
% structure 'gmodel' for multilayered models:
%     vp - P-wave velocities in all layers (m/s)
%     vs - S-wave velocities in all layers (m/s);for fluid layer, vs=0
%     dns - Densities in all layers (kg/m^3)
%     thk - thicknesses of all finite layers (m)
%     Qp - Q values for P waves
%     Qs - Q value for S waves
% structure 'gmodel' for gradient models:
%     vp - P-wave velocities according to depth (m/s)
%     vs - S-wave velocities according to depth (m/s);for fluid layer, vs=0
%     dns - Densities according to depth (kg/m^3)
%     Qp - Q values for P waves
%     Qs - Q value for S waves
%     depth - Vertical coordinate(m)
%     depth_layer - number of nodes in sublayers
%     vp_sub - P-wave velocity of halfspace
%     vs_sub - S-wave velocity of halfspace
%     dns_sub - Density of halfspace
%     Qp - Qp values of halfspace
%     Qs - Qs value of halfspace
% freqs: frequencies (Hz)
% model_type: 1 for multilayered models and 2 for gradient models
% mode_type : 1 for fundamental mode only and 2 for multi modes
% output_v  : whether to output the eigenwavefields (0 or 1)

% Output parameters
% vc: multi-mode phase velocities
% hw: complex horizontal wavenumber

% Copyright 2022 by Caiwang Shi.
%% 
global FC NGLL NGRL;
[zgll,wgll,hgll] = GetGLL(NGLL);
[zgrl,wgrl,hgrl] = GetGRL(NGRL);

vc=ones(1,length(freqs))*nan;
hw=ones(1,length(freqs))*nan;
wavefields=struct([]);

for freno=1:length(freqs)
    omega=2*pi*freqs(freno);
%   Discretization of the model
    if model_type==1
        [fluid,solid] = psv_discrete_homo(gmodel,freqs(freno),zgll,zgrl);
    elseif model_type==2
        [fluid,solid] = psv_discrete_inhomo(gmodel,freqs(freno),zgll,zgrl);
    end
    if isempty(fluid)
       nwater=0; 
    else
       nwater=fluid.n;
    end
%   Complex wave velocities
%     cVs=solid.Vs.*(1+log(freqs(freno)/FC)./(pi*solid.Qs)+1i./(2*solid.Qs));
%     cVp=solid.Vp.*(1+log(freqs(freno)/FC)./(pi*solid.Qp)+1i./(2*solid.Qp));
    cVs=solid.Vs./(1+log(FC/freqs(freno))/pi./solid.Qs)...
        .*(1+sqrt(1+1./solid.Qs.^2)+1i./solid.Qs)./sqrt(1+1./solid.Qs.^2)/2;
    cVp=solid.Vp./(1+log(FC/freqs(freno))/pi./solid.Qp)...
        .*(1+sqrt(1+1./solid.Qp.^2)+1i./solid.Qp)./sqrt(1+1./solid.Qp.^2)/2;
    
%     cVs_ife=solid.Vs_ife.*(1+log(freqs(freno)/FC)./(pi*solid.Qs_ife)+1i./(2*solid.Qs_ife));
%     cVp_ife=solid.Vp_ife.*(1+log(freqs(freno)/FC)./(pi*solid.Qp_ife)+1i./(2*solid.Qp_ife));
    cVs_ife=solid.Vs_ife./(1+log(FC/freqs(freno))/pi./solid.Qs_ife)...
            .*(1+sqrt(1+1./solid.Qs_ife.^2)+1i./solid.Qs_ife)./sqrt(1+1./solid.Qs_ife.^2)/2;
    cVp_ife=solid.Vp_ife./(1+log(FC/freqs(freno))/pi./solid.Qp_ife)...
            .*(1+sqrt(1+1./solid.Qp_ife.^2)+1i./solid.Qp_ife)./sqrt(1+1./solid.Qp_ife.^2)/2;

%   Matrices for solid layers (to model displacements)
    Mu=solid.Dns.*cVs.^2;
    Lamda=solid.Dns.*(cVp.^2-2*cVs.^2);
    
    ndisp=solid.n+NGRL-1+NGLL-1;
    NES=sum(solid.e);
    MS_1 = zeros(ndisp,1);MS_2 = zeros(ndisp,1);
    KS2_1= zeros(ndisp,1); KS2_2= zeros(ndisp,1); 
    KS1_1 = zeros(NGLL,NGLL,NES+1);KS1_2 = zeros(NGLL,NGLL,NES+1);
    ES_1 = zeros(NGLL,NGLL,NES+1);ES_2 = zeros(NGLL,NGLL,NES+1);
    
    for eid=1:NES
        dze = solid.coor(solid.idx(end,eid))-solid.coor(solid.idx(1,eid));
        dz_dzi = 0.5*dze;

        MS_1(solid.idx(:,eid)) = MS_1(solid.idx(:,eid)) + wgll.*solid.Dns(:,eid)*omega^2*dz_dzi;
        MS_2(solid.idx(:,eid)) = MS_2(solid.idx(:,eid)) + wgll.*solid.Dns(:,eid)*omega^2*dz_dzi;
        KS2_1(solid.idx(:,eid)) = KS2_1(solid.idx(:,eid)) + wgll.*(Lamda(:,eid)+2*Mu(:,eid))*dz_dzi;
        KS2_2(solid.idx(:,eid)) = KS2_2(solid.idx(:,eid)) + wgll.*(Mu(:,eid))*dz_dzi;
        
        W1 = Lamda(:,eid).*wgll;W2 = Mu(:,eid).*wgll;
        KS1_1(:,:,eid) = repmat(W1,1,NGLL).* hgll';
        KS1_2(:,:,eid) = repmat(W2,1,NGLL).* hgll';
        
        W1 = (Lamda(:,eid)+2*Mu(:,eid)).*wgll/dz_dzi;W2 = Mu(:,eid).*wgll/dz_dzi;
        ES_1(:,:,eid) = hgll * ( repmat(W2,1,NGLL).* hgll');
        ES_2(:,:,eid) = hgll * ( repmat(W1,1,NGLL).* hgll');
    end

%   Matrices for buffer layer (to model displacements)
    eid=NES+1;
    dze = solid.coor(solid.idx(end,eid))-solid.coor(solid.idx(1,eid));
    dz_dzi = 0.5*dze;
    
    Mu_IFE=solid.Dns_ife.*cVs_ife.^2;
    Lamda_IFE=solid.Dns_ife.*(cVp_ife.^2-2*cVs_ife.^2);
    
    MS_1(solid.idx(:,eid)) = MS_1(solid.idx(:,eid)) + wgll.*solid.Dns_ife*omega^2*dz_dzi;
    MS_2(solid.idx(:,eid)) = MS_2(solid.idx(:,eid)) + wgll.*solid.Dns_ife*omega^2*dz_dzi;
    KS2_1(solid.idx(:,eid)) = KS2_1(solid.idx(:,eid)) + wgll.*(Lamda_IFE+2*Mu_IFE)*dz_dzi;
    KS2_2(solid.idx(:,eid)) = KS2_2(solid.idx(:,eid)) + wgll.*(Mu_IFE)*dz_dzi;

    W1 = Lamda_IFE.*wgll;W2 = Mu_IFE.*wgll;
    KS1_1(:,:,eid) = repmat(W1,1,NGLL).* hgll';
    KS1_2(:,:,eid) = repmat(W2,1,NGLL).* hgll';

    W1 = (Lamda_IFE+2*Mu_IFE).*wgll/dz_dzi;W2 = Mu_IFE.*wgll/dz_dzi;
    ES_1(:,:,eid) = hgll * ( repmat(W2,1,NGLL).* hgll');
    ES_2(:,:,eid) = hgll * ( repmat(W1,1,NGLL).* hgll');
    
%   Matrices for semi-infinite element (to model displacements)
    dz_dzi=solid.scale;
    MS_1(solid.grl_idx) = MS_1(solid.grl_idx) + wgrl.*solid.Dns_ife*omega^2*dz_dzi;
    MS_2(solid.grl_idx) = MS_2(solid.grl_idx) + wgrl.*solid.Dns_ife*omega^2*dz_dzi;
    KS2_1(solid.grl_idx) = KS2_1(solid.grl_idx) + wgrl.*(Lamda_IFE+2*Mu_IFE)*dz_dzi;
    KS2_2(solid.grl_idx) = KS2_2(solid.grl_idx) + wgrl.*Mu_IFE*dz_dzi;
    
    W1 = Lamda_IFE.*wgrl;W2 = Mu_IFE.*wgrl;
    K1_1_IFE=repmat(W1,1,NGRL).* hgrl';
    K1_2_IFE=repmat(W2,1,NGRL).* hgrl';
    
    W1 = (Lamda_IFE+2*Mu_IFE).*wgrl/dz_dzi;W2 = Mu_IFE.*wgrl/dz_dzi;
    E_1_IFE = hgrl * ( repmat(W2,1,NGRL).* hgrl');
    E_2_IFE = hgrl * ( repmat(W1,1,NGRL).* hgrl');
    
    if isempty(fluid)
        M=diag([MS_1;MS_2]);K2=diag([KS2_1;KS2_2]);
        E = zeros(2*ndisp,2*ndisp);

        for eid=1:NES+1
            idx1 = solid.idx(:,eid);idx2=solid.idx(:,eid)+ndisp;
            K2(idx1,idx2)=K2(idx1,idx2)+KS1_1(:,:,eid)-KS1_2(:,:,eid).';   
            E(idx1,idx1)=E(idx1,idx1)+ES_1(:,:,eid); 
            E(idx2,idx2)=E(idx2,idx2)+ES_2(:,:,eid);     
            E(idx2,idx1)=E(idx2,idx1)+KS1_1(:,:,eid).'-KS1_2(:,:,eid);
        end

        idx1 = solid.grl_idx;idx2=solid.grl_idx+ndisp;
        K2(idx1,idx2)=K2(idx1,idx2)+K1_1_IFE-K1_2_IFE.'; 
        E(idx1,idx1)=E(idx1,idx1)+E_1_IFE; 
        E(idx2,idx2)=E(idx2,idx2)+E_2_IFE;
        E(idx2,idx1)=E(idx2,idx1)+K1_1_IFE.'-K1_2_IFE;

    else
       %   Add matrices for fluid layers
       MW = zeros(fluid.n,1);
       KW = zeros(fluid.n,1);
       EW=zeros(NGLL,NGLL,fluid.e);

       for eid=1:sum(fluid.e)
            dze = fluid.coor(fluid.idx(end,eid))-fluid.coor(fluid.idx(1,eid));
            dz_dzi = 0.5*dze;
            MW(fluid.idx(:,eid)) = MW(fluid.idx(:,eid)) ...
                                      + fluid.Dns(:,eid).*wgll.*omega^2*dz_dzi;
            KW(fluid.idx(:,eid)) = KW(fluid.idx(:,eid))...
                                      + fluid.Dns(:,eid).*fluid.Vp(:,eid).^2.*wgll*dz_dzi;
            W = fluid.Dns(:,eid).*fluid.Vp(:,eid).^2.*wgll/dz_dzi;
            EW(:,:,eid) = hgll * ( repmat(W,1,NGLL).* hgll');
       end
        
        M=diag([MW;MS_1;MS_2]);K2=diag([KW;KS2_1;KS2_2]);
        E = zeros(fluid.n+2*ndisp,fluid.n+2*ndisp);
        
        for eid=1:sum(fluid.e)
            idx = fluid.idx(:,eid);
            E(idx,idx)=E(idx,idx)+EW(:,:,eid); 
        end

        for eid=1:NES+1
            idx1 = fluid.n+solid.idx(:,eid);idx2=fluid.n+solid.idx(:,eid)+ndisp;
            K2(idx1,idx2)=K2(idx1,idx2)+KS1_1(:,:,eid)-KS1_2(:,:,eid).';   
            E(idx1,idx1)=E(idx1,idx1)+ES_1(:,:,eid); 
            E(idx2,idx2)=E(idx2,idx2)+ES_2(:,:,eid);     
            E(idx2,idx1)=E(idx2,idx1)+KS1_1(:,:,eid).'-KS1_2(:,:,eid);
        end

        idx1 = fluid.n+solid.grl_idx;idx2=fluid.n+solid.grl_idx+ndisp;
        K2(idx1,idx2)=K2(idx1,idx2)+K1_1_IFE-K1_2_IFE.'; 
        E(idx1,idx1)=E(idx1,idx1)+E_1_IFE; 
        E(idx2,idx2)=E(idx2,idx2)+E_2_IFE;
        E(idx2,idx1)=E(idx2,idx1)+K1_1_IFE.'-K1_2_IFE;

        
        % add free sea surface condition  
        M(1,1)=0;E(1,:)=0;E(:,1)=0;
        % conective condition between water and solid layers
        idx_water_edge=fluid.n;idx_solid_edge=fluid.n+ndisp+1;
        M(idx_water_edge,idx_solid_edge)=-fluid.Vp(end).^2*fluid.Dns(end);
        M(idx_solid_edge,idx_water_edge)=-fluid.Dns(end)*omega^2;
    end
    
    % eigenvalue decomposition
    [eig_U,D]=eig(K2\(M-E));
    K=sqrt(diag(D));

    % mode filter
    
    vs_range0=[min(min(cVs.*conj(cVs)./real(cVs))),max(max(cVs.*conj(cVs)./real(cVs)))];
    vs_range1=[min(min(cVs_ife.*conj(cVs_ife)./real(cVs_ife))),max(max(cVs_ife.*conj(cVs_ife)./real(cVs_ife)))];
    vr_range(1)=0.85*min(vs_range0(1),vs_range1(1));
    vr_range(2)=max(vs_range0(2),vs_range1(2));
    
    flag=ones(length(D),1);
    tol=1e0;
    for index=1:length(D)
        if tol*abs(real(K(index)))<abs(imag(K(index))) || abs(omega/real(K(index)))<vr_range(1) ...
                                                       || abs(omega/real(K(index)))>vr_range(2)...
                                                       %|| abs(imag(K(index)))>tol*abs(real(K(index)))
           flag(index)=0; continue;
        end
    end
    
    index=find(flag);
    hw_freq=K(index);
    eig_U=eig_U(:,index);
    vc_freq=real(omega)./real(hw_freq);
    
    % sorting using real-domain phase velocity and wavenumber
    [~,I] = sort(real(hw_freq),'descend');
    hw_freq=hw_freq(I);
    vc_freq=real(omega)./real(hw_freq);
    eig_U=eig_U(:,I);
    
%     ur=V(1:ndisp,index)./hw_freq.';
%     uz=V(ndisp+1:end,index);

    if length(vc_freq)>size(vc,1)
        vc(length(vc_freq),:)=nan;
        hw(length(vc_freq),:)=nan;
    end

    vc(1:length(vc_freq),freno)=vc_freq;
    hw(1:length(hw_freq),freno)=hw_freq;
    
    if output_v==0
       wavefields=[];continue; 
    end
    
    % output wavefields
    wavefields(freno).ur=eig_U(nwater+(1:ndisp),:)./hw_freq.';
    wavefields(freno).uz=eig_U(nwater+ndisp+(1:ndisp),:);
    wavefields(freno).solid_coor=[solid.coor(1:end-1);solid.coor_s_grl];
    if ~isempty(fluid)
        wavefields(freno).p=eig_U(1:nwater,:);
        wavefields(freno).fluid_coor=fluid.coor;
    end
end
return;

