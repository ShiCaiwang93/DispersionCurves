% Modal analysis for viscoelastic medal via the semi-analytical spectral
% element method

clear;
%frequency parameters
fmin=0.51; % minimum frequency
fmax=3; % maximum frequency
df=0.03;%  interval
freqs=fmin:df:fmax;% frequency array

% model/grid parameters
modelfile='offshore_gradient_model.csv';
model_type= 2;% 1 for multilayered models and 2 for gradient models 
mode_type = 2;% 1 for fundamental mode only and 2 for multi modes
output_v=0;  % whether to output the eigenwavefields (0 or 1)

global FC PPW NGLL NGRL;

FC=1.0;
PPW = 8; % polynomial degree of GLL elements (should be greater than 5; ...
       % higher degree results in higher accuracy and more calculation)
NGLL = 8; % number of GLL nodes per finite element
NGRL = 20; % number of GRL nodes in the infinite element (between 10 to 20)

% forward modeling
if model_type==1
    gmodel=load_layered_model(modelfile);
elseif model_type==2
    gmodel=load_grad_model(modelfile);
end
tic;
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
toc;
%% plot
figure();
set(gcf,'unit','centimeters','position',[10,10,14,6]);
subplot(1,4,1);set(gca,'position',[0.06 0.18 0.13 0.73])
h1=plot(gmodel.vp,[gmodel.depth(1:end-1),800],'b');
hold on;h2=plot(gmodel.vs,[gmodel.depth(1:end-1),800],'r');
hold on;plot([0,1e4],[50,50],'k:');hold on;plot([0,1e4],[730,730],'k:')
legend([h1,h2],'\alpha','\beta');
axis([0,3500,0,inf]);box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
set(gca,'ytick',[0:200:800],'ydir','reverse','xminortick','on')
ylabel('Depth (m)');xlabel('Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);
set(gca,'ydir','reverse','xminortick','on')
title ('(a)')

subplot(1,4,2);set(gca,'position',[0.21 0.18 0.13 0.73])
h1=plot(gmodel.Qp,[gmodel.depth(1:end-1),800],'b');
hold on;h2=plot(gmodel.Qs,[gmodel.depth(1:end-1),800],'r');
hold on;plot([0,1e4],[50,50],'k:');hold on;plot([0,1e4],[730,730],'k:');
legend([h1,h2],'Qp','Qs');
axis([0,400,0,inf]);
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Q');
set(gca,'fontname','times new roman','fontsize',8);
set(gca,'xtick',[100,300],'ytick',[0:200:800],'YTickLabel',[],'ydir','reverse','xminortick','on')

subplot(1,4,3);set(gca,'position',[0.36 0.18 0.13 0.73])
plot(gmodel.dns,[gmodel.depth(1:end-1),800],'b')
hold on;plot([0,1e4],[50,50],'k:');hold on;plot([0,1e4],[730,730],'k:')
axis([900,2500,0,inf]);
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Density (kg/m^3)');
set(gca,'fontname','times new roman','fontsize',8);
set(gca,'xtick',[1000,2000],'ytick',[0:200:800],'YTickLabel',[],'ydir','reverse','xminortick','on')

subplot(1,4,4);set(gca,'position',[0.59 0.18 0.305 0.73])
% hold on;plot(freqs,vc,'k.','markersize',3);axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
hold on;scatter(reshape(repmat(freqs,size(vc,1),1),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
axis([-inf,inf,0.9*min(gmodel.vs(find(gmodel.vs,1,'first'))),1.05*max(gmodel.vs)]);colormap(jet);
% hold on;scatter3(reshape(repmat(freqs,size(hw,1),1),[],1),...
%     reshape(-imag(hw),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
% axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

co=colorbar;
set(get(co,'XLabel'),'string','-Im (k) (m^{-1})','Fontsize',8,'FontName','times new roman');
set(co,'Fontsize',8,'FontName','times new roman');
set(co,'Location','EastOutside');
Po=get(co,'Position');Po(1)=0.91;Po(2)=0.18;Po(3)=0.02;Po(4)=0.73;
set(co,'Position',Po);
title('(b)')