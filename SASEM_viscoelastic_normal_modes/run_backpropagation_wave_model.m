% Modal analysis for viscoelastic medal via the semi-analytical spectral
% element method

clear;
%frequency parameters
fmin=1; % minimum frequency
fmax=100; % maximum frequency
df=0.5;%  interval
freqs=fmin:df:fmax;% frequency array

% model/grid parameters
modelfile='backpropagation_wave_model.csv';
model_type=1;% 1 for multilayered models and 2 for gradient models 
mode_type = 2;% 1 for fundamental mode only and 2 for multi modes
output_v=0;  % whether to output the eigenwavefields (0 or 1)

global FC PPW NGLL NGRL;

FC=15.0;
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
%% plot 3D complete modes
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
figure();
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);
% plot positive-propagating modes
index=find(imag(hw)>0);
vc1=vc;vc1(index)=NaN;hw1=hw;hw1(index)=NaN;
hold on;scatter3(reshape(repmat(freqs,size(hw1,1),1),[],1),...
    reshape(-imag(hw1),[],1),reshape(vc1,[],1),2,reshape(-imag(hw1),[],1),'filled');

% plot negative-propagating modes
index=find(imag(hw)<0);
vc1=vc;vc1(index)=NaN;hw1=-hw;hw1(index)=NaN;
hold on;scatter3(reshape(repmat(freqs,size(hw1,1),1),[],1),...
    reshape(-imag(hw1),[],1),reshape(vc1,[],1),3,reshape(-imag(hw1),[],1));
axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');zlabel('Absolute Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

co=colorbar;
set(get(co,'XLabel'),'string','-Im (k) (m^{-1})','Fontsize',8,'FontName','times new roman');
set(co,'Fontsize',8,'FontName','times new roman');
set(co,'Location','EastOutside');
Po=get(co,'Position');Po(1)=0.81;Po(2)=0.18;Po(3)=0.04;Po(4)=0.73;
set(co,'Position',Po);
%% plot 3D positive-propagating modes
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
figure();
vc(imag(hw)>0)=NaN;
hw(imag(hw)>0)=NaN;
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);

hold on;scatter3(reshape(repmat(freqs,size(hw,1),1),[],1),...
    reshape(-imag(hw),[],1),reshape(vc,[],1),2,reshape(-imag(hw),[],1),'filled');
axis([-inf,inf,-inf,inf,0.8*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');zlabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

co=colorbar;
set(get(co,'XLabel'),'string','-Im (k) (m^{-1})','Fontsize',8,'FontName','times new roman');
set(co,'Fontsize',8,'FontName','times new roman');
set(co,'Location','EastOutside');
Po=get(co,'Position');Po(1)=0.81;Po(2)=0.18;Po(3)=0.04;Po(4)=0.73;
set(co,'Position',Po);
%% plot 2D
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
figure();
vc(imag(hw)<0)=NaN;
hw(imag(hw)<0)=NaN;
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);

hold on;plot(freqs,vc,'r.','MarkerSize',2);
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');zlabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

%% plot frequency-wavenumber
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
figure();
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);

hold on;plot(2*pi*freqs./vc,freqs,'r.');
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');zlabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;
