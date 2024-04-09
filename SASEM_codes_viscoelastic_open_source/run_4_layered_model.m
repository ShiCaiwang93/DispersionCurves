% Modal analysis for viscoelastic medal via the semi-analytical spectral
% element method

clear;
%frequency parameters
fmin=1; % minimum frequency
fmax=100; % maximum frequency
df=1;%  interval
freqs=fmin:df:fmax;% frequency array

% model/grid parameters
modelfile='4_layered_model.csv';
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
%% plot
load Muller_4_layer_dispersion.mat
figure();
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);

hold on;h1=plot(freq,cr_real,'r.','markersize',6);
hold on;h2=plot(freqs,vc,'k.','markersize',3);
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

legend([h1(1),h2(1)],{'Muller','SASEM'})

load Muller_4_layer_dispersion_incomplete.mat
figure();
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.18 0.18 0.61 0.73],'color',[255 255 255]/255);

hold on;h1=plot(freq,cr_real,'r.','markersize',6);
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
box on;
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;

legend(h1(1),'Muller')