% Modal analysis for viscoelastic medal via the semi-analytical spectral
% element method

clear;
%frequency parameters
fmin=1; % minimum frequency
fmax=60; % maximum frequency
df=0.5;%  interval
freqs=fmin:df:fmax;% frequency array

% model/grid parameters
modelfile='LVL_viscoelastic_model.csv';
model_type=1;% 1 for multilayered models and 2 for gradient models 
mode_type = 2;% 1 for fundamental mode only and 2 for multi modes
output_v=1;  % whether to output the eigenwavefields (0 or 1)

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
figure();
set(gcf,'unit','centimeters','position',[10,10,14,6]);
subplot(1,2,1);set(gca,'position',[0.1 0.18 0.35 0.73])
hold on;plot3(freqs,-imag(hw),vc,'k.','markersize',3);axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');zlabel('Phase Velocity (m/s)');ylabel('Attenuation (m^{-1})')
set(gca,'fontname','times new roman','fontsize',8);box on;grid on;
view(-12,12);title('(a)')

subplot(1,2,2);set(gca,'position',[0.565 0.18 0.305 0.73])
% hold on;plot(freqs,vc,'k.','markersize',3);axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
hold on;scatter(reshape(repmat(freqs,size(vc,1),1),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);
% hold on;scatter3(reshape(repmat(freqs,size(hw,1),1),[],1),...
%     reshape(-imag(hw),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
% axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;
title('(b)')

co=colorbar;
set(get(co,'XLabel'),'string','-Im (k) (m^{-1})','Fontsize',8,'FontName','times new roman');
set(co,'Fontsize',8,'FontName','times new roman');
set(co,'Location','EastOutside');
Po=get(co,'Position');Po(1)=0.89;Po(2)=0.18;Po(3)=0.02;Po(4)=0.73;
set(co,'Position',Po);
%% model tracing
%frequency parameters
fmin=1; % minimum frequency
fmax=60; % maximum frequency
df=0.125;%  interval
freqs=fmin:df:fmax;% frequency array
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
[vc,hw,wavefields] = mode_tracing(vc,hw,wavefields);

figure;
set(gcf,'unit','centimeters','position',[10,10,7,6]);
set(gca,'position',[0.12 0.18 0.7 0.73])
for no=1:size(vc,1)
    hold on;plot3(freqs,-imag(hw(no,:)),real(vc(no,:)),'-','LineWidth',1.5)
end
axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);view(3)
box on;
set(gca,'TickDir','in','TickLength',[0.01 0.01])
xlabel('Frequency (Hz)');zlabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);

%% plot eigen wavefields (uz)
plot_wavefields=1;
if plot_wavefields==1
title_array={'(a)','(b)','(c)','(d)'};
for modeno=0:3
    figure(modeno+15);
    set(gcf,'unit','centimeters','position',[10,10,5,6]);
    set(gca,'position',[0.18 0.18 0.71 0.73])
    count=0;
    scale=3;
    for freno=1:8:length(freqs)
        if find(~isnan(vc(:,freno)),1,'last')>modeno 
            if count==0
                uz=scale*real(wavefields(freno).uz(:,modeno+1))./max(abs(wavefields(freno).uz(:,modeno+1)));
                hold on;plot(freqs(freno)+uz,wavefields(freno).solid_coor,'k-');
                uz_old=uz;coor_old=wavefields(freno).solid_coor;
                count=1;
            else            
                uz_new=scale*real(wavefields(freno).uz(:,modeno+1))./max(abs(wavefields(freno).uz(:,modeno+1)));
                coor_new=wavefields(freno).solid_coor;

                coor=linspace(0,60,100);
                coef=corrcoef(interp1(coor_old,uz_old,coor),interp1(coor_new,uz_new,coor));
                if coef(2,1)<0
                    uz_new=uz_new*-1;%disp(freno)
                end
                uz_old=uz_new;coor_old=coor_new;
                if mod(freno,2)==1
                    hold on;plot(freqs(freno)+uz_new,wavefields(freno).solid_coor,'k-');
                end
            end
        end
    end
    axis([min(freqs)-scale,max(freqs)+scale,0,60]);set(gca,'ydir','reverse')
    title(sprintf('Uz of Mode %d',modeno));
    for layerno=1:length(gmodel.thk)
        hold on;plot([min(freqs)-scale,max(freqs)+scale],ones(1,2)*sum(gmodel.thk(1:layerno)),'r--')
    end
    box on;
    
    set(gca,'TickDir','in','TickLength',[0.02 0.02])
    xlabel('Frequency (Hz)');ylabel('Depth (m)');
    set(gca,'fontname','times new roman','fontsize',8);box on;
    text(-1.7,-3.7,title_array(modeno+1),'fontname','times new roman','fontsize',9)
end

% plot eigen wavefields (ur)
for modeno=0:3
    figure(modeno+10);
    set(gcf,'unit','centimeters','position',[10,10,5,6]);
    set(gca,'position',[0.18 0.18 0.71 0.73])
    count=0;
    scale=2;
    for freno=1:8:length(freqs)
        if find(~isnan(vc(:,freno)),1,'last')>modeno 
            if count==0
                ur=scale*real(wavefields(freno).ur(:,modeno+1))./max(abs(wavefields(freno).ur(:,modeno+1)));
                hold on;plot(freqs(freno)+ur,wavefields(freno).solid_coor,'k-');
                ur_old=ur;coor_old=wavefields(freno).solid_coor;
                count=1;
            else            
                ur_new=scale*real(wavefields(freno).ur(:,modeno+1))./max(abs(wavefields(freno).ur(:,modeno+1)));
                coor_new=wavefields(freno).solid_coor;

                coor=linspace(0,60,100);
                coef=corrcoef(interp1(coor_old,ur_old,coor),interp1(coor_new,ur_new,coor));
                if coef(2,1)<0
                    ur_new=ur_new*-1;%disp(freno)
                end
                ur_old=ur_new;coor_old=coor_new;
                if mod(freno,2)==1
                    hold on;plot(freqs(freno)+ur_new,wavefields(freno).solid_coor,'k-');
                end
            end
        end
    end
    axis([min(freqs)-scale,max(freqs)+scale,0,60]);set(gca,'ydir','reverse')
    title(sprintf('Ur of Mode %d',modeno));
    for layerno=1:length(gmodel.thk)
        hold on;plot([min(freqs)-scale,max(freqs)+scale],ones(1,2)*sum(gmodel.thk(1:layerno)),'r--')
    end
    box on;
    
    set(gca,'TickDir','in','TickLength',[0.02 0.02])
    xlabel('Frequency (Hz)');ylabel('Depth (m)');
    set(gca,'fontname','times new roman','fontsize',8);box on;
    text(-1.7,-3.7,title_array(modeno+1),'fontname','times new roman','fontsize',9)
end
end

%% elastic model
clear;
%frequency parameters
fmin=1; % minimum frequency
fmax=60; % maximum frequency
df=0.5;%  interval
freqs=fmin:df:fmax;% frequency array

% model/grid parameters
modelfile='LVL_elastic_model.csv';
model_type=1;% 1 for multilayered models and 2 for gradient models 
mode_type = 2;% 1 for fundamental mode only and 2 for multi modes
output_v=1;  % whether to output the eigenwavefields (0 or 1)

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

figure();
set(gcf,'unit','centimeters','position',[10,10,14,6]);
subplot(1,2,1);set(gca,'position',[0.1 0.18 0.305 0.73])
% hold on;plot(freqs,vc,'k.','markersize',3);axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)])
hold on;scatter(reshape(repmat(freqs,size(vc,1),1),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);caxis([0,0.0331])
% hold on;scatter3(reshape(repmat(freqs,size(hw,1),1),[],1),...
%     reshape(-imag(hw),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
% axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;
title('(a)')
co=colorbar;
set(get(co,'XLabel'),'string','-Im (k) (m^{-1})','Fontsize',8,'FontName','times new roman');
set(co,'Fontsize',8,'FontName','times new roman');
set(co,'Location','EastOutside');
Po=get(co,'Position');Po(1)=0.89;Po(2)=0.18;Po(3)=0.02;Po(4)=0.73;
set(co,'Position',Po);

subplot(1,2,2);set(gca,'position',[0.565 0.18 0.305 0.73])
fmin=35; % minimum frequency
fmax=60; % maximum frequency
df=0.05;%  interval
freqs=fmin:df:fmax;% frequency array
tic;
[vc,hw,wavefields]=sasem_psv(gmodel,freqs,model_type,mode_type,output_v);
toc;
hold on;scatter(reshape(repmat(freqs,size(vc,1),1),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
axis([-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);caxis([0,0.0331])
% hold on;scatter3(reshape(repmat(freqs,size(hw,1),1),[],1),...
%     reshape(-imag(hw),[],1),reshape(vc,[],1),3,reshape(-imag(hw),[],1),'filled');
% axis([-inf,inf,-inf,inf,0.9*min(gmodel.vs),1.05*max(gmodel.vs)]);colormap(jet);view(3)
set(gca,'TickDir','in','TickLength',[0.02 0.02])
xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
set(gca,'fontname','times new roman','fontsize',8);box on;
title('(b)')
axis([38,45,360,370])


