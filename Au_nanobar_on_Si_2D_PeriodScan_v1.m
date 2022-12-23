% Au nanobars on on Si substrate
% H.Altug, et al, Optical Express 19(2011), 11207
clear
clc
nk0 = readtable('Au-nk-Weaver-2015.csv','VariableNamingRule','preserve');
%%
figure();
plot(nk0.("Wavelength, µm"), nk0.k,'bo',nk0.("Wavelength, µm"), nk0.n,'ro')
%%
%wv=linspace(5.4, 7.6, 10)';
%wv=(5.4:0.16:6.9)';
wv=10000./[1250,1300,1350,1400,1450,1500,1550,1580,1600,1610,1620,1630,1640,1645,1650,1655,1660,1665,1670,1675,1680,1690,1700,1720,1750,1800,1850]';
nk_N=interp1(nk0.("Wavelength, µm"),nk0.n,wv);
nk_K=interp1(nk0.("Wavelength, µm"),nk0.k, wv);
%%
clear nk0
%nk=nk(31:71,:);
figure();
plot(wv, nk_K,'bo', wv,nk_N, 'ro')
%%
figure();
plot(10000./wv, nk_K,'bo', 10000./wv,nk_N, 'ro')
%%
%period=[1.75,2];% same unit as wavelength
n_incident_medium=1.0;% refractive index of the top layer
n_transmitted_medium=3.41;% refractive index of the bottom layer
angle_theta0=0; % in degrees
angle_delta=90;
k_parallel=n_incident_medium*sin(angle_theta0*pi/180);
parm=res0;  %  polarization, i.e. E field component is perpendicular to grating direction
parm.not_io=1;% prevents writing of intermediate files onto the disk
parm.res1.champ=1;% the electromagnetic field is calculated accurately
nn=[20,10];% Fourier harmonics run from [-nn,nn]
%%
% textures for all layers including the top and bottom layers
textures=cell(length(wv),3);
textures(:,1)= {n_incident_medium}; % uniform texture
textures(:,3)= {n_transmitted_medium}; % uniform texture
%txt=cell(length(wv),2);
%%
for m=1:length(wv)
    n_gold=nk_N(m)+1i*nk_K(m);
    textures{m,2}={n_incident_medium,[0,0,0.87,0.23,n_gold,1]};
end
%%
% setup profile
profile ={[0,0.07,0],[1,2,3]}; % thickness and texture labels
eff_te_ref=zeros(length(wv),1);
%eff_tm_ref=zeros(length(wv),1);
%eff_te_trans=zeros(length(wv),1);
%eff_tm_trans=zeros(length(wv),1);
%%
%fileID = fopen('exptable.txt','w');
pool = parpool('local');
%%
Py=[1.75,1.81,1.88,1.94,2];
dat=zeros(length(wv),length(Py));
dat(:,1)=wv;
%%
tic
for l=1:length(Py)
    parfor m=1:length(wv)
        [aa,neff]=res1(wv(m),[1.75,Py(l)],textures(m,:),nn,k_parallel,angle_delta,parm);
        result=res2(aa,profile);
        eff_te_ref(m)=result.TEinc_top_reflected.efficiency{0,0};
        %eff_tm_ref(m)=result.TMinc_top_reflected.efficiency{0,0};
        %eff_te_trans(m)=result.TEinc_top_transmitted.efficiency{0,0};
        %eff_tm_trans(m)=result.TMinc_top_transmitted.efficiency{0,0};
        
    end
    fprintf('the Y period is %1.2f\n',Py(l))
    
    dat(:,l)=eff_te_ref;

end
toc
%%
dat=[wv,dat];
%%
%delete(pool);
set(0,'defaultlinelinewidth',1,'defaultTextFontSize',12,'DefaultAxesFontSize',...
    12,'DefaultFigureColor','white','DefaultLineLineSmoothing','on',...
'defaultlegendFontSize',10);
figure('InvertHardcopy','off','PaperUnits', 'centimeters','PaperType',...
     'usletter','PaperPosition', [0 0 12 12], 'color','white');
newcolors={'#0072BD','#FF0000','#EDB120','#77AC30','#000000'};
colororder(newcolors)
ang=0;
for k=2:6
    
    plot(10000./wv,dat(:,k),'-o',LineWidth=1,MarkerSize=4,DisplayName=num2str(Py(k-1)));
    hold on
    ang=ang+5;
end
legend(location='northwest');
%%
ylim([0.25,0.9]);
xlim([1250,1850]);
xlabel('wavenumber(cm^{-1})');ylabel('reflectance');
xticks([1250,1450,1650,1850]);
yticks(0.3:0.1:0.9);
grid on
saveas(gcf,'Au_nanobarL870W230H70_Px1750Py1750-2000_te_reflectance-efficiency-theta0nn-20-10.tiff')
%%
T=array2table(dat,'VariableNames',{'wavelength','Py1.75','Py1.81','Py1.88','Py1.94','Py2.0'});
writetable(T,'Au_nanobarL870W230H70_Px1750Py1750-2000_te_reflectance-efficiency-theta0nn-20-10.csv')


