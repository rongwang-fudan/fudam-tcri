% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2024.1.1
tic
clear;
load('H:\Transfer of carbon-tax revenue\inputs\Table_PW_test1same_county.mat','Table_county')
load('H:\Transfer of carbon-tax revenue\inputs\Table_PW_area_test1same.mat','Table_PW_area')
% km2

% load('inputs\Table_all0311.mat','Table_all'); % 232187x28
load('H:\Transfer of carbon-tax revenue\inputs\Table_all0311_samecapital.mat','Table_all'); % 232187x28
Table_county(:,[7 9])-Table_all(:,[7 9]);
figure;plot(ans,'DisplayName','ans')

% load('inputs\Table_all0105.mat','Table_all'); % 232187x28
cn192id=load('H:\Transfer of carbon-tax revenue\inputs\cn192id.txt'); % 192x4
plots=0;
lifetimelcp=[25 25 25 40 40 25 25]; % lifetime 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
% 2100年全球
% 1.电厂编号；2.电厂类型（1 PV; 2 Onshore wind; 3 Offshore wind; 4 Hydroelectricity;  5 Geothermal; 6 Nuclear; 7 Biomass; 8 CCS）；
% 3.已建/未建（1为已建，2为未建）；4.lat；5.lon；6.总发电量 (TWh/y)；
% 7.power use efficency（%）；8.CP(MW)；9.Country（1-192）；10.电厂-非土地初始投资（million $）；11.电厂-土地初始投资（million $）；12.电厂-年运维成本（million $/y）；
% 13.storage类型（1代表battery,2代表mechanical）；14.storage初始投资（million $）；15.storage年运维成本（million $/y）；
% 16.trans初始投资（million $）；17.trans年运维成本（million $/y）；
% 18.Revenue（million$/yr）根据各国平均燃料价格计算
% 19.年CO2减排量（Mton CO2/y）
% 20.Primary energy (EJ/y)
% 21-28.制造光伏面板和风机的矿物质需求(Mton)（1 Copper, 2 Zinc, 3 Nickel, 4 Silicon, 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths）
load('H:\Transfer of carbon-tax revenue\inputs\FID_countryID.mat')
% 1 FID; 2 CountyID
[m,n]=find(Table_all(:,2)==3);
ccc = Table_county(m,10); % FID
Table_county(m,10) = FID_countryID(ccc+1,2);


% define the type of energy - 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
retype=[1 2 2 5 6 4 3 7];
for i=1:size(Table_all,1)
    Table_all(i,2)=retype(Table_all(i,2));
    Table_all(i,24)=Table_all(i,9); % country id
    Table_all(i,9)=cn192id(Table_all(i,9),3); % region id
end

mac=(Table_all(:,10)+Table_all(:,11)+Table_all(:,14)+Table_all(:,16))./Table_all(:,19)./lifetimelcp(1); % $/tCO2
idx=find(mac<500 & Table_all(:,19)>0 & Table_all(:,9)>0 & Table_all(:,2)<=2);
[v,id]=sort(mac(idx));

punits2=zeros(size(idx,1),31);
punits2(:,1:28)=Table_all(idx(id),1:28);
punits2_county=Table_county(idx(id),:);
punits2_area=Table_PW_area(idx(id),:);
% 1.id; 2. type of unit (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 3.status (1 built, 2 to be built); 4.lat; 5.lon; 6.power generation TWh/y
% 7.power use efficency %; 8.capacity MW; 9.region ID (2-12); 10.non-land investment m$; 11.land investment m$; 12.OM cost m$/y
% 13.storage (1 battery, 2 hydro pump); 14. storage investment m$; 15.storage OM cost m$/y
% 16.transmission investment m$; 15.transmission OM cost m$/y
% 18.Revenue of replacing fossil fuels m$/y
% 19.Abated emission Mt CO2/y
% 20.Primary energy EJ/y
% 21-28 consumption of 1 Copper, 2 Zinc, 3 Nickel, 4 Silicon, 5 Manganese, 6 Chromium, 7 Molybdenum, 8 Rare earths
% 29 total initial investment m$
% 30 fraction of non-land costs in total costs
% 31 total O&M cost m$/y
punits2(:,29)=punits2(:,10)+punits2(:,11)+punits2(:,14)+punits2(:,16); % initial investment million $
punits2(:,30)=(punits2(:,10)+punits2(:,14)+punits2(:,16))./punits2(:,29); % fraction of non-land costs in total costs
punits2(:,31)=punits2(:,12)+punits2(:,15)+punits2(:,17); % O&M cost million $

punits=zeros(size(idx,1),11);
% 1 year of building this power unit; 2.dynamic year; 3 type of unit 4 region id; 5 power generation PJ/yr; 6 total investment t$; 7 Abated emission Mt CO2/y;
% 8 fraction of non-land costs in total costs; 9 lat; 10 longitude; 11 marginal abatement cost $/tCO2
punits(:,1)=2400;
punits(:,2)=punits(:,1); % dynamic year
punits(:,3)=punits2(:,2); % type of unit 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
punits(:,4)=punits2(:,9); % region id
punits(:,5)=punits2(:,6).*punits2(:,7).*3.6; % power generation PJ/yr
punits(:,6)=punits2(:,29)*1e-6; % total investment t$
punits(:,7)=punits2(:,19)*1e-3; % abated emission GtCO2/yr
punits(:,8)=punits2(:,30); % fraction of non-land costs in total costs
punits(:,9)=punits2(:,4); % latitude
punits(:,10)=punits2(:,5); % longitude
punits(:,11)=punits(:,6)./punits(:,7).*(1000/lifetimelcp(1)); % marginal abatement cost $/tCO2
punits(:,12)=punits2_county(:,10); % county
punits(:,13)=punits2(:,24); % country
punits(:,14)=punits2(:,8); % CP(MW)
punits(:,15)=punits2_area; % km2
clear punits2;

load('H:\Transfer of carbon-tax revenue\inputs\ScenarioID_Type.mat'); % C1 to C8
enepro=zeros(12,17,1291,6); % final energy PWh/yr 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Solar_IIASA.mat');
enepro(:,:,:,1)=PrimaryEnergy_Solar_IIASA(1:12,21:5:101,:).*(0.15/3.6); clear PrimaryEnergy_Solar_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Wind_IIASA.mat');
enepro(:,:,:,2)=PrimaryEnergy_Wind_IIASA(1:12,21:5:101,:).*(0.17/3.6); clear PrimaryEnergy_Wind_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Biomass_IIASA.mat');
enepro(:,:,:,3)=PrimaryEnergy_Biomass_IIASA(1:12,21:5:101,:).*(0.3/3.6); clear PrimaryEnergy_Biomass_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Nuclear_IIASA.mat');
enepro(:,:,:,4)=PrimaryEnergy_Nuclear_IIASA(1:12,21:5:101,:).*(0.8/3.6); clear PrimaryEnergy_Nuclear_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Hydro_IIASA.mat');
enepro(:,:,:,5)=PrimaryEnergy_Hydro_IIASA(1:12,21:5:101,:).*(0.44/3.6); clear PrimaryEnergy_Hydro_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Geothermal_IIASA.mat');
enepro(:,:,:,6)=PrimaryEnergy_Geothermal_IIASA(1:12,21:5:101,:).*(0.75/3.6); clear PrimaryEnergy_Geothermal_IIASA;
emax=zeros(1,6);
renemax=zeros(12,17,6);
linecolor=jet(8);
for j=1:6
    a1=sum(sum(enepro(:,:,:,j),2),1);
    idx=find(1-isnan(a1));
    [v,id]=sort(-enepro(1,17,idx,j));
    renemax(:,:,j)=enepro(:,:,idx(id(2)),j); % maximal capacity of six types of renewable energy
    if plots==1
        subplot(3,3,j);
        for i=1:size(idx,1)
            a=sum(enepro(2:12,1:17,idx(i),j),1); plot(a,'LineStyle','-','LineWidth',0.1,'Color',linecolor(ScenarioID_Type(idx(i),2),1:3)); hold on;
        end
        a=sum(renemax(2:12,1:17,j),1); plot(a,'LineStyle','-','LineWidth',1.5,'Color',[0 0 0]); hold on;
    end
end

% determine power units for 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
Nj=50;
Nj2=2500;
totalccs=10; % 10 GtCO2/yr for CCS in each region
macccs=500; % $500-1000/tCO2 for CCS
totalhydrogen=10; % 10 PWh/yr for other renewable in each region
machydrogen=150; % $150-300/tCO2 for other renewable (e.g., sea-water hydrogen)
punitsother=zeros(Nj*17*12*5,15);
j1=1;
for i=2:12
    idxpv=find(punits(:,3)==1 & punits(:,4)==i);
    idxwd=find(punits(:,3)==2 & punits(:,4)==i);
    for t=1:17
        epv=renemax(i,t,1)*3600; % PJ/yr
        jp=1;
        while epv>0 && jp<=size(idxpv,1)
            if t==1
                punits(idxpv(jp),1)=2020-mod(jp,20);
            end
            epv=epv-punits(idxpv(jp),5); % power generation PJ/yr
            jp=jp+1;
        end
        ewind=renemax(i,t,2)*3600; % PJ/yr
        jw=1;
        while ewind>0 && jw<=size(idxwd,1)
            if t==1
                punits(idxwd(jw),1)=2020-mod(jw,20);
            end
            ewind=ewind-punits(idxwd(jw),5); % power generation PJ/yr
            jw=jw+1;
        end
        maxmac=max(punits(idxpv(jp-1),11),punits(idxwd(jw-1),11)); % % marginal abatement cost $/tCO2
        idx=find(punits(:,1)<=2020 & punits(:,3)==2 & punits(:,4)==i);
        efco2=mean(punits(idx,7),1)/mean(punits(idx,5),1); % GtCO2 per PJ energy
        % 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal and others
        for j=3:6
            j2=j1+Nj-1;
            if t==1
                ej=renemax(i,t,j)*3600; % PJ/yr
                for j0=j1:j2
                    punitsother(j0,1)=2020-mod(j0,20);
                end
            else
                ej=max(0,renemax(i,t,j)-renemax(i,t-1,j))*3600; % PJ/yr
                punitsother(j1:j2,1)=2400;
            end
            if ej==0
                continue;
            end
            punitsother(j1:j2,2)=punitsother(j1:j2,1);
            punitsother(j1:j2,3)=j; % type of unit 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
            punitsother(j1:j2,4)=i; % region id
            punitsother(j1:j2,5)=ej/Nj; % power generation PJ/yr
            punitsother(j1:j2,6)=ej/Nj*efco2*maxmac*(lifetimelcp(j)/1000); % total investment t$
            punitsother(j1:j2,7)=ej/Nj*efco2; % abated emission GtCO2/yr
            punitsother(j1:j2,11)=maxmac; % marginal abatement cost $/tCO2
            j1=j2+1;
        end
    end
    
    % CCS
    j2=j1+Nj2-1;
    punitsother(j1:j2,1)=2400;
    punitsother(j1:j2,2)=2400;
    punitsother(j1:j2,3)=7; % type of unit 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
    punitsother(j1:j2,4)=i; % region id
    punitsother(j1:j2,5)=totalccs/Nj2/efco2; % power generation PJ/yr
    punitsother(j1:j2,7)=totalccs/Nj2; % abated emission GtCO2/yr
    for j3=j1:j2
        punitsother(j3,11)=macccs*(1+(j3-j1)/Nj2); % marginal abatement cost $/tCO2
        punitsother(j3,6)=totalccs/Nj2*lifetimelcp(7)*punitsother(j3,11)/1000; % total investment t$
    end
    j1=j2+1;
    
    % 10 PWh for other renewable (e.g., sea-water hydrogen) in each region in 2100
    j2=j1+Nj2-1;
    punitsother(j1:j2,1)=2400;
    punitsother(j1:j2,2)=2400;
    punitsother(j1:j2,3)=6; % type of unit 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
    punitsother(j1:j2,4)=i; % region id
    punitsother(j1:j2,5)=totalhydrogen*3600/Nj2; % power generation PJ/yr
    punitsother(j1:j2,7)=totalhydrogen*3600/Nj2*efco2; % abated emission GtCO2/yr
    for j3=j1:j2
        punitsother(j3,11)=machydrogen*(1+(j3-j1)/Nj2); % marginal abatement cost $/tCO2
        punitsother(j3,6)=punitsother(j3,7)*lifetimelcp(6)*punitsother(j3,11)/1000; % total investment t$
    end
    j1=j2+1;
end

% idx=find(punitsother(:,3)==4);
% totalnuclear=sum(punitsother(idx,5),1)./3600;

punits2=zeros(size(punits,1)+j2,15);
punits2(1:size(punits,1),1:end)=punits;
punits2((size(punits,1)+1):(size(punits,1)+j2),1:end)=punitsother(1:j2,1:end);
[v,id]=sort(punits2(:,11));
punits=punits2(id,1:end); clear punits2;
if plots==1
    %     idx=find(punits(:,11)<500 & punits(:,3)<3);
    idx=find(punits(:,11)<150);
    a=punits(idx,1:11);
    n=size(a,1);
    q=zeros(n,2);
    q(1,1)=a(1,7); % GtCO2/yr
    q(1,2)=a(1,5)/3600; % PWh/yr
    for i=2:n
        q(i,1)=q(i-1,1)+a(i,7); % GtCO2/yr
        q(i,2)=q(i-1,2)+a(i,5)/3600; % PWh/yr
    end
    subplot(3,3,7); plot(q(:,1),a(:,11));
    subplot(3,3,8); plot(q(:,2),a(:,11));
    subplot(3,3,9); plot(punits(:,end))
end


save('H:\Transfer of carbon-tax revenue\Ans\fig3_punits.dat','punits');

