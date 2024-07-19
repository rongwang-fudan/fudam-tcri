tic
clear;
clear global;
Initialset_econ;
Initialset_clim;
fracinv0=[0 300 52 45 39 19 325 175 10 14 10 10]; % Figure 20 in Global Landscape of Climate Finance 2023, Climate Policy Initiative, https://www.climatepolicyinitiative.org/publication/global-landscape-of-climate-finance-2023/
fracinv0(1,1)=sum(fracinv0(1,2:12),2); fracinv0=fracinv0./fracinv0(1,1);
inv15=[0.340 0.263 0.351 0.322 0.329 0.348 0.430 0.499; 75.28 76.52 81.48 86.54 87.78 85.26 97.53 101.33; 0.00452 0.00344 0.00431 0.00372 0.00375 0.00408 0.00441 0.00492]; % 1 renewable energy investment t$; 2 gdp; 3 ratio from 2015 to 2022
t15=[2015:2022]; [sR,lrv0,bb0] = regression(t15,inv15(3,:));
lrv=lrv0*2; % rate of growth in renewable energy investment

load('H:\Transfer of carbon-tax revenue\inputs\global_carbon_budget_1959_2019.txt'); % % 36.5 GtCO2/yr Friedlingstein, et al., Global Carbon Budget 2020, https://doi.org/10.5194/essd-12-3269-2020.
global_carbon_budget_1959_2019(:,2)*44/12; % 36.5 GtCO2/yr

load('H:\Transfer of carbon-tax revenue\inputs\punits.dat','-mat');
% 1 year of building this power unit; 2.dynamic year; 3 type of unit( % type of unit 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS)
% 4 region id; 5 effective power generation PJ/yr; 6 total initial investment t$; 7 Abated emission Gt CO2/y;
% 8 fraction of non-land costs in total costs; 9 lat; 10 longitude; 11 marginal abatement cost $/tCO2
punits2 = zeros(size(punits,1),23);
punits2(:,2) = punits(:,3); % type,1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
% define the type of energy - 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
punits2(:,23) = punits(:,4); % region id
EF = [0.15 0.17 0.3 0.8 0.44 0.75 0.46];
for i = 1:7
    [m,n]=find(punits2(:,2)==i);
    punits2(m,20) = punits(m,5)./1000/EF(i); % Primary energy, EJ/y
end
punits2(:,5) = punits(:,5)./1000; % effective power generation, EJ/y
punits2(:,6) = punits(:,5)/3.6; % effective power generation, TWh/y
punits2(:,19) = punits(:,7)*1000; % Abated emission Mt CO2/y;
punits2(:,21) = punits(:,6)*10^6; % total initial investment m$
pplifetime=[25,25,25,40,40,25,25];
for i = 1:7
    [m,n]=find(punits2(:,2)==i);
    punits2(m,22) = punits2(m,21)./pplifetime(i); % O&M cost, m$/y
    %   punits2(m,1) = pplifetime(i); % y
end
% mac=punits2(:,21)./punits2(:,19)./punits2(:,1);
clear punits
punits = punits2;

load('inputs\ScenarioID_Type.mat');
% Type=1→C1: Limiting warming to 1.5℃ (>50%) with no or limited overshoot
% Type=2→C2: return warming to 1.5℃ (>50%) after a high overshoot
% Type=3→C3: limit warming to 2℃ (>67%)
% Type=4→C4: limit warming to 2℃ (>50%)
% Type=5→C5: limit warming to 2.5℃ (>50%)
% Type=6→C6: limit warming to 3℃ (>50%)
% Type=7→C7: limit warming to 4℃ (>50%)
% Type=8→C8: exceed warming to 4℃ (>=50%)
% 整理后的IIASA数据均为 42*101*1291
% 其中42为区域ID，101为年份（2000-2100），1291为模型情景ID（ScenarioID_Type.mat中第一列为模型情景ID，
% 第二列为对应的升温范围），当模型情景没有数据时赋值为NaN，区域ID和模型情景ID对应的名字在“ScenarioID_Type & Region ID.xlsx”中，
enepro=zeros(12,17 ,1291,12);
% 1 coal, 2 oil, 3 gas, 4 biomass, 5 solar, 6 wind, 7 nuclear, 8 hydropower, 9 geothermal
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Coal_IIASA.mat');
% Primary Energy|Coal, EJ/yr
enepro(:,:,:,1)=PrimaryEnergy_Coal_IIASA(1:12,21:5:101,:).*(0.46); clear PrimaryEnergy_Coal_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Oil_IIASA.mat')
% Primary Energy|Oil, EJ/yr
enepro(:,:,:,2)=PrimaryEnergy_Oil_IIASA(1:12,21:5:101,:).*(0.46); clear PrimaryEnergy_Oil_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Gas_IIASA.mat')
% Primary Energy|Gas, EJ/yr
enepro(:,:,:,3)=PrimaryEnergy_Gas_IIASA(1:12,21:5:101,:).*(0.46); clear PrimaryEnergy_Gas_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Biomass_IIASA.mat')
% Primary Energy|Biomass, EJ/yr
enepro(:,:,:,4)=PrimaryEnergy_Biomass_IIASA(1:12,21:5:101,:).*0.3; clear PrimaryEnergy_Biomass_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Solar_IIASA.mat')
% Primary Energy|Solar, EJ/yr
enepro(:,:,:,5)=PrimaryEnergy_Solar_IIASA(1:12,21:5:101,:).*0.15; clear PrimaryEnergy_Solar_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Wind_IIASA.mat')
% Primary Energy|Wind, EJ/yr
enepro(:,:,:,6)=PrimaryEnergy_Wind_IIASA(1:12,21:5:101,:).*0.17; clear PrimaryEnergy_Wind_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Nuclear_IIASA.mat')
% Primary Energy|Nuclear, EJ/yr
enepro(:,:,:,7)=PrimaryEnergy_Nuclear_IIASA(1:12,21:5:101,:).*0.8; clear PrimaryEnergy_Nuclear_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Hydro_IIASA.mat')
% Primary Energy|Hydro, EJ/yr
enepro(:,:,:,8)=PrimaryEnergy_Hydro_IIASA(1:12,21:5:101,:).*0.44; clear PrimaryEnergy_Hydro_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_Geothermal_IIASA.mat')
% Primary Energy|Geothermal, EJ/yr
enepro(:,:,:,9)=PrimaryEnergy_Geothermal_IIASA(1:12,21:5:101,:).*0.75; clear PrimaryEnergy_Geothermal_IIASA;
load('H:\Transfer of carbon-tax revenue\inputs\PrimaryEnergy_total_IIASA.mat')
% Primary Energy|total, EJ/yr
enepro(:,:,:,10)=PrimaryEnergy_total_IIASA(1:12,21:5:101,:); clear PrimaryEnergy_total_IIASA;
load('inputs\CO2_Energy_IIASA.mat')
% Emissions|CO2, Mt CO2/yr
enepro(:,:,:,11)=CO2_Energy_IIASA(1:12,21:5:101,:); clear CO2_Energy_IIASA;
load('inputs\GDP_PPP_IIASA.mat')
% GDP|PPP, billion US$2010/yr
enepro(:,:,:,12)=GDP_PPP_IIASA(1:12,21:5:101,:);
GDP_iiasa = GDP_PPP_IIASA(1:12,21:1:101,:)/1000; % trillion $
clear GDP_PPP_IIASA;


load('H:\Transfer of carbon-tax revenue\Ans\Fig2data_FUMCCM_Prienergy_reg_CUR.mat','Prienergy_reg')
% energy PWh
% FUMCCM中CUR情景所用到的一次能源量，作为基础进行校正
Prienergy_reg_FUMCCM_CUR = Prienergy_reg(1:5:81,:)';
enepro2a = enepro;
for i = 1:size(enepro,3)
    for j = 1:9
        enepro2a(:,:,i,j) = enepro(:,:,i,10)*0.3263./sum(enepro(:,:,i,1:9),4).*enepro(:,:,i,j);
    end;
    enepro(:,:,i,:) = enepro2a(:,:,i,:);
    
    rrr(:,:,i) = Prienergy_reg_FUMCCM_CUR./(enepro(:,:,i,10)/3.6*0.3263);
    for j = 1:10
        enepro2(:,:,i,j) = enepro(:,:,i,j).*rrr(:,:,i);
    end
    %     enepro2(:,:,i,11) = enepro(:,:,i,11)+enepro(:,:,i,11).*(rrr(:,:,i)-1).*sum(enepro(:,:,i,1:3),4)./sum(enepro(:,:,i,1:9),4);
    enepro2(:,:,i,11) = enepro(:,:,i,11).*rrr(:,:,i);
    enepro2(:,:,i,12) = enepro(:,:,i,12);
end
clear enepro
enepro = enepro2;
clear enepro2

% find a regional profile for (globe + 11 regions) x (10 energies + co2 emission + gdp)
energyprofile=zeros(12,12);
sce_yesss = zeros(1291,1);
nn = 1;
for sce=1:1291
    %     if ScenarioID_Type(sce,2)>2 % warming >1.5℃
    for i=1:10
        energyprofile(1:12,i)=enepro(1:12,17,sce,i); % 2100
    end
    energyprofile(1:12,11)=enepro(1:12,1,sce,11); % 2020, Emissions|CO2, Mt CO2/yr
    energyprofile(1:12,12)=enepro(1:12,1,sce,12); % 2020, GDP|PPP, billion US$2010/yr
    idx=find(isnan(energyprofile));
    if size(idx,1)==0
        sce_comp(nn,1) = sce;
        sce_comp(nn,2) = ScenarioID_Type(sce,2);
        nn= nn+1;
        sce_yesss(sce,1) = 1;
        %             break;
    end
    %     end
end % choose one scenario
sce_yesss00 = sce_yesss;

%
for i = [1 6]
    [m,n]=find(ScenarioID_Type(:,2)==i);
    a = reshape(enepro(1,(2050-2015)/5,m,11),[size(m,1) 1]);
    z=find(~isnan(a));
    amedian=median(a(z));
    amin=min(a(z));
    amax=max(a(z));
    
    [m,n]=find(ScenarioID_Type(:,2)==i & sce_yesss==1);
    sce_id = ScenarioID_Type(m,:);
    a = reshape(enepro(1,(2050-2015)/5,m,11),[size(m,1) 1]);
    
    C=abs(a-amedian);
    [m2,n2]=find(C==min(C));
    scee_s8(i,1) = sce_id(m2(1),1); % 2050年所有iiasa情景能源总需求 PWh/y的Median
    
    C=abs(a-amin);
    [m2,n2]=find(C==min(C));
    scee_s8min(i,1) = sce_id(m2(1),1); % 2050年所有iiasa情景能源总需求 PWh/y的Median
    
    C=abs(a-amax);
    [m2,n2]=find(C==min(C));
    scee_s8max(i,1) = sce_id(m2(1),1); % 2050年所有iiasa情景能源总需求 PWh/y的Median
end

% find a regional profile for (globe + 11 regions) x (10 energies + co2 emission + gdp)
energyprofile=zeros(12,12);
for sce=1:1291
    if ScenarioID_Type(sce,2)>2
        for i=1:10
            energyprofile(1:12,i)=enepro(1:12,17,sce,i);
        end
        energyprofile(1:12,11)=enepro(1:12,1,sce,11);
        energyprofile(1:12,12)=enepro(1:12,1,sce,12);
        idx=find(isnan(energyprofile));
        if size(idx,1)==0
            break;
        end
    end
end


n3=0;
%     learningrate=0.23; % learning rate Victoria M , Haegel N , Peters I M ,et al.Solar photovoltaics is ready to power a sustainable future. 2021
learningrate=[0.2,0.2,0.1,0.01,0.1,0.1,0.01];
financeexchange=1; % 1 allowing carbon tax revenue to flow from one region to other regions; 2 no fund flow
option_omcost=2; % 1 using carbon tax revenue to cover OM cost; 2 assume that the income of electricity cover OM cost
%     pplifetime=25; % lifetime of low-carbon power units (100: assume that the income of electricity cover OM cost and reconstruction / 25: assume that the income of electricity cannot cover OM cost and reconstruction)
pplifetime=[25,25,25,40,40,25,25];
fflifetime=40; % 40 years for power plants and industrial boilers Davis S J,et al.Committed Emissions from Existing Energy Infrastructure Jeopardize 1.5 °C Climate Target
%     carbonpricetype=2; % 1 global carbon price; 2 carbon tax as a percentage of GDP in 2020
carbonpricetype=3; % 1 global carbon price; 2 carbon tax as a percentage of GDP in 2020
frac_unsub=0; % fraction of fossil fuel that cannot be substituded

[mbmb,n]=find(sce_yesss==1);
save('H:\Transfer of carbon-tax revenue\Ans\iiasa_AGR_data.mat')

%%
tic
clear;

load('H:\Transfer of carbon-tax revenue\Ans\iiasa_AGR_data.mat')
enepro22= enepro*0;
invREpath_sce_reg = zeros(17,11,1291);
invREiiasa_sce_reg = zeros(17,11,1291);
invREpath_sce = zeros(1291,17);
invREiiasa_sce = zeros(1291,17);
invREpath_sce_ori = zeros(1291,17);
invREiiasa_sce_ori = zeros(1291,17);
bb=zeros(8,16);
for exp=1%:1
    co2path=zeros(12,17,17,1291); % 12 regions; 17 years from 2020, 2025 to 2100; 16 carbon prices
    energypath=zeros(12,17,17,10,1291); % fossil fuel energy EJ/yr
    costREpath=zeros(12,17,17,10,1291); % fossil fuel energy EJ/yr
    costRE_iiasa=zeros(17,10,1291); % million $
    costRE_iiasa_reg=zeros(17,12,10,1291); % million $
    % 1 coal, 2 oil, 3 gas, 4 solar, 5 wind, 6 biomass, 7 nuclear, 8 hydropower, 9 geothermal, 10 CCS
    n3=0;
    %     learningrate=0.23; % learning rate Victoria M , Haegel N , Peters I M ,et al.Solar photovoltaics is ready to power a sustainable future. 2021
    learningrate=[0.2,0.2,0.1,0.01,0.1,0.1,0.01];
    financeexchange=1; % 1 allowing carbon tax revenue to flow from one region to other regions; 2 no fund flow
    option_omcost=2; % 1 using carbon tax revenue to cover OM cost; 2 assume that the income of electricity cover OM cost
    %     pplifetime=25; % lifetime of low-carbon power units (100: assume that the income of electricity cover OM cost and reconstruction / 25: assume that the income of electricity cannot cover OM cost and reconstruction)
    pplifetime=[25,25,25,40,40,25,25];
    fflifetime=40; % 40 years for power plants and industrial boilers Davis S J,et al.Committed Emissions from Existing Energy Infrastructure Jeopardize 1.5 °C Climate Target
    %     carbonpricetype=2; % 1 global carbon price; 2 carbon tax as a percentage of GDP in 2020
    carbonpricetype=3; % 1 global carbon price; 2 carbon tax as a percentage of GDP in 2020
    frac_unsub=0; % fraction of fossil fuel that cannot be substituded
    
    
    %     [m1,n]=find(ScenarioID_Type(:,2)==1 & sce_yesss==1);
    %     [m2,n]=find(ScenarioID_Type(:,2)==6 & sce_yesss==1);
    %     mbmb = [m1;m2];
    [mbmb,n]=find(sce_yesss==1);
    %     for sss = 1%:size(mbmb,1)
    for sss = 1:size(mbmb,1)
        load('H:\Transfer of carbon-tax revenue\Ans\iiasa_AGR_data.mat')
    if exp==2
        financeexchange=2; % no fund flow
    elseif exp==3
        fflifetime=1; % rapid retirement of fossil fuel
    elseif exp==4
        fflifetime=1; % rapid retirement of fossil fuel
        financeexchange=2; % no fund flow
    elseif exp==5
        learningrate=0; % no learning
    elseif exp==6
        pplifetime=100; % long lifetime of low-carbon power plants
    elseif exp==7
        frac_unsub=0.1; % all fossil fuel can be substituded
    elseif exp==8
        carbonpricetype=1; % carbon prices from 0, $5, $10, ... to $75 per tCO2
    end
        rrr = ones(size(punits,1),1);
        sce=mbmb(sss); % 768 for type 4, 963 for type 5 and 1209 for type 7
%         display(sce);
        lev=ScenarioID_Type(sce,2);
        idx1=find(isnan(enepro(1,1:17,sce,1:12))); % global total
        idx2=find(isnan(sum(enepro(:,:,sce,1:3),4))); % regional coal, oil and gas
        if size(idx1,1)>0 || size(idx2,1)>0
            continue;
        end
        n3=n3+1;
        
        cp=17;
        idyear=ones(size(punits,1),2).*2200; % year of building each power unit by 1 IIASA, 2 carbon tax revenue
        rrr_finer = zeros(size(punits,1),17);
        % find id of power units that have been considered in the IIASA model
        idyear(:,2)=2200; % initializing the year to build each power unit
        pff=zeros(17,12,3); % fossil fuel energy EJ/yr that can be replaced by low-carbon energy after 2025
        co2path(1:12,1:17,cp,sce)=enepro(1:12,1:17,sce,11)./1e3; % Gt CO2/yr initializing for the starting year 2020
        for j=1:9
            energypath(:,:,cp,j,sce)=enepro(:,:,sce,j); % initializing total energy EJ/yr
        end
        
        punits1 = punits;
        inv_supply_all = zeros(12,17);
        
        for t= 1:17
            tyear = 2011+t*5; % t=1为起始值，t=2是2021-2025
            for i=1:12
                clear inv_supply
                if t>2
                    inv_a(1)=(inv15(1,2022-2014)*fracinv0(i)/GDP_iiasa(i,2022-2019,sce)+lrv*(tyear-1-2022))*GDP_iiasa(i,tyear-1-2019,sce);
                    inv_supply = sum(inv_a)*10^6*5; % million $
                    inv_supply_all(i,t) = inv_supply/2;
                end
                
                
                for j=1:6
                    if t==1
                        newpowerdemand=max(0,enepro(i,t,sce,j+3)); % EJ/yr
                    else
                        newpowerdemand=max(0,enepro(i,t,sce,j+3)-max(enepro(i,1:t-1,sce,j+3))); % EJ/yr
                    end
                    enepro22(i,t,sce,j+3) = min(enepro(i,t,sce,j+3), enepro22(i,t,sce,j+3));
                    
                    %                     idx=find(punits(:,2)==j & punits(:,23)==i & idyear(:,1)>2100 & rrr~=0); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal
                    idx=find(punits(:,2)==j & punits(:,23)==i & rrr~=0); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal
                    if size(idx,1)==0
                        continue;
                    end
                    %                     idx=find(punits(:,2)==j & punits(:,23)==i & rrr==1); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal
                    te=0; % total energy EJ/yr
                    for k=1:size(idx,1)
                        if te+rrr(idx(k),1)*punits(idx(k),5)<=newpowerdemand
                            te=te+rrr(idx(k),1)*punits(idx(k),5); % EJ/yr
                            idyear(idx(k),1)=t*5+2015; % year of building this power unit
                            rrr_finer(idx(k),t) = rrr(idx(k),1);
                            enepro22(i,t:end,sce,j+3) = enepro22(i,t:end,sce,j+3) + rrr(idx(k),1).*punits(idx(k),5);
                            rrr(idx(k),1) = 0;
                        else if te+rrr(idx(k),1)*punits(idx(k),5)>newpowerdemand && te<newpowerdemand
                                rrr(idx(k),1) = rrr(idx(k),1)-(newpowerdemand-te)/punits(idx(k),5);
                                idyear(idx(k),1)=t*5+2015; % year of building this power unit
                                rrr_finer(idx(k),t) = (newpowerdemand-te)/punits(idx(k),5);
                                % enepro22(i,t,sce,j+3) = enepro22(i,t,sce,j+3) + punits(idx(k),5)*(1-rrr(idx(k)));
                                enepro22(i,t:end,sce,j+3) = enepro22(i,t:end,sce,j+3) + newpowerdemand-te;
                                te=newpowerdemand; % EJ/yr
                                %                             break;
                            end
                        end
                    end
                end
                
                for j = 1:6
                    %                     aa2 = enepro(i,t,sce,j+3)-enepro22(i,t,sce,j+3);
                    aa2 = energypath(i,t,cp,j+3,sce)-enepro22(i,t,sce,j+3);
                    aa2(aa2<0)=0;
                    an = sum(aa2);
                    
                    if an>0.001
                        %                         idx=find(punits(:,23)==i & rrr==1); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal
                        idx=find(punits(:,23)==i & rrr~=0); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal
                        if size(idx,1)~=0
                            te=0; % total energy EJ/yr
                            for k=1:size(idx,1)
                                if te+rrr(idx(k),1)*punits(idx(k),5)<=an
                                    te=te+rrr(idx(k),1)*punits(idx(k),5); % EJ/yr
                                    idyear(idx(k),1)=t*5+2015; % year of building this power unit
                                    rrr_finer(idx(k),t) = rrr_finer(idx(k),t)+rrr(idx(k),1);
                                    %                                     rrr_finer(idx(k),t) = rrr_finer(idx(k),t)+rrr(idx(k),1);
                                    
                                    enepro22(i,t:end,sce,punits(idx(k),2)+3) = enepro22(i,t:end,sce,punits(idx(k),2)+3) + rrr(idx(k),1)*punits(idx(k),5);
                                    energypath(i,t:end,cp,punits(idx(k),2)+3,sce)=energypath(i,t:end,cp,punits(idx(k),2)+3,sce)+ rrr(idx(k),1)*punits(idx(k),5);
                                    energypath(i,t:end,cp,j+3,sce)=energypath(i,t:end,cp,j+3,sce) - rrr(idx(k),1)*punits(idx(k),5);
                                    rrr(idx(k),1) = 0;
                                else if te+rrr(idx(k),1)*punits(idx(k),5)>an && te<an
                                        rrr(idx(k),1) = rrr(idx(k),1)-(an-te)/punits(idx(k),5);
                                        idyear(idx(k),1)=t*5+2015; % year of building this power unit
                                        % enepro22(i,t,sce,j+3) = enepro22(i,t,sce,j+3) + punits(idx(k),5)*(1-rrr(idx(k)));
                                        rrr_finer(idx(k),t) = rrr_finer(idx(k),t)+(an-te)/punits(idx(k),5);
                                        enepro22(i,t:end,sce,punits(idx(k),2)+3) = enepro22(i,t:end,sce,punits(idx(k),2)+3) + an-te;
                                        energypath(i,t:end,cp,punits(idx(k),2)+3,sce)=energypath(i,t:end,cp,punits(idx(k),2)+3,sce)+ an-te;
                                        energypath(i,t:end,cp,j+3,sce)=energypath(i,t:end,cp,j+3,sce) - (an-te);
                                        te=an; % EJ/yr
                                        %                             break;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            costlearning=ones(7,1);
            if t>=3
                for j=1:6
                    costlearning(j)=max((sum(energypath(2:12,t-1,cp,j+3,sce),1)./sum(energypath(2:12,2,cp,j+3,sce),1)),1)^(log2(1-learningrate(j)));
                end
            end
            
            % 比较投资供应是否满足需求
            for i = 1:12
                if t ==1
                    for eerre = 1:6
                        %                     [m,n]=find(punits(:,2)==eerre & idyear(:,1)==t*5+2015 & punits(:,23)==i);
                        [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                        %                     costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*(1-rrr(m))); % total initial investment m$
                        %                     costRE_iiasa_reg(t,eerre+3,sce,i)=costRE_iiasa_reg(t,eerre+3,sce,i)+sum(punits(m,21).*(1-rrr(m))); % total initial investment m$
                        costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).* rrr_finer(m,t)); % total initial investment m$
                        costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                    end
                end
                
                
                if t==2
                    for eerre = 1:6
                        %                     [m,n]=find(punits(:,2)==eerre & idyear(:,1)==t*5+2015 & punits(:,23)==i);
                        %                     costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*(1-rrr(m))); % total initial investment m$
                        %                     costRE_iiasa_reg(t,eerre+3,sce,i)=costRE_iiasa_reg(t,eerre+3,sce,i)+sum(punits(m,21).*(1-rrr(m))); % total initial investment m$
                        [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                        costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                        costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                    end
                else if t>2
                        for eerre = 1:6
                            %                         [m,n]=find(punits(:,2)==eerre & idyear(:,1)==t*5+2015 & punits(:,23)==i);
                            %                         costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*(1-rrr(m)))*costlearning(eerre); % total initial investment m$
                            %                         costRE_iiasa_reg(t,eerre+3,sce,i)=costRE_iiasa_reg(t,eerre+3,sce,i)+sum(punits(m,21).*(1-rrr(m)))*costlearning(eerre); % total initial investment m$
                            [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                            costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t))*costlearning(eerre); % total initial investment m$
                            costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t))*costlearning(eerre); % total initial investment m$
                        end
                    end
                end
                
                
                if t>=3
                    if sum(costRE_iiasa_reg(t,i,4:9,sce))<=inv_supply_all(i,t) % 投资供应超过安装需求
                        ctr(t,i) = inv_supply_all(i,t) - sum(costRE_iiasa_reg(t,i,4:9,sce));
                    else
                        %                     ctr(t,i) = 0;
                        aaa = zeros(size(punits,1),1);
                        for eerre = 1:6
                            [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                            aaa(m,1) = (punits(m,21).*rrr_finer(m,t)).*costlearning(eerre);
                        end
                        [m,n] = find(cumsum(aaa)<=inv_supply_all(i,t) & punits(:,23)==i & rrr_finer(:,t)~=0);
                        ctr(t,i) = inv_supply_all(i,t)-sum(aaa(m));
                        [m,n] = find(cumsum(aaa)>inv_supply_all(i,t) & punits(:,23)==i & rrr_finer(:,t)~=0);
                        idyear(m,1) = 2200;
                        te=te-sum(rrr_finer(m,t).*punits(m,5)); % EJ/yr
                        rrr(m,1) = rrr(m,1)+rrr_finer(m,t);
                        for eerre = 1:6
                            %                             aaa_reg = aaa*0;
                            %                             [m,n] = find(punits(:,2)==eerre);
                            %                             aaa_reg(m) = aaa(m);
                            [m,n] = find(punits(:,2)==eerre & cumsum(aaa)>inv_supply_all(i,t) & punits(:,23)==i & rrr_finer(:,t)~=0);
                            if ~isempty(m)
                                enepro22(i,t:end,sce,eerre+3) = enepro22(i,t:end,sce,eerre+3) - sum(rrr_finer(m,t).*(punits(m,5)));
                                energypath(i,t:end,cp,1,sce)=energypath(i,t:end,cp,1,sce)+ sum(rrr_finer(m,t).*(punits(m,6)*3.6e-3)); % primary energy oil EJ/yr
                                energypath(i,t:end,cp,eerre+3,sce)=energypath(i,t:end,cp,eerre+3,sce)- sum(rrr_finer(m,t).*punits(m,5)); % primary energy EJ/yr
                                co2path(i,t:end,cp,sce)=co2path(i,t:end,cp,sce)+ sum(rrr_finer(m,t).*(punits(m,19)))/1e3; % Gt CO2/yr
                                co2path(1,t:end,cp,sce)=co2path(1,t:end,cp,sce)+ sum(rrr_finer(m,t).*(punits(m,19)))/1e3; % Gt CO2/yr
                            end
                            rrr_finer(m,t) = 0;
                        end
                    end
                end
            end
            energypath(energypath<0)=0;
        end
        ctr(:,1) = sum(ctr(:,2:12),2);
        
        
        rrr_ori = rrr;
        % 16 carbon prices
        for cp=17%1:16
            idyear(:,2)=2200; % initializing the year to build each power unit
            %             ctr=zeros(17,12); % remaining carbon tax revenue million $
            %             pff=zeros(17,12,3); % fossil fuel energy EJ/yr that can be replaced by low-carbon energy after 2025
            %             co2path(1:12,1,cp,sce)=enepro(1:12,1,sce,11)./1e3; % Gt CO2/yr initializing for the starting year 2020
            %             for j=1:9
            %                 energypath(:,:,cp,j,sce)=enepro(:,:,sce,j); % initializing total energy EJ/yr
            %             end
            costRE_iiasa1 = costRE_iiasa;
            costRE_iiasa_reg1 = costRE_iiasa_reg;
            t = 1;
            costRE_iiasa=zeros(17,10,1291); % million $
            costRE_iiasa_reg=zeros(17,12,10,1291); % million $
            for i = 1:12
                for eerre = 1:6
                    [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                    costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                    costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                end
            end

            
            for t=2:17
                rrr1 = rrr;
                % effects of technological change on the costs of low-carbon energy based on the Wright's law
                costlearning=ones(7,1);
                if t>=3
                    for j=1:6
                        costlearning(j)=max(sum(energypath(2:12,t-1,cp,j+3,sce),1)/sum(energypath(2:12,2,cp,j+3,sce),1),1)^(log2(1-learningrate(j)));
                    end
                end
                
                for i = 2:12
                    if t==2
                        for eerre = 1:6
                            [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                            costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                            costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t)); % total initial investment m$
                        end
                    else
                        for eerre = 1:6
                            [m,n]=find(punits(:,2)==eerre & rrr_finer(:,t)~=0 & punits(:,23)==i);
                            costRE_iiasa(t,eerre+3,sce)=costRE_iiasa(t,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t))*costlearning(eerre); % total initial investment m$
                            costRE_iiasa_reg(t,i,eerre+3,sce)=costRE_iiasa_reg(t,i,eerre+3,sce)+sum(punits(m,21).*rrr_finer(m,t))*costlearning(eerre); % total initial investment m$
                        end
                    end
                    if t>2
                    ctr(t,i) = inv_supply_all(i,t) - sum(costRE_iiasa_reg(t,i,4:9,sce));
                    end
                end
                ctr(t,1) = sum(ctr(t,2:12),2);
                
                % year of building each power plant in a scenario
                for i=2:12
                    % consider the power generation by fossil fuel infrastructure built by 2025 that has not been been naturally retired
                    pff(t,i,1)=energypath(i,t,cp,1,sce)-enepro(i,1,sce,1)*max(0,(1-(t-1)*5/fflifetime));
                    pff(t,i,2)=energypath(i,t,cp,2,sce)-enepro(i,1,sce,2)*max(0,(1-(t-1)*5/fflifetime));
                    pff(t,i,3)=energypath(i,t,cp,3,sce)-enepro(i,1,sce,3)*max(0,(1-(t-1)*5/fflifetime));
                    
                    % consider the power generation by new fossil fuel infrastructure built after 2025 that has not been been naturally retired
                    for j=1:(fflifetime/5-1)
                        if (t-j)>2
                            pff(t,i,1)=pff(t,i,1)-pff(t-j,i,1); % new built fossil fuel power units in year (t-j)
                            pff(t,i,2)=pff(t,i,2)-pff(t-j,i,2); % new built fossil fuel power units in year (t-j)
                            pff(t,i,3)=pff(t,i,3)-pff(t-j,i,3); % new built fossil fuel power units in year (t-j)
                        end
                    end
                    
                    % consider that a fraction of fossil fuel cannot be substitutued (e.g. aviation)
                    pff(t,i,1)=max(0,pff(t,i,1)-energypath(i,t,cp,1,sce)*frac_unsub); % final energy coal EJ/yr
                    pff(t,i,2)=max(0,pff(t,i,2)-energypath(i,t,cp,2,sce)*frac_unsub); % final energy oil EJ/yr
                    pff(t,i,3)=max(0,pff(t,i,3)-energypath(i,t,cp,3,sce)*frac_unsub); % final energy gas EJ/yr
                    
                    % build power plants using carbon tax revenue
                    decommissions=0;
                    if cp>1 && t>=3
                        % decommissions of old power units after their lifetime
                        for jjjt = 1:6
                            idx=find(punits(:,23)==i & idyear(:,2)==((t-1)*5+2015-pplifetime(jjjt)) & punits(:,2)==jjjt);
                            if size(idx,1)>0
                                idyear(idx,2)=2200; % natural retirement
                            end
                        end
                        
                        % decommissions of old power units that exceed power demand
                        idx=find(punits(:,23)==i & idyear(:,2)<2100);
                        if size(idx,1)>0
                            for j=1:size(idx,1)
                                e=punits(idx(j),2); % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
                                pj=rrr_ori(idx(j))*punits(idx(j),6)*3.6e-3; % final energy TWh/y -> EJ/yr
                                % replace coal, oil and gas sequentially
                                if pff(t,i,1)>pj
                                    pff(t,i,1)=pff(t,i,1)-pj; % replace final energy coal EJ/yr
                                    energypath(i,t,cp,1,sce)=energypath(i,t,cp,1,sce)-pj; % primary energy coal EJ/yr
                                    energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+rrr_ori(idx(j))*punits(idx(j),5); % primary energy EJ/yr
                                else
                                    if pff(t,i,2)>pj
                                        pff(t,i,2)=pff(t,i,2)-pj; % replace final energy oil EJ/yr
                                        energypath(i,t,cp,2,sce)=energypath(i,t,cp,2,sce)-pj; % primary energy oil EJ/yr
                                        energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+rrr_ori(idx(j))*punits(idx(j),5); % primary energy EJ/yr
                                    else
                                        if pff(t,i,3)>pj
                                            pff(t,i,3)=pff(t,i,3)-pj; % replace final energy gas EJ/yr
                                            energypath(i,t,cp,3,sce)=energypath(i,t,cp,3,sce)-pj; % primary energy gas EJ/yr
                                            energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+rrr_ori(idx(j))*punits(idx(j),5); % primary energy EJ/yr
                                        else
                                            % decommissions of power plants due to limits of power demand
                                            decommissions=1;
                                            idyear(idx(j),2)=2200;
                                            continue;
                                        end
                                    end
                                end
                            end
                        end
                        
                        % build new power units using carbon tax revenue in this region
                        %                         idx=find(punits(:,23)==i & idyear(:,1)>2100 & idyear(:,2)>2100);
                        idx=find(punits(:,23)==i & rrr~=0 & idyear(:,2)>2100);
                        if size(idx,1)>0 && ctr(t,i)>0 && decommissions==0
                            j=1;
                            while ctr(t,i)>(punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)))
                                e=punits(idx(j),2); % % 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
                                pj=punits(idx(j),6)*3.6e-3*rrr(idx(j)); % final energy TWh/y -> EJ/yr
                                % replace coal, oil and gas sequentially using power units built in the last time period
                                if pff(t,i,1)>pj
                                    pff(t,i,1)=pff(t,i,1)-pj; % replace final energy coal EJ/yr
                                    energypath(i,t,cp,1,sce)=energypath(i,t,cp,1,sce)-pj; % primary energy coal EJ/yr
                                    energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                    costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                                else
                                    if pff(t,i,2)>pj
                                        pff(t,i,2)=pff(t,i,2)-pj; % replace final energy oil EJ/yr
                                        energypath(i,t,cp,2,sce)=energypath(i,t,cp,2,sce)-pj; % primary energy oil EJ/yr
                                        energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                        costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                                    else
                                        if pff(t,i,3)>pj
                                            pff(t,i,3)=pff(t,i,3)-pj; % replace final energy gas EJ/yr
                                            energypath(i,t,cp,3,sce)=energypath(i,t,cp,3,sce)-pj; % primary energy gas EJ/yr
                                            energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                            costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                                        else
                                            break; % cancel all remaining power units due to the limit of power demand
                                        end
                                    end
                                end
                                idyear(idx(j),2)=(t)*5+2015; % year of building this power unit
                                ctr(t,i)=ctr(t,i)-punits(idx(j),21)*rrr(idx(j))*costlearning(e); % consume carbon tax revenue million $
                                rrr(idx(j))=0;
                                if j<size(idx,1)
                                    j=j+1;
                                else
                                    break;
                                end
                            end
                        end
                    end
                end
                
                % allowing carbon tax revenue to be distributed from one region to other regions
                if financeexchange==1 && cp>0 && t>=3
                    %                     ctr(t,1)=ctr(t-1,1)+sum(ctr(t,2:12),2); % the remaining carbon tax revenue becomes global finance million $
                    ctr(t,1)=sum(ctr(t,2:12),2); % the remaining carbon tax revenue becomes global finance million $
                    ctr(t,2:12)=ctr(t,2:12).*0; % the remaining carbon tax revenue becomes zero in each region
                    %                     idx=find(idyear(:,1)>2100 & idyear(:,2)>2100);
                    idx=find(rrr~=0 & idyear(:,2)>2100);
                    % find the power units that have been built using carbon tax revenue in a region
                    if size(idx,1)>0 && ctr(t,1)>0
                        j=1;
                        while ctr(t,1)>(punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)))
                            i=punits(idx(j),23);
                            if i==0
                                if j<size(idx,1)
                                    j=j+1; % skip this power unit
                                    continue;
                                else
                                    break;
                                end
                            end
                            e=punits(idx(j),2); % % 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
                            pj=punits(idx(j),6)*3.6e-3*rrr(idx(j)); % final energy TWh/y -> EJ/yr
                            % replace coal, oil and gas sequentially
                            if pff(t,i,1)>pj
                                pff(t,i,1)=pff(t,i,1)-pj; % replace final energy coal EJ/yr
                                energypath(i,t,cp,1,sce)=energypath(i,t,cp,1,sce)-pj; % primary energy coal EJ/yr
                                energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                            else
                                if pff(t,i,2)>pj
                                    pff(t,i,2)=pff(t,i,2)-pj; % replace final energy oil EJ/yr
                                    energypath(i,t,cp,2,sce)=energypath(i,t,cp,2,sce)-pj; % primary energy oil EJ/yr
                                    energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                    costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                                else
                                    if pff(t,i,3)>pj
                                        pff(t,i,3)=pff(t,i,3)-pj; % replace final energy gas EJ/yr
                                        energypath(i,t,cp,3,sce)=energypath(i,t,cp,3,sce)-pj; % primary energy gas EJ/yr
                                        energypath(i,t,cp,e+3,sce)=energypath(i,t,cp,e+3,sce)+punits(idx(j),5)*rrr(idx(j)); % primary energy EJ/yr
                                        costREpath(i,t,cp,e+3,sce)=costREpath(i,t,cp,e+3,sce)+punits(idx(j),21)*rrr(idx(j))*costlearning(punits(idx(j),2)); % primary energy EJ/yr
                                    else
                                        if j<size(idx,1)
                                            j=j+1; % skip this power unit
                                            continue;
                                        else
                                            break;
                                        end
                                    end
                                end
                            end
                            % build this power unit by replacing coal, oil or gas
                            idyear(idx(j),2)=(t-1)*5+2015; % year of building this power unit by carbon tax revenue
                            ctr(t,1)=ctr(t,1)-(punits(idx(j),21)*rrr(idx(j))*costlearning(e)); % remaining carbon tax revenue million $
                            rrr(idx(j))=0;
                            if j<size(idx,1)
                                j=j+1;
                            else
                                break;
                            end
                        end
                    end
                end
                
                % determine co2 emissions abated by power plants built using carbon tax revenue
                for i=2:12
                    totalco2=0;
                    totalomcost=0;
                    if cp>1 && t>=3
                        idx=find(punits(:,23)==i & idyear(:,2)<2100);
                        if size(idx,1)>0
                            totalco2=sum(punits(idx,19).*rrr_ori(idx)); % Mt CO2/yr
                            if option_omcost==1
                                totalomcost=sum(punits(idx,22).*rrr_ori(idx)); % total O&M cost million $/y
                            end
                        end
                    end
                    % determine carbon tax revenue based on CO2 emissions and carbon prices
                    %                     co2path(i,t,cp,sce)=enepro(i,t,sce,11)/1e3-totalco2/1e3; % Gt CO2/yr
                    co2path(i,t,cp,sce)=co2path(i,t,cp,sce)-totalco2/1e3; % Gt CO2/yr
                    co2path(1,t,cp,sce)=co2path(1,t,cp,sce)-totalco2/1e3; % Gt CO2/yr
                    %                     if carbonpricetype==1
                    %                         ctr(t,i)=ctr(t-1,i)+max(0,((cp-1)*10*co2path(i,t,cp,sce)*1000-totalomcost)*5); % carbon tax revenue Gt CO2/yr * $/tCO2 * 5 years * 1000 -> million $
                    %                     else if carbonpricetype==2
                    %                             cpc=enepro(i,1,sce,12)*(cp-1)*2/enepro(i,1,sce,11); % GDP|PPP, billion US$2010/yr
                    %                             ctr(t,i)=ctr(t-1,i)+max(0,(cpc*co2path(i,t,cp,sce)*1000-totalomcost)*5); % carbon tax revenue Gt CO2/yr * $/tCO2 * 5 years * 1000 -> million $
                    %                         else if carbonpricetype==3
                    %                                 if i==2 || i==8 || i==7 || i==10
                    %                                     ctr(t,i)=ctr(t-1,i)+max(0,(100*co2path(i,t,cp,sce)*1000-totalomcost)*5); % carbon tax revenue Gt CO2/yr * $/tCO2 * 5 years * 1000 -> million $
                    %                                 else
                    %                                     if t>=5
                    %                                         ctr(t,i)=ctr(t-1,i)+max(0,(100*co2path(i,t,cp,sce)*1000-totalomcost)*5); % carbon tax revenue Gt CO2/yr * $/tCO2 * 5 years * 1000 -> million $
                    %                                     end
                    %                                 end
                    %
                    %                             end
                    %                         end
                    %                     end
                end
            end
        end
        for j=17 % 1:16
            co2path_sce(sce,:) = co2path(1,1:17,j,sce); % Gt CO2/yr, 最终需要
            %             invREpath_sce(sce,:) = sum(sum(costREpath(2:12,1:7,j,4:10,sce),4));% million $
            for typee = 4:1:10
                invREpath_sce_reg(:,:,sce) = invREpath_sce_reg(:,:,sce)+cumsum((costREpath(2:12,1:17,j,typee,sce))'/pplifetime(typee-3)); % 最终需要
                invREiiasa_sce_reg(:,:,sce) = invREiiasa_sce_reg(:,:,sce)+cumsum(costRE_iiasa_reg(1:17,2:12,typee,sce)/pplifetime(typee-3)); % 最终需要
                
                invREpath_sce(sce,:) = invREpath_sce(sce,:)+cumsum(sum(costREpath(2:12,1:17,j,typee,sce),1)/pplifetime(typee-3)); % 最终需要
                invREiiasa_sce(sce,:) = invREiiasa_sce(sce,:)+cumsum(costRE_iiasa(1:17,typee,sce)'/pplifetime(typee-3)); % 最终需要
                invREpath_sce_ori(sce,:) = invREpath_sce_ori(sce,:)+sum(costREpath(2:12,1:17,j,typee,sce),1);
                invREiiasa_sce_ori(sce,:) = invREiiasa_sce_ori(sce,:)+costRE_iiasa(1:17,typee,sce)';
            end
            
            
            energyREpath_sce(sce,:) = sum(sum(energypath(2:12,1:17,j,4:10,sce),4))*1000/3.6;% EJ/yr→TWh/y % 最终需要
            energyALLpath_sce(sce,:) = sum(sum(energypath(2:12,1:17,j,1:10,sce),4))*1000/3.6;% EJ/yr→TWh/y % 最终需要
        end
        sce
    clearvars -except invREpath_sce_reg invREiiasa_sce_reg exp energypath co2path costREpath costRE_iiasa costRE_iiasa_reg costRE_iiasa costREpath costRE_iiasa_reg enepro22 invREpath_sce invREiiasa_sce_ori invREpath_sce_ori invREiiasa_sce invREpath_sce co2path_sce energyREpath_sce energyALLpath_sce
    end
end
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREiiasa_sce_reg_finalenergy0424_AGR.mat','invREiiasa_sce_reg','-v7.3'); % million/y
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREpath_sce_reg_finalenergy0424_AGR.mat','invREpath_sce_reg','-v7.3');

save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREpath_sce_noiiasa_finalenergy0424_AGR.mat','invREpath_sce','-v7.3');
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREiiasa_sce_ori_finalenergy0424_AGR.mat','invREiiasa_sce_ori','-v7.3'); % million/y
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREpath_sce_ori_finalenergy0424_AGR.mat','invREpath_sce_ori','-v7.3');
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREiiasa_sce_finalenergy0424_AGR.mat','invREiiasa_sce','-v7.3'); % million/y
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_invREpath_sce_finalenergy0424_AGR.mat','invREpath_sce','-v7.3');
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_co2path_sce_finalenergy0424_AGR.mat','co2path_sce','-v7.3');
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_energyREpath_sce_finalenergy0424_AGR.mat','energyREpath_sce','-v7.3');
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_energyALLpath_sce_finalenergy0424_AGR.mat','energyALLpath_sce','-v7.3');
energyREpath_sce_iiasa = (reshape(sum(sum(enepro22(2:12,1:7,:,4:9),1),4),[7,1291]))'/3.6;
save('H:\Transfer of carbon-tax revenue\Ans\fig2_iiasa_energyREpath_sce_iiasa_finalenergy0424_AGR.mat','energyREpath_sce_iiasa','-v7.3');
