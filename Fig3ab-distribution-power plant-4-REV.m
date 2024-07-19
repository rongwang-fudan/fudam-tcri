%% TAX
tic
clear
load('H:\Transfer of carbon-tax revenue\Ans\fig3_punits.dat','-mat');
% 2.dynamic year; 3.type; 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2
load('H:\Transfer of carbon-tax revenue\Ans\idxxx_2020to2200_TAX.mat','idxxx_2020to2200');
% 1.year; 2.reg; 3.type; 4.id
idxxx_2020to2200_TAX = idxxx_2020to2200;
% [m,n]=find(idxxx_2020to2200_TAX(:,1) == 2020);
% idxxx_2020to2200_TAX(m,:) = [];
[m,n]=find(idxxx_2020to2200_TAX(:,1)>2050);
idxxx_2020to2200_TAX(m,:) = [];
idxxx_2020to2200_TAX(:,1) = [];
idxxx_2020to2200_TAX = unique(idxxx_2020to2200_TAX,'rows','stable');

load('H:\Transfer of carbon-tax revenue\Ans\idxxx_2020to2200_CUR.mat','idxxx_2020to2200');
% 1.year; 2.reg; 3.type; 4.id
idxxx_2020to2200_CUR = idxxx_2020to2200;
% [m,n]=find(idxxx_2020to2200_CUR(:,1) == 2020);
% idxxx_2020to2200_CUR(m,:) = [];
[m,n]=find(idxxx_2020to2200_CUR(:,1)>2050);
idxxx_2020to2200_CUR(m,:) = [];
idxxx_2020to2200_CUR(:,1) = [];
idxxx_2020to2200_CUR = unique(idxxx_2020to2200_CUR,'rows','stable');
% 1.reg; 2.type; 3.id
clear idxxx_2020to2200

% sum(punits(idxxx_2020to2200_TAX(:,3),5)/3.6/1000)
% sum(punits(idxxx_2020to2200_CUR(:,3),5)/3.6/1000)

punits_TAX = punits(idxxx_2020to2200_TAX(:,3),:);
[m,n] = find(punits_TAX(:,5)/3.6>10);
[m,n] = find(punits_TAX(:,5)/3.6>10 & punits_TAX(:,3)>2);
sum(punits_TAX(:,15))*100/10^6 % Mha

% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2
for reg = 2:12
    [m,n] = find(punits_TAX(:,4)==reg);
    aaaaa0(reg,1) = sum(punits_TAX(m,5))/3.6;
end
load('H:\Transfer of carbon-tax revenue\Ans\Power_renewable_TAX.mat','Power_renewable_TAX');
aaaaa = Power_renewable_TAX(31,:)';
r = aaaaa(2:12)./aaaaa0(2:12);
for reg = 2:12
    [m,n] = find(punits_TAX(:,4)==reg);
    punits_TAX(m,5) = punits_TAX(m,5)*r(reg-1);
end


load('H:\Transfer of carbon-tax revenue\Ans\Power_all_TAX.mat','Power_all_TAX');
% load('H:\world\code\region_ID12_0228.mat'); %
cn192id=load('inputs\cn192id.txt'); % 192x4
region_ID=cn192id(:,3);
load('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_fumccm.mat','Power_county2050'); % TWh/y
for iii = 2:1:12
    [i2,n] = find(region_ID(:,1)==iii);
    ttti(iii,1) = sum(sum(Power_county2050(i2,:)));
end
[ttti Power_all_TAX(31,:)']
Power_renewable_TAX(31,:)'./Power_all_TAX(31,:)'



%
[m,n]=find(punits_TAX(:,3)<=2);
punits_pw = punits_TAX(m,:);
[m,n]=find(punits_TAX(:,3)>2);
punits_others = punits_TAX(m,:);

[m,n]=find(punits_TAX(:,3)==3);
punits_bioenergy = punits_TAX(m,:);
[m,n]=find(punits_TAX(:,3)==4);
punits_nuclear = punits_TAX(m,:);
[m,n]=find(punits_TAX(:,3)==5);
punits_hydropower = punits_TAX(m,:);
[m,n]=find(punits_TAX(:,3)==6);
punits_geothermal = punits_TAX(m,:);
[m,n]=find(punits_TAX(:,3)==7);
punits_CCS = punits_TAX(m,:);
% 12.county; 13.country; 14.CP(MW);

load('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_fumccm.mat','Power_county2050'); % TWh/y
load('H:\Transfer of carbon-tax revenue\inputs\geothermal2020.dat','-mat');
[m,n] = find(geothermal2020(:,1)>0 & geothermal2020(:,2)>0);
geothermal2020(m,5)= Power_county2050(sub2ind(size(Power_county2050), geothermal2020(m,1), geothermal2020(m,2)));
[B,IX] = sort(geothermal2020(:,5),'descend');
geothermal2020_IX = geothermal2020(IX,:);
clear geothermal2020
load('H:\Transfer of carbon-tax revenue\inputs\nuclear2020.dat','-mat');
[m,n] = find(nuclear2020(:,1)>0 & nuclear2020(:,2)>0);
nuclear2020(m,5)= Power_county2050(sub2ind(size(Power_county2050), nuclear2020(m,1), nuclear2020(m,2)));
[B,IX] = sort(nuclear2020(:,5),'descend');
nuclear2020_IX = nuclear2020(IX,:);
clear nuclear2020
load('H:\Transfer of carbon-tax revenue\inputs\bioenergy2020.dat','-mat');
[m,n] = find(bioenergy2020(:,1)>0 & bioenergy2020(:,2)>0);
bioenergy2020(m,5)= Power_county2050(sub2ind(size(Power_county2050), bioenergy2020(m,1), bioenergy2020(m,2)));
[B,IX] = sort(bioenergy2020(:,5),'descend');
bioenergy2020_IX = bioenergy2020(IX,:);
clear bioenergy2020
load('H:\Transfer of carbon-tax revenue\inputs\hydropower2020.dat','-mat');
[m,n] = find(hydropower2020(:,1)>0 & hydropower2020(:,2)>0);
hydropower2020(m,5)= Power_county2050(sub2ind(size(Power_county2050), hydropower2020(m,1), hydropower2020(m,2)));
[B,IX] = sort(hydropower2020(:,5),'descend');
hydropower2020_IX = hydropower2020(IX,:);
clear hydropower2020
% 1.Country; 2.County; 3.Power (TWh/y); 4.reg; 5.power demand of this county (TWh/y); 6.latitude; 7.longitude

% punits_bioenergy_xz
plant_new = [];
for i = 2:1:12
    plant_new2= [];
    [m2020,n2020] = find(bioenergy2020_IX(:,4)==i);
    plantaa2020 = bioenergy2020_IX(m2020,:);
    % 1.Country; 2.County; 3.Power (TWh/y); 4.reg; 5.power demand of this county (TWh/y); 6.latitude; 7.longitude
    [m,n] = find(punits_bioenergy(:,4)==i);
    plantaa = punits_bioenergy(m,:);
    plantaa_o = plantaa;
    if ~isempty(m)
        if sum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6)
            [m2,n2] = find(cumsum(plantaa(:,5)/3.6)<=sum(plantaa2020(:,3)));% TWh/y
            aaaw = sum(plantaa2020(:,3))-sum(plantaa(m2,5)/3.6); % 还可以指定位置的电厂
            plant_new2(:,13) = plantaa2020(:,1);
            plant_new2(:,15) =0;
            plant_new2(:,12) = plantaa2020(:,2);
            plant_new2(:,5) = plantaa2020(:,3)*3.6;
            plant_new2(:,4) = plantaa2020(:,4);
            plant_new2(:,9) = plantaa2020(:,6);
            plant_new2(:,10) = plantaa2020(:,7);
            plant_new2(:,3) = plantaa(1,3);
            plant_new2(:,11) = mean(plantaa(m2,11));
            plantaa(m2,:)=[];
            plantaa(1,5) = plantaa(1,5)-aaaw*3.6;
            plant_new = [plant_new;plant_new2;plantaa];
        else if sum(plantaa2020(:,3))>sum(plantaa(:,5)/3.6)
                [m2,n2] = find(cumsum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6));% TWh/y
                aaaaw = sum(plantaa(:,5))-sum(plantaa2020(m2,3)*3.6);
                if ~isempty(m2)
                    m2 = [m2;m2(end)+1];
                    plant_new2(:,13) = plantaa2020(m2,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(m2,2);
                    plant_new2(:,5) = plantaa2020(m2,3)*3.6;
                    plant_new2(:,4) = plantaa2020(m2,4);
                    plant_new2(:,9) = plantaa2020(m2,6);
                    plant_new2(:,10) = plantaa2020(m2,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new2(end,5) = aaaaw;
                    plant_new = [plant_new;plant_new2];
                else
                    plant_new2(:,13) = plantaa2020(1,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(1,2);
                    plant_new2(:,5) = sum(plantaa(:,5));
                    plant_new2(:,4) = plantaa2020(1,4);
                    plant_new2(:,9) = plantaa2020(1,6);
                    plant_new2(:,10) = plantaa2020(1,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new = [plant_new;plant_new2];
                end
            end
        end
    end
end
punits_bioenergy_xz=plant_new;
sum(punits_bioenergy_xz(:,5))/3.6
sum(punits_bioenergy(:,5))/3.6
clear punits_bioenergy

% punits_nuclear_xz
plant_new = [];
for i = 2:1:12
    plant_new2= [];
    [m2020,n2020] = find(nuclear2020_IX(:,4)==i);
    plantaa2020 = nuclear2020_IX(m2020,:);
    % 1.Country; 2.County; 3.Power (TWh/y); 4.reg; 5.power demand of this county (TWh/y)
    [m,n] = find(punits_nuclear(:,4)==i);
    plantaa = punits_nuclear(m,:);
    plantaa_o = plantaa;
    if ~isempty(m)
        if sum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6)
            [m2,n2] = find(cumsum(plantaa(:,5)/3.6)<=sum(plantaa2020(:,3)));% TWh/y
            aaaw = sum(plantaa2020(:,3))-sum(plantaa(m2,5)/3.6); % 还可以指定位置的电厂
            plant_new2(:,13) = plantaa2020(:,1);
            plant_new2(:,15) =0;
            plant_new2(:,12) = plantaa2020(:,2);
            plant_new2(:,5) = plantaa2020(:,3)*3.6;
            plant_new2(:,4) = plantaa2020(:,4);
            plant_new2(:,9) = plantaa2020(:,6);
            plant_new2(:,10) = plantaa2020(:,7);
            plant_new2(:,3) = plantaa(1,3);
            plant_new2(:,11) = mean(plantaa(m2,11));
            plantaa(m2,:)=[];
            plantaa(1,5) = plantaa(1,5)-aaaw*3.6;
            plant_new = [plant_new;plant_new2;plantaa];
        else if sum(plantaa2020(:,3))>sum(plantaa(:,5)/3.6)
                [m2,n2] = find(cumsum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6));% TWh/y
                aaaaw = sum(plantaa(:,5))-sum(plantaa2020(m2,3)*3.6);
                if ~isempty(m2)
                    m2 = [m2;m2(end)+1];
                    plant_new2(:,13) = plantaa2020(m2,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(m2,2);
                    plant_new2(:,5) = plantaa2020(m2,3)*3.6;
                    plant_new2(:,4) = plantaa2020(m2,4);
                    plant_new2(:,9) = plantaa2020(m2,6);
                    plant_new2(:,10) = plantaa2020(m2,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new2(end,5) = aaaaw;
                    plant_new = [plant_new;plant_new2];
                else
                    plant_new2(:,13) = plantaa2020(1,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(1,2);
                    plant_new2(:,5) = sum(plantaa(:,5));
                    plant_new2(:,4) = plantaa2020(1,4);
                    plant_new2(:,9) = plantaa2020(1,6);
                    plant_new2(:,10) = plantaa2020(1,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new = [plant_new;plant_new2];
                end
            end
        end
    end
end
punits_nuclear_xz=plant_new;
clear punits_nuclear

% punits_hydropower_xz
plant_new = [];
for i = 2:1:12
    plant_new2= [];
    [m2020,n2020] = find(hydropower2020_IX(:,4)==i);
    plantaa2020 = hydropower2020_IX(m2020,:);
    % 1.Country; 2.County; 3.Power (TWh/y); 4.reg; 5.power demand of this county (TWh/y)
    [m,n] = find(punits_hydropower(:,4)==i);
    plantaa = punits_hydropower(m,:);
    plantaa_o = plantaa;
    if ~isempty(m)
        if sum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6)
            [m2,n2] = find(cumsum(plantaa(:,5)/3.6)<=sum(plantaa2020(:,3)));% TWh/y
            aaaw = sum(plantaa2020(:,3))-sum(plantaa(m2,5)/3.6); % 还可以指定位置的电厂
            plant_new2(:,13) = plantaa2020(:,1);
            plant_new2(:,15) =0;
            plant_new2(:,12) = plantaa2020(:,2);
            plant_new2(:,5) = plantaa2020(:,3)*3.6;
            plant_new2(:,4) = plantaa2020(:,4);
            plant_new2(:,9) = plantaa2020(:,6);
            plant_new2(:,10) = plantaa2020(:,7);
            plant_new2(:,3) = plantaa(1,3);
            plant_new2(:,11) = mean(plantaa(m2,11));
            plantaa(m2,:)=[];
            plantaa(1,5) = plantaa(1,5)-aaaw*3.6;
            plant_new = [plant_new;plant_new2;plantaa];
        else if sum(plantaa2020(:,3))>sum(plantaa(:,5)/3.6)
                [m2,n2] = find(cumsum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6));% TWh/y
                aaaaw = sum(plantaa(:,5))-sum(plantaa2020(m2,3)*3.6);
                if ~isempty(m2)
                    m2 = [m2;m2(end)+1];
                    plant_new2(:,13) = plantaa2020(m2,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(m2,2);
                    plant_new2(:,5) = plantaa2020(m2,3)*3.6;
                    plant_new2(:,4) = plantaa2020(m2,4);
                    plant_new2(:,9) = plantaa2020(m2,6);
                    plant_new2(:,10) = plantaa2020(m2,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new2(end,5) = aaaaw;
                    plant_new = [plant_new;plant_new2];
                else
                    plant_new2(:,13) = plantaa2020(1,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(1,2);
                    plant_new2(:,5) = sum(plantaa(:,5));
                    plant_new2(:,4) = plantaa2020(1,4);
                    plant_new2(:,9) = plantaa2020(1,6);
                    plant_new2(:,10) = plantaa2020(1,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new = [plant_new;plant_new2];
                end
            end
        end
    end
end
punits_hydropower_xz=plant_new;
sum(punits_hydropower_xz(:,5))/3.6
sum(punits_hydropower(:,5))/3.6
clear punits_hydropower

% punits_geothermal_xz
plant_new = [];
for i = 2:1:12
    plant_new2= [];
    [m2020,n2020] = find(geothermal2020_IX(:,4)==i);
    plantaa2020 = geothermal2020_IX(m2020,:);
    % 1.Country; 2.County; 3.Power (TWh/y); 4.reg; 5.power demand of this county (TWh/y)
    [m,n] = find(punits_geothermal(:,4)==i);
    plantaa = punits_geothermal(m,:);
    plantaa_o = plantaa;
    if ~isempty(m)
        if sum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6)
            [m2,n2] = find(cumsum(plantaa(:,5)/3.6)<=sum(plantaa2020(:,3)));% TWh/y
            aaaw = sum(plantaa2020(:,3))-sum(plantaa(m2,5)/3.6); % 还可以指定位置的电厂
            plant_new2(:,13) = plantaa2020(:,1);
            plant_new2(:,15) =0;
            plant_new2(:,12) = plantaa2020(:,2);
            plant_new2(:,5) = plantaa2020(:,3)*3.6;
            plant_new2(:,4) = plantaa2020(:,4);
            plant_new2(:,9) = plantaa2020(:,6);
            plant_new2(:,10) = plantaa2020(:,7);
            plant_new2(:,3) = plantaa(1,3);
            plant_new2(:,11) = mean(plantaa(m2,11));
            plantaa(m2,:)=[];
            plantaa(1,5) = plantaa(1,5)-aaaw*3.6;
            plant_new = [plant_new;plant_new2;plantaa];
        else if sum(plantaa2020(:,3))>sum(plantaa(:,5)/3.6)
                [m2,n2] = find(cumsum(plantaa2020(:,3))<=sum(plantaa(:,5)/3.6));% TWh/y
                aaaaw = sum(plantaa(:,5))-sum(plantaa2020(m2,3)*3.6);
                if ~isempty(m2)
                    m2 = [m2;m2(end)+1];
                    plant_new2(:,13) = plantaa2020(m2,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(m2,2);
                    plant_new2(:,5) = plantaa2020(m2,3)*3.6;
                    plant_new2(:,4) = plantaa2020(m2,4);
                    plant_new2(:,9) = plantaa2020(m2,6);
                    plant_new2(:,10) = plantaa2020(m2,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new2(end,5) = aaaaw;
                    plant_new = [plant_new;plant_new2];
                else
                    plant_new2(:,13) = plantaa2020(1,1);
                    plant_new2(:,15) =0;
                    plant_new2(:,12) = plantaa2020(1,2);
                    plant_new2(:,5) = sum(plantaa(:,5));
                    plant_new2(:,4) = plantaa2020(1,4);
                    plant_new2(:,9) = plantaa2020(1,6);
                    plant_new2(:,10) = plantaa2020(1,7);
                    plant_new2(:,3) = plantaa(1,3);
                    plant_new2(:,11) = mean(plantaa(:,11));
                    plant_new = [plant_new;plant_new2];
                end
            end
        end
    end
end
punits_geothermal_xz=plant_new;
sum(punits_geothermal_xz(:,5))/3.6
sum(punits_geothermal(:,5))/3.6
clear punits_geothermal
% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2

% 有明确位置的电厂
[m1,n] = find(punits_bioenergy_xz(:,13)>0 & punits_bioenergy_xz(:,12)>0);
[m2,n] = find(punits_nuclear_xz(:,13)>0 & punits_nuclear_xz(:,12)>0);
[m3,n] = find(punits_hydropower_xz(:,13)>0 & punits_hydropower_xz(:,12)>0);
[m4,n] = find(punits_geothermal_xz(:,13)>0 & punits_geothermal_xz(:,12)>0);
aaaf = [punits_pw; punits_bioenergy_xz(m1,:);punits_nuclear_xz(m2,:);punits_hydropower_xz(m3,:);punits_geothermal_xz(m4,:)];
rrrryui = [punits_pw; punits_bioenergy_xz;punits_nuclear_xz;punits_hydropower_xz;punits_geothermal_xz;punits_CCS];
[m,n] = find(rrrryui(:,5)/3.6>5);
[m,n] = find(rrrryui(:,5)/3.6>5 & rrrryui(:,3)>2);

for uu = 1:7
    [m,n] = find(rrrryui(:,5)/3.6>5 & rrrryui(:,3)==uu);
    awer(uu,1) = size(m,1);
end

%
load('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_fumccm.mat','Power_county2050'); % TWh/y
Power_generation_county2050 = Power_county2050*0;
for i = 1:size(aaaf,1)
    Power_generation_county2050(aaaf(i,13), aaaf(i,12)) = Power_generation_county2050(aaaf(i,13), aaaf(i,12)) + aaaf(i,5)/3.6; % country
end
power_gap = Power_county2050-Power_generation_county2050;
power_gap(power_gap<0)=0; % 缺口

power_sur = Power_generation_county2050-Power_county2050;
[m,n] = find(power_sur>0);
Power_generation_county2050(sub2ind(size(Power_generation_county2050), m,n)) = Power_county2050(sub2ind(size(Power_county2050), m,n));
power_sur(power_sur<0)=0; % 多余电量
power_sur_cou = sum(power_sur,2);

for i = 1:192
    [m,n1] = find(power_sur(i,:)>0);
    power_sur_cou(i,1) = sum(power_sur(i,:));
    [m,n2] = find(power_gap(i,:)>0); % 缺电
    [B,IX] = sort(Power_county2050(i,n2)','descend');
    for j = 1:size(IX,1)
        if power_sur_cou(i,1)>=power_gap(i,n2(IX(j)))
            Power_generation_county2050(i,n2(IX(j))) = Power_generation_county2050(i,n2(IX(j)))+power_gap(i,n2(IX(j)));
            power_sur_cou(i,1) = power_sur_cou(i,1)-power_gap(i,n2(IX(j)));
            power_gap(i,n2(IX(j)))=0;
        else
            Power_generation_county2050(i,n2(IX(j))) = Power_generation_county2050(i,n2(IX(j)))+power_sur_cou(i,1);
            power_sur_cou(i,1) = power_sur_cou(i,1) - power_sur_cou(i,1);
            power_gap(i,n2(IX(j)))=power_gap(i,n2(IX(j)))-power_sur_cou(i,1);
        end
    end
    i
end

% 9.latitude; 10.longitude;
% 12.county; 13.country
% load('H:\world\code\region_ID12_0228.mat'); %
cn192id=load('H:\Transfer of carbon-tax revenue\inputs\cn192id.txt'); % 192x4
region_ID=cn192id(:,3);
power_sur_reg = zeros(12,1);
for reg = 2:1:12
    [m,n] = find(region_ID==reg);
    power_sur_reg(reg,1) = sum(power_sur_cou(m));
end

power_gap = Power_county2050-Power_generation_county2050;
power_gap(power_gap<0)=0; % 缺口
power_sur = Power_generation_county2050-Power_county2050;
[m,n] = find(power_sur>0);
Power_generation_county2050(sub2ind(size(Power_generation_county2050), m,n)) = Power_county2050(sub2ind(size(Power_county2050), m,n));

for iii = 2:12 % 2:1:12
    [i2,n] = find(region_ID(:,1)==iii);
    [m,n2] = find(power_gap(i2,:)>0); % 缺电
    if ~isempty(m)
        xx = i2(m);
        ASD = reshape(Power_county2050(sub2ind(size(Power_county2050), xx, n2)),[max(size(n2)) 1]);
        [B,IX] = sort(ASD,'descend');
        for j = 1:size(IX,1)
            i = xx(IX(j));
            if power_sur_reg(iii,1)>=power_gap(i,n2(IX(j)))
                Power_generation_county2050(i,n2(IX(j))) = Power_generation_county2050(i,n2(IX(j)))+power_gap(i,n2(IX(j)));
                power_sur_reg(iii,1) = power_sur_reg(iii,1)-power_gap(i,n2(IX(j)));
                power_gap(i,n2(IX(j)))=0;
            else
                Power_generation_county2050(i,n2(IX(j))) = Power_generation_county2050(i,n2(IX(j)))+power_sur_reg(iii,1);
                power_sur_reg(iii,1) = power_sur_reg(iii,1) - power_sur_reg(iii,1);
                power_gap(i,n2(IX(j)))=power_gap(i,n2(IX(j)))-power_sur_reg(iii,1);
            end
        end
    end
    iii
end

[m1,n] = find(punits_bioenergy_xz(:,13)>0 & punits_bioenergy_xz(:,12)>0);
[m2,n] = find(punits_nuclear_xz(:,13)>0 & punits_nuclear_xz(:,12)>0);
[m3,n] = find(punits_hydropower_xz(:,13)>0 & punits_hydropower_xz(:,12)>0);
[m4,n] = find(punits_geothermal_xz(:,13)>0 & punits_geothermal_xz(:,12)>0);
punits_bioenergy_xz(m1,:)=[];
punits_nuclear_xz(m2,:)=[];
punits_hydropower_xz(m3,:)=[];
punits_geothermal_xz(m4,:)=[];
aaaf22 = [punits_bioenergy_xz;punits_nuclear_xz;punits_hydropower_xz;punits_geothermal_xz;punits_CCS];

for reg = 2:1:12
    [m,n] = find(aaaf22(:,4)==reg);
    power_sur_reg(reg,1) = power_sur_reg(reg,1)+sum(aaaf22(m,5))/3.6;
end

power_gap = Power_county2050-Power_generation_county2050;
power_gap(power_gap<0)=0; % 缺口
power_sur = Power_generation_county2050-Power_county2050;
[m,n] = find(power_sur>0);
Power_generation_county2050(sub2ind(size(Power_generation_county2050), m,n)) = Power_county2050(sub2ind(size(Power_county2050), m,n));

% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2
[B,IX] = sort(aaaf22(:,11));
aaaf22_IX = aaaf22(IX,:);
clear aaaf22
aaaf22 = aaaf22_IX;
clear aaaf22_IX

for i2cc = 1:size(aaaf22,1)
    reg = aaaf22(i2cc,4);
    [i2,n] = find(region_ID(:,1)==reg);
    [m,n2] = find(power_gap(i2,:)>0); % 缺电
    if ~isempty(m)
        xx = i2(m);
        ASD = reshape(Power_county2050(sub2ind(size(Power_county2050), xx, n2)),[max(size(n2)) 1]);
        [B,IX] = sort(ASD,'descend');
        i = xx(IX(1));
        countyyy = n2(IX(1));
        aaaf22(i2cc,13) = i;
        aaaf22(i2cc,12) = countyyy;
        power_gap(i,n2(IX(1)))= power_gap(i,n2(IX(1)))-aaaf22(i2cc,5)/3.6;
        power_sur_reg(reg,1) = power_sur_reg(reg,1)-aaaf22(i2cc,5)/3.6;
    end
    i2cc
end
renewable_power_pant_fix = aaaf;
save('H:\Transfer of carbon-tax revenue\Ans\renewable_power_pant_TAX_posifix.mat','renewable_power_pant_fix','-v7.3');

renewable_power_pant = aaaf22;
save('H:\Transfer of carbon-tax revenue\Ans\renewable_power_pant_TAX_no.mat','renewable_power_pant','-v7.3');
% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2

%% 对已确定国家和县级 但没有精确经纬度的电厂确定经纬度
tic
clear;
load('H:\Transfer of carbon-tax revenue\Ans\renewable_power_pant_TAX_no.mat','renewable_power_pant');
% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2
load('H:\Transfer of carbon-tax revenue\inputs\lat_county_1.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lon_county_1.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lat_county_2.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lon_county_2.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lat_county_4.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lon_county_4.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lat_county_16.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lon_county_16.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lat_county_64.mat')
load('H:\Transfer of carbon-tax revenue\inputs\lon_county_64.mat')
% 前64是第一个countyde 64个位置
load('H:\Transfer of carbon-tax revenue\inputs\ID_county_64.mat')
lat_county_1(find(isnan(lat_county_1)==1))=0;
lon_county_1(find(isnan(lon_county_1)==1))=0;
lat_county_2(find(isnan(lat_county_2)==1))=0;
lon_county_2(find(isnan(lon_county_2)==1))=0;
lat_county_4(find(isnan(lat_county_4)==1))=0;
lon_county_4(find(isnan(lon_county_4)==1))=0;
lat_county_16(find(isnan(lat_county_16)==1))=0;
lon_county_16(find(isnan(lon_county_16)==1))=0;
lat_county_64(find(isnan(lat_county_64)==1))=0;
lon_county_64(find(isnan(lon_county_64)==1))=0;

renewable_power_pant_a = renewable_power_pant;
clear renewable_power_pant
% 2.dynamic year; 3.type (1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS);
% 4.region id; 5.power generation PJ/yr; 6.total investment t$
% 7.abated emission GtCO2/yr; 8.fraction of non-land costs in total costs;
% 9.latitude; 10.longitude; 11.marginal abatement cost $/tCO2;
% 12.county; 13.country; 14.CP(MW); 15.km2

cou_p_num = zeros(size(lat_county_1));
for i = 1:size(renewable_power_pant_a,1)
    county = renewable_power_pant_a(i,12);
    country = renewable_power_pant_a(i,13);
    cou_p_num(county,country) = cou_p_num(county,country)+ 1;
end

[mc,nc] = find(cou_p_num>0);
for i = 1:size(mc,1)
    [m,n] = find(renewable_power_pant_a(:,12)==mc(i) & renewable_power_pant_a(:,13)==nc(i));
    aq1 = [lat_county_1(mc(i),nc(i)),lon_county_1(mc(i),nc(i))];
    aq2 = [lat_county_2(mc(i)*2-1:mc(i)*2,nc(i)),lon_county_2(mc(i)*2-1:mc(i)*2,nc(i))];
    aq4 = [lat_county_4(mc(i)*4-3:mc(i)*4,nc(i)),lon_county_4(mc(i)*4-3:mc(i)*4,nc(i))];
    aq16 = [lat_county_16(mc(i)*16-15:mc(i)*16,nc(i)),lon_county_16(mc(i)*16-15:mc(i)*16,nc(i))];
    [mkk,nkk] = find(floor(ID_county_64(:,nc(i)))==mc(i));
    aq64 = [lat_county_64(mkk,nc(i)), lon_county_64(mkk,nc(i))];
    [mko,nk0] = find(aq1(:,1)==0 | aq1(:,2)==0);
    aq1(mko,:) = [];
    [mko,nk0] = find(aq2(:,1)==0 | aq2(:,2)==0);
    aq2(mko,:) = [];
    [mko,nk0] = find(aq4(:,1)==0 | aq4(:,2)==0);
    aq4(mko,:) = [];
    [mko,nk0] = find(aq16(:,1)==0 | aq16(:,2)==0);
    aq16(mko,:) = [];
    [mko,nk0] = find(aq64(:,1)==0 | aq64(:,2)==0);
    aq64(mko,:) = [];
    if size(aq64,1)>=size(aq4,1)
        aq6422 = [];
        for uuy = 1:1000
            aq6422 = [aq6422;aq64];
        end
        aq64 = aq6422;
        clear aq6422
    else
        aq6422 = [];
        for uuy = 1:1000
            aq6422 = [aq6422;aq16];
        end
        aq64 = aq6422;
        clear aq6422
    end
    
    if size(aq1,1)>=size(m,1)
        renewable_power_pant_a(m,9) = 90-aq1(1:size(m,1),1)/120; % latitude
        renewable_power_pant_a(m,10) = aq1(1:size(m,1),2)/30-180;
    else if size(aq2,1)>=size(m,1)
            renewable_power_pant_a(m,9) = 90-aq2(1:size(m,1),1)/120; % latitude
            renewable_power_pant_a(m,10) = aq2(1:size(m,1),2)/30-180;
        else if size(aq4,1)>=size(m,1)
                renewable_power_pant_a(m,9) = 90-aq4(1:size(m,1),1)/120; % latitude
                renewable_power_pant_a(m,10) = aq4(1:size(m,1),2)/30-180;
            else if size(aq16,1)>=size(m,1)
                    renewable_power_pant_a(m,9) = 90-aq16(1:size(m,1),1)/120; % latitude
                    renewable_power_pant_a(m,10) = aq16(1:size(m,1),2)/30-180;
                else if size(aq64,1)>=size(m,1)
                        renewable_power_pant_a(m,9) = 90-aq64(1:size(m,1),1)/120; % latitude
                        renewable_power_pant_a(m,10) = aq64(1:size(m,1),2)/30-180;
                    else if size(aq64,1)<size(m,1)
                            aq6422 = [];
                            for uuy = 1:1000
                                aq6422 = [aq6422;aq64];
                            end
                            aq64 = aq6422;
                            clear aq6422
                            renewable_power_pant_a(m,9) = 90-aq64(1:size(m,1),1)/120; % latitude
                            renewable_power_pant_a(m,10) = aq64(1:size(m,1),2)/30-180;
                        end
                    end
                end
            end
        end
    end
    i
end
[m,n] = find(renewable_power_pant_a(:,10)==-180);
figure; plot(renewable_power_pant_a(m,9))
[m,n] = find(renewable_power_pant_a(:,9)==0 & renewable_power_pant_a(:,10)==0);
[m,n] = find(isnan(renewable_power_pant_a(:,9)));
renewable_power_pant_a(m,5)/3.6;
figure;plot(ans)
renewable_power_pant_TAX_posinofix = renewable_power_pant_a;
save('H:\Transfer of carbon-tax revenue\Ans\renewable_power_pant_TAX_posinofix.mat','renewable_power_pant_TAX_posinofix','-v7.3');



