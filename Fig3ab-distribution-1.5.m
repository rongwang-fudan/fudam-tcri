%% power consumption distribution in 2050
% 仅有各国前10%的格点的能源需求量有增长
tic
clear;
load('H:\Transfer of carbon-tax revenue\Ans\Power_all_TAX.mat','Power_all_TAX');
Power_reg2050 = Power_all_TAX(31,:)';
% load('H:\world\code\region_ID12_0228.mat'); %
cn192id=load('H:\Transfer of carbon-tax revenue\inputs\cn192id.txt'); % 192x4
region_ID=cn192id(:,3);

load('H:\Transfer of carbon-tax revenue\inputs\Ph_country2020.mat')
Ph_country2050 = Ph_country2020*0;
for i = 2:1:12
    [m,n] = find(region_ID==i);
    Ph_country2050(m) = Ph_country2020(m)./sum(Ph_country2020(m)).*Power_reg2050(i); % TWh/y
end

%
load('H:\Transfer of carbon-tax revenue\inputs\GDP_PPP2015_2.mat'); % 2011US dollar
GDP_PPP2015 = GDP_PPP2015_2;
clear GDP_PPP2015_2
load('H:\Transfer of carbon-tax revenue\inputs\GADM_country120_xz2.mat')
load('H:\Transfer of carbon-tax revenue\inputs\E_120_2020.mat', 'E_120')  % billion kWh / year
E_120_2020 = E_120;

power_country2050more = Ph_country2050-Ph_country2020;
[m_more_cou,n] = find(power_country2050more<0);
power_country2050less = power_country2050more*0;
power_country2050less(m_more_cou) = Ph_country2050(m_more_cou)./Ph_country2020(m_more_cou);
power_country2050more(power_country2050more<0)=0;

a = unique(GADM_country120);
a(1)=[];

E_120_1=E_120_2020;
for i = 1:size(m_more_cou,1)
    [m,n]=find(GADM_country120==m_more_cou(i));
    E_120_1(sub2ind(size(E_120_1), m, n)) = E_120_1(sub2ind(size(E_120_1), m, n))*power_country2050less(m_more_cou(i));
end
for i =1:1:size(a,1)
    [m1,n1]=find(GADM_country120==i);
    [B,IX] = sort(E_120_2020(sub2ind(size(E_120_2020), m1, n1)),'descend');
    m = m1(IX(1:floor(size(m1,1)/10)));
    n = n1(IX(1:floor(size(n1,1)/10)));
    GDP_SUM(i,1)=sum(sum(GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n))));
    GDP_CN2015_120coun=GDP_PPP2015(sub2ind(size(GDP_PPP2015), m, n));
    E_120_1(sub2ind(size(E_120_1), m, n))=E_120_1(sub2ind(size(E_120_1), m, n))+GDP_CN2015_120coun/GDP_SUM(i,1)*power_country2050more(i,1); % Consumption, billion kWh / year
    i
end

E_120_1(find(isnan(E_120_1)==1))=0; % Consumption, billion kWh / year
E_120=E_120_1;
save('H:\Transfer of carbon-tax revenue\Ans\E_120_ele_2500_fumccm.mat', 'E_120', '-v7.3')  % billion kWh / year
sum(sum(E_120))
Power_reg2050(1)

%%
tic
clear;
load('H:\Transfer of carbon-tax revenue\Ans\E_120_ele_2500_fumccm.mat')  % billion kWh / year
load('H:\Transfer of carbon-tax revenue\inputs\GADM_county120_xz2.mat')
load('H:\Transfer of carbon-tax revenue\inputs\GADM_country120_xz2.mat')
power_clean_s = zeros(21600,10800);
for i = 1:192
    County = GADM_county120*0;
    [m,n]=find(GADM_country120==i);
    County(sub2ind(size(County), m, n))=GADM_county120(sub2ind(size(GADM_county120), m, n));
    County2 = County(min(m):max(m),min(n):max(n));
    %     [mmm,nnn]=find(power_cy1(:,i)~=0);
    mmm=unique(County2);
    mmm(mmm==0)=[];
    if ~isempty(mmm)
        for i2 = 1:size(mmm,1)
            [m2,n2]=find(County2==mmm(i2));
            Power_county2050(i,mmm(i2)) = sum(sum(E_120(sub2ind(size(E_120), min(m)-1+m2, min(n)-1+n2))));
        end
    end
    i
end

save('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_fumccm.mat','Power_county2050'); % TWh/y

%%
tic
clear;
load('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_fumccm.mat','Power_county2050'); % TWh/y

load('H:\Transfer of carbon-tax revenue\inputs\GADM_county120_xz2.mat')
load('H:\Transfer of carbon-tax revenue\inputs\GADM_country120_xz2.mat')
Power_county2050_s = zeros(21600,10800);
for i = 1:192
    County = GADM_county120*0;
    [m,n]=find(GADM_country120==i);
    County(sub2ind(size(County), m, n))=GADM_county120(sub2ind(size(GADM_county120), m, n));
    County2 = County(min(m):max(m),min(n):max(n));
    %     mmm=unique(County2);
    %     mmm(mmm==0)=[];
    [nnnn,mmm] = find(Power_county2050(i,:)>0);
    if ~isempty(mmm)
        for i2 = 1:size(mmm,2)
            [m2,n2]=find(County2==mmm(i2));
            Power_county2050_s(sub2ind(size(Power_county2050_s), min(m)-1+m2, min(n)-1+n2)) = Power_county2050(i,mmm(i2));
        end
    end
    i
end
save('H:\Transfer of carbon-tax revenue\Ans\Power_county2050_s.mat','Power_county2050_s','-v7.3'); % TWh/y


