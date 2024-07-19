% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2024.1.5
tic
clear;
clear global;
% Initialset_econ;
Initialset_econ_test; % idx_punits0为punits中在2020年及之前建设的电厂的id
idxxx_2020to2200(:,2:4) = idx_punits0;
idxxx_2020to2200(:,1) = 2020;
Initialset_clim;
fracinv0=[0 300 52 45 39 19 325 175 10 14 10 10]; % Figure 20 in Global Landscape of Climate Finance 2023, Climate Policy Initiative, https://www.climatepolicyinitiative.org/publication/global-landscape-of-climate-finance-2023/
fracinv0(1,1)=sum(fracinv0(1,2:12),2); fracinv0=fracinv0./fracinv0(1,1);
inv15=[0.340 0.263 0.351 0.322 0.329 0.348 0.430 0.499; 75.28 76.52 81.48 86.54 87.78 85.26 97.53 101.33; 0.00452 0.00344 0.00431 0.00372 0.00375 0.00408 0.00441 0.00492]; % 1 renewable energy investment t$; 2 gdp; 3 ratio from 2015 to 2022
t15=[2015:2022]; [sR,lrv0,bb0] = regression(t15,inv15(3,:));

% set-up of experiments
fcpc=1; % fraction of carbon tax revenue that is used for renewable energy investments
tinv=1; % transfer of carbon tax revenue between regions: 0 no; 1 yes
lrv=lrv0; % rate of growth in renewable energy investment
cpc0=100; % carbon price $/tCO2

% simulations
Ty2=2200;
S=zeros(Ty2-Ty0+1,Ne*12+46); % ne economic variables for 12 regions + 46 climatic variables for the world
S(1,1:(Ne*12))=econo0(1,1:(Ne*12)); % economic variable in 2020
S(1,(Ne*12+1):(Ne*12+11))=clim0(1,1:11); % climatic variable in 2020
for t=1:(Ty2-Ty0) % 2020 to 2303
    tyear=t+Ty0-1;
    display(tyear);
    for i=2:12
        i2=Ne*(i-1);
        ffout=0;
        for j=(t-lifetimeff(2)):(t-lifetimeff(1))
            ffout=ffout+S(max(j,1),i2+20)/(lifetimeff(2)-lifetimeff(1)+1); % fossil-fuel retired PJ/yr
        end
        if i==2 || i==7 || i==8 || i==10
            cpc=cpc0*max(0,min(1,(tyear-2024)/5)); % carbon price $/tCO2
        else
            cpc=cpc0*max(0,min(1,(tyear-2035)/5)); % carbon price $/tCO2
        end
        if tyear<=2022
            finv=inv15(1,tyear-2014)*fracinv0(i)/S(tyear-2019,i2+7); % climate investments % of GDP
        else
            finv=inv15(1,2022-2014)*fracinv0(i)/S(2022-2019,i2+7)+lrv*(min(2100,tyear)-2022);
        end
        S(t+1,(i2+1):(i2+Ne))=econdyn(t+1,i,S(t,(i2+1):(i2+Ne)),finv,cpc,fcpc,ffout,S(t,8+Ne*12),S(t,46+Ne*12));
    end
    %     S(t+1,1:(Ne*12))=invdyn(t+1,S(t+1,1:(Ne*12)));
    clear idxxx
    [S(t+1,1:(Ne*12)),idxxx]=invdyn_test(t+1,S(t+1,1:(Ne*12)));
    idxxx_2020to2200 = [idxxx_2020to2200;[ones(size(idxxx,1),1)*(t+2020) idxxx]];
    emico2=S(t+1,17); % emissions GtCO2
    rff=(S(t+1,12)-sum(S(t+1,44:49),2))/(S(1,12)-sum(S(1,44:49),2)); % ratio of fossil fuel change
    S(t+1,(Ne*12+1):(Ne*12+46))=climdyn(t+1, S(t,(Ne*12+1):(Ne*12+46)), emico2, rff );
end
% idxxx_2020to2200
% 1.year; 2.reg; 3.type; 4.id
save('H:\Transfer of carbon-tax revenue\Ans\idxxx_2020to2200_TAX.mat','idxxx_2020to2200');

Power_all_TAX = S(:,Ne*([1:12]-1)+12)./3600*1000; % energy TWh/y
Power_renewable_TAX = S(:,Ne*([1:12]-1)+21)./3600*1000; % reneable energy TWh/y
save('H:\Transfer of carbon-tax revenue\Ans\Power_all_TAX.mat','Power_all_TAX');
save('H:\Transfer of carbon-tax revenue\Ans\Power_renewable_TAX.mat','Power_renewable_TAX');


%%
tic
clear;
clear global;
% Initialset_econ;
Initialset_econ_test; % idx_punits0为punits中在2020年及之前建设的电厂的id
idxxx_2020to2200(:,2:4) = idx_punits0;
idxxx_2020to2200(:,1) = 2020;
Initialset_clim;
fracinv0=[0 300 52 45 39 19 325 175 10 14 10 10]; % Figure 20 in Global Landscape of Climate Finance 2023, Climate Policy Initiative, https://www.climatepolicyinitiative.org/publication/global-landscape-of-climate-finance-2023/
fracinv0(1,1)=sum(fracinv0(1,2:12),2); fracinv0=fracinv0./fracinv0(1,1);
inv15=[0.340 0.263 0.351 0.322 0.329 0.348 0.430 0.499; 75.28 76.52 81.48 86.54 87.78 85.26 97.53 101.33; 0.00452 0.00344 0.00431 0.00372 0.00375 0.00408 0.00441 0.00492]; % 1 renewable energy investment t$; 2 gdp; 3 ratio from 2015 to 2022
t15=[2015:2022]; [sR,lrv0,bb0] = regression(t15,inv15(3,:));

% set-up of experiments
fcpc=1; % fraction of carbon tax revenue that is used for renewable energy investments
tinv=1; % transfer of carbon tax revenue between regions: 0 no; 1 yes
lrv=lrv0; % rate of growth in renewable energy investment
cpc0=0; % carbon price $/tCO2

% simulations
Ty2=2200;
S=zeros(Ty2-Ty0+1,Ne*12+46); % ne economic variables for 12 regions + 46 climatic variables for the world
S(1,1:(Ne*12))=econo0(1,1:(Ne*12)); % economic variable in 2020
S(1,(Ne*12+1):(Ne*12+11))=clim0(1,1:11); % climatic variable in 2020
for t=1:(Ty2-Ty0) % 2020 to 2303
    tyear=t+Ty0-1;
    display(tyear);
    for i=2:12
        i2=Ne*(i-1);
        ffout=0;
        for j=(t-lifetimeff(2)):(t-lifetimeff(1))
            ffout=ffout+S(max(j,1),i2+20)/(lifetimeff(2)-lifetimeff(1)+1); % fossil-fuel retired PJ/yr
        end
        if i==2 || i==7 || i==8 || i==10
            cpc=cpc0*max(0,min(1,(tyear-2024)/5)); % carbon price $/tCO2
        else
            cpc=cpc0*max(0,min(1,(tyear-2035)/5)); % carbon price $/tCO2
        end
        if tyear<=2022
            finv=inv15(1,tyear-2014)*fracinv0(i)/S(tyear-2019,i2+7); % climate investments % of GDP
        else
            finv=inv15(1,2022-2014)*fracinv0(i)/S(2022-2019,i2+7)+lrv*(min(2100,tyear)-2022);
        end
        S(t+1,(i2+1):(i2+Ne))=econdyn(t+1,i,S(t,(i2+1):(i2+Ne)),finv,cpc,fcpc,ffout,S(t,8+Ne*12),S(t,46+Ne*12));
    end
    %     S(t+1,1:(Ne*12))=invdyn(t+1,S(t+1,1:(Ne*12)));
    clear idxxx
    [S(t+1,1:(Ne*12)),idxxx]=invdyn_test(t+1,S(t+1,1:(Ne*12)));
    idxxx_2020to2200 = [idxxx_2020to2200;[ones(size(idxxx,1),1)*(t+2020) idxxx]];
    emico2=S(t+1,17); % emissions GtCO2
    rff=(S(t+1,12)-sum(S(t+1,44:49),2))/(S(1,12)-sum(S(1,44:49),2)); % ratio of fossil fuel change
    S(t+1,(Ne*12+1):(Ne*12+46))=climdyn(t+1, S(t,(Ne*12+1):(Ne*12+46)), emico2, rff );
end
% idxxx_2020to2200
% 1.year; 2.reg; 3.type; 4.id
save('H:\Transfer of carbon-tax revenue\Ans\idxxx_2020to2200_CUR.mat','idxxx_2020to2200');
Power_all_CUR = S(:,Ne*([1:12]-1)+12)./3600*1000; % energy TWh/y
Power_renewable_CUR = S(:,Ne*([1:12]-1)+21)./3600*1000; % reneable energy TWh/y
save('H:\Transfer of carbon-tax revenue\Ans\Power_all_CUR.mat','Power_all_CUR');
save('H:\Transfer of carbon-tax revenue\Ans\Power_renewable_CUR.mat','Power_renewable_CUR');

