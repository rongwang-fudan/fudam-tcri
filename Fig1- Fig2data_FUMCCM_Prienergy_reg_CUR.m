tic
clear;
clear global;
Initialset_econ;
Initialset_clim;
fracinv0=[0 300 52 45 39 19 325 175 10 14 10 10]; % Figure 20 in Global Landscape of Climate Finance 2023, Climate Policy Initiative, https://www.climatepolicyinitiative.org/publication/global-landscape-of-climate-finance-2023/
fracinv0(1,1)=sum(fracinv0(1,2:12),2); fracinv0=fracinv0./fracinv0(1,1);
inv15=[0.340 0.263 0.351 0.322 0.329 0.348 0.430 0.499; 75.28 76.52 81.48 86.54 87.78 85.26 97.53 101.33; 0.00452 0.00344 0.00431 0.00372 0.00375 0.00408 0.00441 0.00492]; % 1 renewable energy investment t$; 2 gdp; 3 ratio from 2015 to 2022
t15=[2015:2022]; [sR,lrv0,bb0] = regression(t15,inv15(3,:));
rlearning0=rlearning;
lifetimeff0=lifetimeff;
econo00=econo0;

Ty2=2200;
S2=zeros(Ty2-Ty0+1,Ne*12+46,8);
for sce=1:8
    % set-up of experiments
    fcpc=1; % fraction of carbon tax revenue that is used for renewable energy investments
    tinv=1; % transfer of carbon tax revenue between regions: 0 no; 1 yes
    lrv=lrv0; % rate of growth in renewable energy investment
    cpc0=100; % carbon price $/tCO2
    rlearning=rlearning0;
    lifetimeff=lifetimeff0;
    econo0=econo00;
    if sce==1
        cpc0=0;
    elseif sce==3
        tinv=0;
    elseif sce==4
        rlearning=rlearning0*0;
    elseif sce==5
        lifetimeff=lifetimeff0*2;
        econo0(1,(Ne+20):Ne:(Ne*11+20))=econo0(1,(Ne+20):Ne:(Ne*11+20))/2;
    elseif sce==6
        fcpc=0;
    elseif sce==7
        cpc0=0;
        lrv=lrv0*2;
    elseif sce==8
        cpc0=100;
        lrv=lrv0*2;
    end
    % simulations
    S=zeros(Ty2-Ty0+1,Ne*12+46); % ne economic variables for 12 regions + 46 climatic variables for the world
    S(1,1:(Ne*12))=econo0(1,1:(Ne*12)); % economic variable in 2020
    S(1,(Ne*12+1):(Ne*12+11))=clim0(1,1:11); % climatic variable in 2020
    for t=1:(Ty2-Ty0) % 2020 to 2303
        tyear=t+Ty0-1; display(tyear);
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
        S(t+1,1:(Ne*12))=invdyn(t+1,S(t+1,1:(Ne*12)));
        emico2=S(t+1,17); % emissions GtCO2
        rff=(S(t+1,12)-sum(S(t+1,44:49),2))/(S(1,12)-sum(S(1,44:49),2)); % ratio of fossil fuel change
        S(t+1,(Ne*12+1):(Ne*12+46))=climdyn(t+1, S(t,(Ne*12+1):(Ne*12+46)), emico2, rff );
    end
    S2(:,:,sce)=S;
end
save('H:\Transfer of carbon-tax revenue\Ans\fig2_sce5_2200.dat','S2');


%%
tic
clear;
clear global;
Initialset_econ;
Initialset_clim;

load('H:\Transfer of carbon-tax revenue\Ans\fig2_sce5_2200.dat','-mat');
Ty2=2200;
Ne = 63;
ns=7;
output=zeros(Ty2-Ty0+1,4*5*ns+ns);
for sce=1:ns
    S=S2(:,:,sce);
    output(:,sce)=S(:,8+Ne*12); % Warming C
    output(:,ns+sce)=S(:,Ne*(1-1)+5)+S(:,Ne*(1-1)+63); % renewable energy investments t$
    output(:,ns*2+sce)=S(:,Ne*(2-1)+5)+S(:,Ne*(2-1)+63); % renewable energy investments t$
    output(:,ns*3+sce)=S(:,Ne*(7-1)+5)+S(:,Ne*(8-1)+5)+S(:,Ne*(10-1)+5)+S(:,Ne*(7-1)+63)+S(:,Ne*(8-1)+63)+S(:,Ne*(10-1)+63); % renewable energy investments t$
    output(:,ns*4+sce)=S(:,Ne*(4-1)+5)+S(:,Ne*(5-1)+5)+S(:,Ne*(6-1)+5)+S(:,Ne*(4-1)+63)+S(:,Ne*(5-1)+63)+S(:,Ne*(6-1)+63); % renewable energy investments t$
    output(:,ns*5+sce)=S(:,Ne*(1-1)+17); % emission GtCO2/yr, global
    output(:,ns*6+sce)=S(:,Ne*(2-1)+17); % emission GtCO2/yr, East Asia
    output(:,ns*7+sce)=S(:,Ne*(7-1)+17)+S(:,Ne*(8-1)+17)+S(:,Ne*(10-1)+17); % emission GtCO2/yr, North America+Europe+Pacific OECD
    output(:,ns*8+sce)=S(:,Ne*(4-1)+17)+S(:,Ne*(5-1)+17)+S(:,Ne*(6-1)+17); % emission GtCO2/yr, South Asia+Africa+Middle East
    output(:,ns*9+sce)=S(:,Ne*(1-1)+12)./3600; % energy PWh
    output(:,ns*10+sce)=S(:,Ne*(2-1)+12)./3600; % energy PWh
    output(:,ns*11+sce)=(S(:,Ne*(7-1)+12)+S(:,Ne*(8-1)+12)+S(:,Ne*(10-1)+12))./3600; % energy PWh
    output(:,ns*12+sce)=(S(:,Ne*(4-1)+12)+S(:,Ne*(5-1)+12)+S(:,Ne*(6-1)+12))./3600; % energy PWh
    output(:,ns*13+sce)=S(:,Ne*(1-1)+21)./3600; % renewable energy PWh
    output(:,ns*14+sce)=S(:,Ne*(2-1)+21)./3600; % renewable energy PWh
    output(:,ns*15+sce)=(S(:,Ne*(7-1)+21)+S(:,Ne*(8-1)+21)+S(:,Ne*(10-1)+21))./3600; % renewable energy PWh
    output(:,ns*16+sce)=(S(:,Ne*(4-1)+21)+S(:,Ne*(5-1)+21)+S(:,Ne*(6-1)+21))./3600; % renewable energy PWh
    output(:,ns*17+sce)=S(:,Ne*(1-1)+16)-S2(:,Ne*(1-1)+16,1); % Utility change relative to baseline
    output(:,ns*18+sce)=S(:,Ne*(2-1)+16)-S2(:,Ne*(2-1)+16,1); % Utility change relative to baseline
    output(:,ns*19+sce)=S(:,Ne*(7-1)+16)+S(:,Ne*(8-1)+16)+S(:,Ne*(10-1)+16)-S2(:,Ne*(7-1)+16,1)-S2(:,Ne*(8-1)+16,1)-S2(:,Ne*(10-1)+16,1); % Utility change relative to baseline
    output(:,ns*20+sce)=S(:,Ne*(4-1)+16)+S(:,Ne*(5-1)+16)+S(:,Ne*(6-1)+16)-S2(:,Ne*(4-1)+16,1)-S2(:,Ne*(5-1)+16,1)-S2(:,Ne*(6-1)+16,1); % Utility change relative to baseline
end

Prienergy_reg=zeros(Ty2-Ty0+1,12);
for sce=1
    S=S2(:,:,sce);
    for i = 1:12
        Prienergy_reg(:,i)=S(:,Ne*(i-1)+12)./3600; % energy PWh
    end
end
save('H:\Transfer of carbon-tax revenue\Ans\Fig2data_FUMCCM_Prienergy_reg_CUR.mat','Prienergy_reg')