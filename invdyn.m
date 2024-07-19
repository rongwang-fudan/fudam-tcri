% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2024.1.28

function S2  =  invdyn( t1, S1 )
%   tyear: year
%   S1/S2: input and output
%   1 EUE; 2 EPE; 3 ENE; 4 saving rate;
%   5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
%   8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
%   10 capital (trill $); 11 energy capital (trill $); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega; 16 discounted utility
%   17 CO2 emissions Gt CO2; 18 carbon tax revenue t$; 19 climate investment t$; 20 new-built fossil fuel PJ/yr; 21 renewable energy PJ/yr; 22 remaining climate investment t$
%   23-29 abated emission GtCO2/yr; 30-36 new-built investments t$; 37-43 new-built renewable PJ/yr; 44-50 in-use renewable PJ/yr
%   51-61 investment transfer from region 1-11 t$; 62 OM costs t$/yr; 63 net transfer of money from other regions i t$/y

global Ne punits lifetimelcp tinv rlearning econo0 Ty0

S2=S1;
tyear = t1+Ty0-1;

% initializing the time of building all power units at the begining of the simulation
if t1==2
    punits(:,2)=punits(:,1);
end

% effects of technological improvements
punitcosts=zeros(size(punits,1),2);
for j=1:7
    elearning=ones(12,1);
    ene0=max(1000,econo0(1,43+j));
%     for i=2:12
%         elearning(i)=((S1(1,36+j+(i-1)*Ne)+ene0)/ene0)^(log2(1-rlearning(j)))-1; % fraction of low carbon energy cost reduction in a year to that in 2015
%     end
    ene1=max(1000,sum(S1(1,(36+j+Ne):Ne:(36+j+11*Ne)),2)); % new built power
    elearning(1)=((ene1+ene0)/ene0)^(log2(1-rlearning(j)))-1;
    idx=find(punits(:,3)==j);
    if size(idx,1)>0
        punitcosts(idx,1)=punits(idx,6).*(1+punits(idx,8).*elearning(1)); % investment t$
        punitcosts(idx,2)=punitcosts(idx,1)/lifetimelcp(j); % O&M costs t$
    end
end

% decommission of power units
for i=2:12
    extrapower=S2(1,21+(i-1)*Ne)-S2(1,12+(i-1)*Ne);
    if extrapower>0
        idx=find(punits(:,3)<7 & punits(:,4)==i & punits(:,2)<2300);
        if size(idx,1)>0
            k=size(idx,1);
            while extrapower>0 && k>0
                punits(idx(k),2)=2400; % decommission
                extrapower=extrapower-punits(idx(k),5); % reduce power PJ/yr 
                k=k-1;
            end
        end
    end
end

% build new power units in a region using climate investment from this region
pwdtotal=0; % total low-carbon power demand PJ/yr
invtotal=0; % total remaining investment t$
for i=2:12
    if S2(1,20+(i-1)*Ne)>1 && S2(1,22+(i-1)*Ne)>1e-4
        idx=find(punits(:,3)<7 & punits(:,4)==i & punits(:,2)>2300);
        if size(idx,1)>0
            for j=1:size(idx,1)
                j2=idx(j); % id of this power unit
                e=punits(j2,3); % type of energy - 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal
                if S2(1,20+(i-1)*Ne)<1
                    break; % no power demand
                end
                if S2(1,22+(i-1)*Ne)<1e-4
                    break; % no money
                end
                if S2(1,20+(i-1)*Ne)>=punits(j2,5) && S2(1,22+(i-1)*Ne)>=punitcosts(j2,1)
                    punits(j2,2)=tyear; % year of building this power unit
                    S2(1,20+(i-1)*Ne)=S2(1,20+(i-1)*Ne)-punits(j2,5); % reduce power demand PJ/yr 
                    S2(1,22+(i-1)*Ne)=S2(1,22+(i-1)*Ne)-punitcosts(j2,1); % reduce climate investment t$
                    S2(1,29+e+(i-1)*Ne)=S2(1,29+e+(i-1)*Ne)+punitcosts(j2,1); % investment t$ in 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
                    S2(1,36+e+(i-1)*Ne)=S2(1,36+e+(i-1)*Ne)+punits(j2,5); % new-built power generation PJ/yr by 1 solar, 2 wind, 3 biomass, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
                    S2(1,50+i-1+(i-1)*Ne)=S2(1,50+i-1+(i-1)*Ne)+punitcosts(j2,1); % transfer of investment from region i to region i
                end
            end
        end
    end
    invtotal=invtotal+S2(1,22+(i-1)*Ne); % total remaining investment t$
    pwdtotal=pwdtotal+S2(1,20+(i-1)*Ne); % total low-carbon power demand PJ/yr
end

% transferring the remaining carbon tax revenue across regions
if invtotal>0 && pwdtotal>0 && tinv>0
    i2=1;
    idx=find(punits(:,3)<7 & punits(:,2)>2300);
    for j=1:size(idx,1)
        j2=idx(j); % id of this power plant
        i=punits(j2,4); % region id of this power plant
        e=punits(j2,3); % type of energy - 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
        if pwdtotal<1
            break; % no power demand
        end
        if invtotal<1e-4
            break; % no money
        end
        if invtotal>=punitcosts(j2,1) && pwdtotal>=punits(j2,5) && S2(1,20+(i-1)*Ne)>=punits(j2,5)
            % transfer the money from region i2 to region i
            if S2(1,22+i2*Ne)>punitcosts(j2,1)
                S2(1,22+i2*Ne)=S2(1,22+i2*Ne)-punitcosts(j2,1); % reduce climate investment t$
                S2(1,50+i2+(i-1)*Ne)=S2(1,50+i2+(i-1)*Ne)+punitcosts(j2,1); % transfer of investment from region i2 to region i
                i2=i2+1;
                if i2>11
                    i2=1;
                end
            else
                % use the money from the next region
                for i3=1:12
                    i2=i2+1;
                    if i2>11
                        i2=1;
                    end
                    if S2(1,22+i2*Ne)>punitcosts(j2,1)
                        break;
                    end
                end
                if i3==12
                    continue;
                end
                S2(1,22+i2*Ne)=max(0,S2(1,22+i2*Ne)-punitcosts(j2,1)); % reduce climate investment t$
                S2(1,50+i2+(i-1)*Ne)=S2(1,50+i2+(i-1)*Ne)+punitcosts(j2,1); % transfer of investment from region i2 to region i
                i2=i2+1;
                if i2>11
                    i2=1;
                end
            end
            punits(j2,2)=tyear; % year of building this power unit
            pwdtotal=pwdtotal-punits(j2,5); % reduce power demand PJ/yr
            invtotal=invtotal-punitcosts(j2,1); % reduce climate investment t$
            S2(1,20+(i-1)*Ne)=S2(1,20+(i-1)*Ne)-punits(j2,5); % reduce power demand PJ/yr
            S2(1,29+e+(i-1)*Ne)=S2(1,29+e+(i-1)*Ne)+punitcosts(j2,1); % investment t$ in 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
            S2(1,36+e+(i-1)*Ne)=S2(1,36+e+(i-1)*Ne)+punits(j2,5); % new-built power generation PJ/yr by 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
        end
    end
end

% building CCS in a region after replacing 99% fossil fuel
for i=2:12
    if S2(1,22+(i-1)*Ne)>1e-4 && (S2(1,21+(i-1)*Ne)/S2(1,12+(i-1)*Ne))>0.95
        idx=find(punits(:,3)==7 & punits(:,4)==i & punits(:,2)>2300); % CCS
        if size(idx,1)>0
            for j=1:size(idx,1)
                j2=idx(j); % id of this power unit
                if S2(1,22+(i-1)*Ne)<1e-4
                    break; % no money
                end
                if S2(1,22+(i-1)*Ne)>=punitcosts(j2,1)
                    punits(j2,2)=tyear; % year of building this power unit
                    S2(1,22+(i-1)*Ne)=S2(1,22+(i-1)*Ne)-punitcosts(j2,1); % reduce climate investment t$
                    S2(1,29+7+(i-1)*Ne)=S2(1,29+7+(i-1)*Ne)+punitcosts(j2,1); % CCS investment t$
                    S2(1,36+7+(i-1)*Ne)=S2(1,36+7+(i-1)*Ne)+punits(j2,5); % new-built CCS power PJ/yr
                    S2(1,50+i-1+(i-1)*Ne)=S2(1,50+i-1+(i-1)*Ne)+punitcosts(j2,1); % transfer of investment from region i to region i
                end
            end
        end
    end
end

% abated emission GtCO2/yr for 1 biomass, 2 solar, 3 wind, 4 nuclear, 5 hydropower, 6 geothermal, 7 CCS
for i=2:12
    for e=1:7
        idx=find(punits(:,3)==e & punits(:,4)==i & punits(:,2)<=2300);
        if size(idx,1)>0
            S2(1,22+e+(i-1)*Ne)=sum(punits(idx,7),1); % abated emission GtCO2/yr
            S2(1,43+e+(i-1)*Ne)=sum(punits(idx,5),1); % low-carbon power PJ/yr
            S2(1,62+(i-1)*Ne)=S2(1,62+(i-1)*Ne)+sum(punitcosts(idx,2),1); % O&M cost t$/y
        end
    end
    for j=2:12
        if i~=j
%             S2(1,63+(i-1)*Ne)=S2(1,63+(i-1)*Ne)+S2(1,50+j-1+(i-1)*Ne)-S2(1,50+i-1+(j-1)*Ne); % net transfer of money from region j to i t$/y
            S2(1,63+(i-1)*Ne)=S2(1,63+(i-1)*Ne)+S2(1,50+j-1+(i-1)*Ne); % net transfer of money from region j to i t$/y
        end
    end
end

% global total or average
S2(1,1:61)=0;
for i=2:12
    for j=1:63 %62
        if j<=4 || j==6 || j==8 || j==9 || j==14 || j==15
            S2(1,j)=S2(1,j)+S2(1,j+(i-1)*Ne)/11; % average
        else
            S2(1,j)=S2(1,j)+S2(1,j+(i-1)*Ne); % sum
        end
    end
end

end



