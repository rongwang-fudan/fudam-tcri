% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2024.1.30

function econo2  =  econdyn( t1, z, econo1, finv, cpc, fcpc, ffout, deltaT, dam )
%   Finds the next state (t1+1) given the current state and actions at time t1
%   cpc: carbon price $/tCO2
%   finv:  fraction of climate investments in GDP
%   ffout:  fossil-fuel retired PJ/yr
%   deltaT:  atm temperature
%   dam:  damage by tipping

global prtp bhm efco20 econo0 alpha elas elasmu pop
global cesin dk_x dk_e iec T miec fsav falab fainv Ne Ty0
%   econo:  economic variables over time
%   iec:  rates of induced efficiency changes
%   1 EUE; 2 EPE; 3 ENE; 4 Saving rate;
%   5 Abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
%   8 Fraction of labor allocation to energy; 9 fraction of investment allocation to energy
%   10 Capital (trill $); 11 energy capital (trill $); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega; 16 discounted utility
%   17 CO2 emissions Gt CO2; 18 carbon tax revenue t$; 19 climate investment t$; 20 new-built fossil fuel PJ/yr; 21 renewable energy PJ/yr; 22 remaining climate investment t$
%   23-29 Abated emission GtCO2/yr; 30-36 new-built investments t$; 37-43 new-built renewable PJ/yr; 44-50 in-use renewable PJ/yr
%   51-61 Investment transfer from region 1-11 t$; 62 OM costs t$/yr; 63 net transfer of money from other regions i t$/y
%   efco20: CO2 emission factors for fossil fuel only tCO2 / MJ
%   alpha: Elasticity of output to capital
%   elas: Elasticity of substitution between energy and non-energy in output

%Time step (year)
tstep = 1;
tyear = t1+Ty0-1;
tpop  = pop(z,min(T,tyear-1971));
tpop2  = pop(z,min(T,tyear-1970));
ziec = iec((1+(z-1)*3):(3+(z-1)*3),:);

%Economic variables
econo2 = zeros(1,Ne);

%Impact of warming on the productivity: bhm 1-2 for ax+bx2; bhm 3 for preindustrial temperature 
dpwarm = (bhm(1)*deltaT+bhm(2)*((deltaT+bhm(3))^2-bhm(3)^2));

%EUE rate with the inducing effects
euerate = ziec(1,1)+ziec(1,9)*min(1,(t1-2)/75)+ziec(1,2)*log10(econo1(15));

%EPE rate with the inducing effects
eperate = ziec(2,1)+ziec(2,9)*min(1,(t1-2)/75)+ziec(2,2)*log10(econo1(15));

%ENE rate with the inducing effects
enerate = ziec(3,1)+ziec(3,9)*min(1,(t1-2)/75)+ziec(3,2)*log10(1-econo1(15));

%no feedback of mitigation to efficiencies
if miec==0
    euerate = ziec(1,1)+ziec(1,2)*log10(econo0(1,15+(z-1)*Ne));
    eperate = ziec(2,1)+ziec(2,2)*log10(econo0(1,15+(z-1)*Ne));
    enerate = ziec(3,1)+ziec(3,2)*log10(1-econo0(1,15+(z-1)*Ne));
end

%EUE $ / KJ
econo2(1) = econo1(1)*(1+dpwarm+euerate);

%EPE PJ / (trillion $)^0.3 / (billion cap)^0.7
econo2(2) = econo1(2)*(1+dpwarm+eperate);

%ENE (trillion $)^0.7 / (billion people)^0.7
econo2(3) = econo1(3)*(1+dpwarm+enerate);

%Abatement costs t$
econo2(5) = min(econo1(13)*0.5, max(0, econo1(62) + sum(econo1(1,30:36),2) - econo1(63)));

%Damage by climate change as a percentage to output
econo2(6) = dam;

%Economic net output minus damage t$
econo2(7) = econo1(13) * (1-econo2(6)) - econo2(5);

%Investment t$
I = econo2(7) * econo1(4);

%Equivalent energy price $/kJ
ecpc = cpc * max(0, econo1(17)/econo1(12)/1000);

%Parameters for labor/investment allocation: 1ke, 2kx, 3eue, 4epe, 5ene, 6investment, 7labor, 8timestep,9carbonprice
cesin = [ econo1(11),econo1(10)-econo1(11),econo2(1),econo2(2),econo2(3),I,tpop/1000,tstep,ecpc ];

%Optimal allocation of labor and investment
allo = ces_allocation(1/(1+(econo2(3)/econo2(1)/econo2(2))^(elas-1)), 1); % 0 for percentile decimal, 1 for thousand decimal

%Feedback of mitigation to labor reallocation
if falab==0
    allo(1)=econo0(1,8+(z-1)*Ne);
end

%Feedback of mitigation to investment reallocation
if fainv==0
    allo(2)=econo0(1,9+(z-1)*Ne);
end

%Allocation of labor to produce energy
econo2(8) = allo(1);

%Allocation of investment to produce energy
econo2(9) = allo(2);

%Capital trill $
econo2(10) = tstep * I + (1 - dk_x) ^ tstep * (econo1(10)-econo1(11)) + (1 - dk_e) ^ tstep * econo1(11);

%Energy capital trill $
econo2(11) = tstep * I * econo2(9) + (1 - dk_e) ^ tstep * econo1(11);

%Energy PJ/yr
econo2(12) = econo2(2) * econo2(11)^alpha * ( tpop * econo2(8)/1000)^(1-alpha);

%Non-Energy labor fraction
if econo1(8)>econo2(8)
    nlabor=1-econo1(8)+(econo1(8)-econo2(8))/20; % 20 years to transfer labor from energy to non-energy sectors
else
    nlabor=1-econo2(8);
end

%Non-Energy production
X = econo2(3) * (econo2(10) - econo2(11))^alpha * ( tpop * nlabor/1000)^(1-alpha);

%Output trill$/yr
econo2(13) = ((econo2(1)*econo2(12))^((elas-1)/elas) + X^((elas-1)/elas))^(1/((elas-1)/elas));

%Energy price $/kWh final energy
econo2(14) = (1/((econo2(3)/econo2(1)/econo2(2))^(elas-1)+1) * econo2(13) / econo2(12) + ecpc) * 3600;

%Share of energy expenditure in GDP
econo2(15) = econo2(12) * econo2(14) / econo2(13) / 3600;

%population-weighted utility of per capita consumption
utility = ((econo2(7)*(1-econo1(4))/tpop*1000)^(1-elasmu))/(1-elasmu) * tpop/1000/(1+prtp)^max(0,(tyear-Ty0));

%utility denominated in terms of current consumption t$
econo2(16) = utility / (econo0(1,16+(z-1)*Ne)*(1-elasmu)) * econo0(1,7+(z-1)*Ne)*(1-econo0(1,4+(z-1)*Ne));

%industrial emission GtCO2/yr
econo2(17) = efco20(z) * max(0,(econo2(12)-sum(econo1(1,44:49),2)))-econo1(1,29);

%carbon tax revenue t$
econo2(18) = max(0, econo2(17)) * cpc / 1000 * fcpc;

%climate investment t$
econo2(19) = econo2(7) * finv;

%remaining fossil fuel PJ/yr
ff = max(0, econo2(12)-sum(econo1(1,44:49),2));

%power demand PJ/yr
econo2(20) = max(0, econo2(12)-econo1(12)+min(ffout, ff));

%renewable energy PJ/yr
econo2(21) = sum(econo1(1,44:49),2);

%renewable energy investment t$
econo2(22) = max(0, econo1(22) * econo1(4) + econo2(19) + econo2(18) - econo1(62));

%reduce mitigation when global warming or emission is low
if deltaT<1 || econo2(17)<-1
    econo2(22)=econo2(22)/2;
end

%cumulative low-carbon power PJ/yr
econo2(1,37:50) = econo1(1,37:50);

%growth of per capita output
dlny_dt=(econo2(13)/tpop2)/(econo1(13)/tpop)-1;

%change rate of investment based on the rate of per capita consumption: dlnc/dt
dln1s_dt = (alpha*econo1(13)/econo1(10)-dk_x-prtp)/elasmu * 0.8253 - 0.0097 - dlny_dt;

%saving rate
if fsav==1
    econo2(4) = max(0.1,min(0.3,1-(1-econo1(4))*(1-dln1s_dt)));
else
    econo2(4) = econo1(4);
end

end



