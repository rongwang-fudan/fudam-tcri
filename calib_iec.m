% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function iec = calib_iec ( inputs, alpha, elas )
%   iec:  regresssion coefficients
%   alpha:  elasticity of output to capital
%   elas:   elasticity of substitution between energy and non-energy in output
%   inputs 12,8*48: 1971-2018 and 1 Gt CO2/yr; 2 gdp t$; 3 population M; 4 capital t$; 5 energy price ($/kWh); 6 energy PWh; 7 renewable PWh; 8 nuclear PWh

plots = 0;
iec=zeros(3*12,9); % 1:2 linear regression coefficient; 3-6 uncertainty in the linear regression coefficient; 7-8 min and max of x

for z=1:12

tlen=26;
xy_iec=zeros(tlen,8);
for i=1:tlen
    x=[(1970+i):(1989+i)]; % time
    y=inputs(z,(i+192):(i+211)); [sR,lr_pe0,bb0] = regression(x,log(y)); % energy price
    y=inputs(z,(i+240):(i+259))./inputs(z,(i+96):(i+115)); [sR,lr_e0,bb0] = regression(x,log(y)); % e=E/L
    y=inputs(z,(i+48):(i+67))./inputs(z,(i+96):(i+115)); [sR,lr_y0,bb0] = regression(x,log(y)); % y=Y/L
    y=inputs(z,(i+144):(i+163))./inputs(z,(i+96):(i+115)); [sR,lr_k0,bb0] = regression(x,log(y)); % k=K/L
    y=(inputs(z,(i+192):(i+211)).*inputs(z,(i+240):(i+259)))./inputs(z,(i+48):(i+67)); avef0=mean(y,1); startf0=y(1); endf0=y(end); % omega=pe*E/Y
    lr_se0 = lr_pe0 + lr_e0 - lr_y0;
    lr_B = lr_e0 - alpha*lr_k0;
    lr_A = lr_y0 - alpha*lr_k0;
    lr_taue0 = elas/(elas-1)*lr_se0 - lr_e0 + lr_y0;
    lr_we0 = lr_B - lr_se0;
    lr_b0 = (lr_A - avef0*(lr_we0 + lr_taue0))/(1-avef0);
    xy_iec(i,1) = startf0;
    xy_iec(i,2) = lr_taue0;
    xy_iec(i,3) = lr_we0;
    xy_iec(i,4) = lr_b0;
    xy_iec(i,8) = endf0;
end

if z==6
    xy_iec=xy_iec(10:end,:); % incomplete data for middle east before 1989
end
if z==11
    aa=xy_iec;
    xy_iec=aa(10:26,:); 
    xy_iec(1:8,:)=aa(1:8,:); % remove data for Russia from 1988 to 1996
end

% EUE
x=log10(xy_iec(:,1)); y=xy_iec(:,2); tlen=size(xy_iec,1);
[b22,bint22,r22,rint22,statseue] = regress(y,[ones(tlen,1) x]);
iec(1+(z-1)*3,1:2)=b22(1:2,1); iec(1+(z-1)*3,3:4)=bint22(1:2,1); iec(1+(z-1)*3,5:6)=bint22(1:2,2); iec(1+(z-1)*3,7)=min(x,[],1); iec(1+(z-1)*3,8)=max(x,[],1);
[b2eue, Seue]=polyfit(x,y,1);
xeue=[min(x,[],1):0.001:max(x,[],1)]; [yfiteue, deltaeue]=polyconf(b2eue,xeue,Seue,'alpha',0.05,'predopt','curve');
xy_iec(:,5)=polyval(b2eue,x);

% EPE
y=xy_iec(:,3);
[b22,bint22,r22,rint22,statsepe] = regress(y,[ones(tlen,1) x]);
iec(2+(z-1)*3,1:2)=b22(1:2,1); iec(2+(z-1)*3,3:4)=bint22(1:2,1); iec(2+(z-1)*3,5:6)=bint22(1:2,2); iec(2+(z-1)*3,7)=min(x,[],1); iec(2+(z-1)*3,8)=max(x,[],1);
[b2epe, Sepe]=polyfit(x,y,1);
xepe=[min(x,[],1):0.001:max(x,[],1)]; [yfitepe, deltaepe]=polyconf(b2epe,xepe,Sepe,'alpha',0.05,'predopt','curve');
xy_iec(:,6)=polyval(b2epe,x);

% ENE
x=log10(1-xy_iec(:,1)); y=xy_iec(:,4);
[b22,bint22,r22,rint22,statsene] = regress(y,[ones(tlen,1) x]);
iec(3+(z-1)*3,1:2)=b22(1:2,1); iec(3+(z-1)*3,3:4)=bint22(1:2,1); iec(3+(z-1)*3,5:6)=bint22(1:2,2); iec(3+(z-1)*3,7)=min(x,[],1); iec(3+(z-1)*3,8)=max(x,[],1);
[b2ene, Sene]=polyfit(x,y,1);
xene=[min(x,[],1):0.001:max(x,[],1)]; [yfitene, deltaene]=polyconf(b2ene,xene,Sene,'alpha',0.05,'predopt','curve');
xy_iec(:,7)=polyval(b2ene,x);

if plots==1
    subplot(12,3,1+(z-1)*3);
    x=log10(xy_iec(:,1));
    y=xy_iec(:,2); % EUE rates    
    xconf = [xeue xeue(end:-1:1)] ; 
    yconf = [yfiteue+deltaeue yfiteue(end:-1:1)-deltaeue];
    p = fill(xconf,yconf,'red'); hold on;
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none'; 
    plot(xeue,yfiteue,'k-','LineWidth',3); hold on;
    plot(x,y,'o','MarkerEdgeColor',[0.8 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;
    title(['function: ',texlabel(strcat('function:',num2str(round(b2eue,4)),',r2=',num2str(round(statseue(1),4))))]);
    
    subplot(12,3,2+(z-1)*3);
    x=log10(xy_iec(:,1));
    y=xy_iec(:,3); % EPE rates
    xconf = [xepe xepe(end:-1:1)] ; 
    yconf = [yfitepe+deltaepe yfitepe(end:-1:1)-deltaepe];
    p = fill(xconf,yconf,'red'); hold on;
    p.FaceColor = [1 0.8 0.8];
    p.EdgeColor = 'none'; 
    plot(xepe,yfitepe,'k-','LineWidth',3); hold on;
    plot(x,y,'o','MarkerEdgeColor',[0.8 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;
    title(['function: ',texlabel(strcat('function:',num2str(round(b2epe,4)),',r2=',num2str(round(statsepe(1),4))))]);

    subplot(12,3,3+(z-1)*3);
    x=log10(1-xy_iec(:,1));
    y=xy_iec(:,4); % ENE rates
    xconf = [xene xene(end:-1:1)] ; 
    yconf = [yfitene+deltaene yfitene(end:-1:1)-deltaene];
    p = fill(xconf,yconf,'red'); hold on;
    p.FaceColor = [1 0.8 0.8];      
    p.EdgeColor = 'none'; 
    plot(xene,yfitene,'k-','LineWidth',3); hold on;
    plot(x,y,'o','MarkerEdgeColor',[0.8 0 0],'MarkerFaceColor','none','MarkerSize',3); hold on;
    title(['function: ',texlabel(strcat('function:',num2str(round(b2ene,4)),',r2=',num2str(round(statsene(1),4))))]);    
end

end

end



