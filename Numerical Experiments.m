%% Font and Color Setup

set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');

str = '#0065bd';
TUMblue = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#072140';
TUMblue_dark= sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#0E396E';
TUMblue_dark2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#14519A';
TUMblue_dark3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#5E94D4';
TUMblue_light = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#C2D7EF';
TUMblue_light2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#E3EEFA';
TUMblue_light3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#98c6ea';
lighterblue_TUM = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#8F81EA';
brightblue_TUM = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#e37222';
TUMorange = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#EA7237';
TUMred = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#999999';
gray_TUM = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#a2ad00';
Green_TUM = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#A2142F';
darkred = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;


%% Comparison of various Wirtinger gradient descent methods (Figures 2-5)

n = 256;
m = 10*n;

rng(1)

num =  20; 
num2 = 6; 

dose = zeros(1,num2);

dose1 = zeros(1,num2);
dose2 = zeros(1,num2);
dose3 = zeros(1,num2);
dose4 = zeros(1,num2);
dose5 = zeros(1,num2);
doseWF = zeros(1,num2);
doseTWF = zeros(1,num2);
doseoF = zeros(1,num2);
doseoF_without0 = zeros(1,num2);


for k = 1:num2

    average1 = 0;
    average2 = 0;
    average3 = 0;
    average4 = 0;
    average5 = 0;
    averageWF = 0;
    averageTWF = 0;
    averageoF = 0;
    averageoF_without0 = 0;

    for i = 1:num 

    rng(i)

    % Generate data
    x = (randn(n,1) + randn(n,1)*1.0i);
     
    A_mod = randn(m,n) + randn(m,n)*1.0i; 
    for j = 1:m
         A_mod(j,:) = A_mod(j,:)/norm(A_mod(j,:));
    end
    b_noisefree = abs(A_mod*x).^2;
    dose_model = sum(b_noisefree);
    A = A_mod/sqrt(dose_model);
    b_noisefree = abs(A*x).^2;    
    
    dosis = 10^4 * 0.05*k;
    
    A = A*sqrt(dosis);

% Generate noise
b = poissrnd(abs(A*x).^2);

SNR = norm(dosis*b_noisefree)/norm(b-dosis*b_noisefree);

if (k == 1 && i == num)

dose_b = sum(b);
figure; 
histogram(b,'FaceColor','#999999', 'LineWidth',1.5, 'EdgeColor','#999999')
xlabel('Measurement','FontSize',32); ylabel('Count','FontSize',32);
set(gca,'FontSize',32)

SNR = norm(dosis*b_noisefree)/norm(b-dosis*b_noisefree)

elseif (k == num2 && i == num)

dose_b = sum(b);
figure; 
histogram(b,'FaceColor','#999999', 'LineWidth',1.5, 'EdgeColor','#999999')
xlabel('Measurement','FontSize',32); ylabel('Count','FontSize',32);
set(gca,'FontSize',32)

SNR = norm(dosis*b_noisefree)/norm(b-dosis*b_noisefree)
end

dose(k) = sum(b);

% Initialization
y = b;
npower_iter = 10;                           
z0 = randn(n,1); z0 = z0/norm(z0,'fro');    
for tt = 1:npower_iter                      
    z0 = A'*(y.* (A*z0)); z0 = z0/norm(z0,'fro');
end
normest = sqrt(sum(y)/numel(y));           
x0 = normest * z0;

% Compare different Wirtinger gradient descent algorithms
T = 1000;

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,1,i,10^(-3));
average2 = average2 + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,1,i,10^(-1));
average1 = average1 + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,1,i,0.25);
average3 = average3 + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,1,i,0.5);
average4 = average4 + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,1,i,1);
average5 = average5 + rel_err(Tstop);
   
[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,0,i); 
averageWF = averageWF + rel_err(Tstop);
 
[~,~,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,2,i);
averageTWF = averageTWF + rel_err(end);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,4,i);
averageoF = averageoF + rel_err(Tstop); 
  
[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,3,i);
averageoF_without0 = averageoF_without0 + rel_err(Tstop);

end

dose2(k) = average2/num;
dose1(k) = average1/num;
dose3(k) = average3/num;
dose4(k) = average4/num;
dose5(k) = average5/num;
doseWF(k) = averageWF/num;
doseTWF(k) = averageTWF/num;
doseoF(k) = averageoF/num;
doseoF_without0(k) = averageoF_without0/num;

end

dose = 10^4 * 0.05*(1:num2);

figure; hold on
plot(dose,dose2,'x','Color',darkred,'LineWidth',2.5,'MarkerSize',16);
plot(dose,dose1,'^','Color',TUMblue,'LineWidth',2.5,'MarkerSize',16);
plot(dose,dose3,'square','Color',TUMorange,'LineWidth',2.5,'MarkerSize',14);
plot(dose,dose4,'o','Color',Green_TUM,'LineWidth',2.5,'MarkerSize',14);
plot(dose,dose5,'diamond','Color',brightblue_TUM,'LineWidth',2.5,'MarkerSize',14);

plot(dose,dose2,'k--','LineWidth',0.5);
plot(dose,dose1,'k--','LineWidth',0.5); 
plot(dose,dose3,'k--','LineWidth',0.5); 
plot(dose,dose4,'k--','LineWidth',0.5);
plot(dose,dose5,'k--','LineWidth',0.5);

legend({'~PF $\varepsilon_P = 10^{-3}$','~PF $\varepsilon_P = 0.1$','~PF $\varepsilon_P = 0.25$','~PF $\varepsilon_P = 0.5$','~PF $\varepsilon_P = 1$'},'FontSize',26);
xlabel('Dose','FontSize',36); ylabel('Relative reconstruction error','FontSize',36);
ylim([0.2 1.3]); 
xlim([dose(1)-50,dose(end)+50])
set(gca,'FontSize',32)

figure; hold on
plot(dose,doseoF,'-','Color',TUMred,'LineWidth',4);
plot(dose,doseoF_without0,'k--','LineWidth',3);

plot(dose,dose2,'x','Color',gray_TUM,'LineWidth',1.5,'MarkerSize',16);
plot(dose,dose1,'^','Color',gray_TUM,'LineWidth',1.5,'MarkerSize',16);
plot(dose,dose3,'square','Color',gray_TUM,'LineWidth',1.5,'MarkerSize',14);
plot(dose,dose4,'o','Color',gray_TUM,'LineWidth',1.5,'MarkerSize',14);
plot(dose,dose5,'diamond','Color',gray_TUM,'LineWidth',1.5,'MarkerSize',14);

plot(dose,dose2,'k--','LineWidth',0.2);
plot(dose,dose1,'k--','LineWidth',0.2); 
plot(dose,dose3,'k--','LineWidth',0.2); 
plot(dose,dose4,'k--','LineWidth',0.2);
plot(dose,dose5,'k--','LineWidth',0.2);

legend({'~FIVS',sprintf('~FIVS-WA'),'~PF $\varepsilon_P = 10^{-3}$','~PF $\varepsilon_P = 0.1$','~PF $\varepsilon_P = 0.25$','~PF $\varepsilon_P = 0.5$','~PF $\varepsilon_P = 1$'},'FontSize',26);
xlabel('Dose','FontSize',26); ylabel('Relative reconstruction error','FontSize',26);
ylim([0.2 1.3]); 
xlim([dose(1)-50,dose(end)+50])
set(gca,'FontSize',32)

figure; hold on
plot(dose,doseoF,'-','Color',TUMred,'LineWidth',4);
plot(dose,doseWF,'-','Color',lighterblue_TUM,'LineWidth',4);
plot(dose,doseTWF,'-','Color',Green_TUM,'LineWidth',4);

legend({'~FIVS','~WF','~TWF'},'FontSize',26);
xlabel('Dose','FontSize',26); ylabel('Relative reconstruction error','FontSize',26);
ylim([0.2 1.3]); 
xlim([dose(1)-50,dose(end)+50])
set(gca,'FontSize',32)



%% Oversampling - Dose - TradeOff (Figure 6)

rng(1)

num = 5; % number of trials
num1 = 15; % number of different dose parameters
num2 = 15; % number of different oversampling ratios

for run = 1:3

dose = zeros(1,num1*num2);
oversampling = zeros(1,num1*num2);

dose2 = zeros(1,num1*num2);
dose1 = zeros(1,num1*num2);
dose3 = zeros(1,num1*num2);
dose4 = zeros(1,num1*num2);
dose5 = zeros(1,num1*num2);
doseAF = zeros(1,num1*num2);
doseoF = zeros(1,num1*num2);


for l = 1:num1
for k = 1:num2

    average2 = 0;
    average1 = 0;
    average3 = 0;
    average4 = 0;
    average5 = 0;  
    averageAF = 0;
    averageoF = 0;
    

    for i = 1:num 
    rng(i) 

    n = 2^(run + 4);
    fun_m = @(x) (4+x-1);
    m = fun_m(k)*n; 
    
    % Generate data
    x = (randn(n,1) + randn(n,1)*1.0i);
    
    A_mod = randn(m,n) + randn(m,n)*1.0i;
    b_noisefree = abs(A_mod*x).^2;
    dose_model = sum(b_noisefree);
    A = A_mod/sqrt(dose_model);
    b_noisefree = abs(A*x).^2;
    
    dosis =  l*n;
    A = A*sqrt(dosis);

% Generate noise
b = poissrnd(abs(A*x).^2);

dose((l-1) * num2 + k) = dosis;
oversampling((l-1) * num2 + k) = m/n;


% Initialization 
y = b;
npower_iter = 10;                           
z0 = randn(n,1); z0 = z0/norm(z0,'fro');   
for tt = 1:npower_iter                      
    z0 = A'*(y.* (A*z0)); z0 = z0/norm(z0,'fro');
end
normest = sqrt(sum(y)/numel(y));            
x0 = normest * z0;

T = 1000;

[~,Tstop,rel_err_AF,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,10^(-3));
average2 = average2 + rel_err_AF(Tstop);

[~,Tstop,rel_err_AF,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,10^(-1));
average1 = average1 + rel_err_AF(Tstop);

[~,Tstop,rel_err_AF,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,0.25);
average3 = average3 + rel_err_AF(Tstop);

[~,Tstop,rel_err_AF,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,0.5);
average4 = average4 + rel_err_AF(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,1);
average5 = average5 + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,5,i,10^(-6));
averageAF = averageAF + rel_err(Tstop);

[~,Tstop,rel_err,~] = WirtingerFlowAlgorithms(A,b,x0,T,x,4,i);
averageoF = averageoF + rel_err(Tstop);

end

dose1((l-1) * num2 + k) = average1/num;
dose2((l-1) * num2 + k) = average2/num;
dose3((l-1) * num2 + k) = average3/num;
dose4((l-1) * num2 + k) = average4/num;
dose5((l-1) * num2 + k) = average5/num;
doseAF((l-1) * num2 + k) = averageAF/num;
doseoF((l-1) * num2 + k) = averageoF/num;

end
end


dose = dose/n;

figure(2); hold on;
subplot(1, 3, run);

h1 = plot3(dose,oversampling,doseAF,'square','Color',TUMblue_dark,'MarkerFaceColor',TUMblue_dark,'LineWidth',1,'MarkerSize',50); hold on;
h2 = plot3(dose,oversampling,dose2,'square','Color',TUMblue_dark2,'MarkerFaceColor',TUMblue_dark2,'LineWidth',1,'MarkerSize',50); 
h3 = plot3(dose,oversampling,dose1,'square','Color',TUMblue_dark3,'MarkerFaceColor',TUMblue_dark3,'LineWidth',1,'MarkerSize',50);
h4 = plot3(dose,oversampling,dose3,'square','Color',TUMblue_light,'MarkerFaceColor',TUMblue_light,'LineWidth',1,'MarkerSize',50);
h5 = plot3(dose,oversampling,dose4,'square','Color',TUMblue_light2,'MarkerFaceColor',TUMblue_light2,'LineWidth',1,'MarkerSize',50);
h6 = plot3(dose,oversampling,dose5,'square','Color',TUMblue_light3,'MarkerFaceColor',TUMblue_light3,'LineWidth',1,'MarkerSize',50);
h7 = plot3(dose,oversampling,doseoF,'square','Color',TUMorange,'MarkerFaceColor',TUMorange,'LineWidth',1,'MarkerSize',50);

xlabel('Dose per object pixel $(\frac{d}{n})$')
ylabel('Oversampling ratio $(\frac{m}{n})$')
zlabel('rel.rec.err.')
title(['$n~=~$', num2str(n)]);

set(gca,'FontSize',32)
xlim([dose(1)-0.6,dose(end)+0.6])
ylim([oversampling(1)-0.6,oversampling(end)+0.6])

view(0,-90)
axis vis3d

pos = get(gca, 'Position');  
pos(3) = 1.1 * pos(3);       
set(gca, 'Position', pos);   

allHandles = [h1,h2,h3,h4,h5,h6,h7];
legend(allHandles, {'AF $\varepsilon_A = 10^{-6}$','AF $\varepsilon_A = 10^{-3}$','AF $\varepsilon_A = 0.1$','AF $\varepsilon_A = 0.25$','AF $\varepsilon_A = 0.5$','AF $\varepsilon_A = 1$','FIVS'}, ...
    'NumColumns',8,'FontSize',32, 'Position', [0.46, 0.1, 0.1, 0.1], 'Orientation', 'horizontal');

end

