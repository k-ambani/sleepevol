%% Creating data structure
%This section should have all the matrices that are required for other
%analysis code files. 
data_struct=1;
tic

paramsets = readtable('0_paramsetrec.xlsx');
paramsets = table2array(paramsets(:,1:3));  
num_a_coefs = length(unique(paramsets(:,1)));
num_b_coefs = length(unique(paramsets(:,2)));
num_c_coefs = length(unique(paramsets(:,3)));

a_coef_str = unique(paramsets(:,1));
b_coef_str = unique(paramsets(:,2));
c_coef_str = unique(paramsets(:,3));

genvec = readtable('0_gen_vec.xlsx');
genvec = table2array(genvec);
gennum = genvec(end);
transient = 0.1*gennum;
notexist=[];
varnames = { 'Median_prey_awake_dur', 'IQR_prey_awake_dur', 'Median_pred_awake_dur', 'IQR_pred_awake_dur', ...
    'Median_prey_awake_bouts', 'IQR_prey_awake_bouts', 'Median_pred_awake_bouts', 'IQR_pred_awake_bouts', ...
    'PreyD0med','PreyD0iqr','PredD0med','PredD0iqr','PreyCHImed','PreyCHIiqr','PredCHImed','PredCHIiqr', ...
    'PreyAmed','PreyAiqr','PredAmed','PredAiqr','PreyTAUmed','PreyTAUiqr','PredTAUmed','PredTAUiqr'};
if data_struct==1
all_data_struct = [];

for i = 1:length(paramsets)
    filename = strcat(num2str(i),'_paramsetdata.xlsx');
    if exist(filename, 'file')==2; % Check that file exists
        dataset = table2struct(readtable(filename));
         all_data_struct = [all_data_struct, dataset];
        else 
        fprintf(strcat(filename,' does not exist','\n\r'));
        notexist = [notexist, i];
        T1 = table( zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),...
        zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),...
        zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),...
        zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),zeros(gennum,1),...
        'VariableNames',varnames);
        dataset = table2struct(T1); 
        all_data_struct = [all_data_struct, dataset];
    end

end

all_data_table = struct2cell(all_data_struct); %Converting whole structure to table
all_data_table = cell2mat(all_data_table);
toc
else
end

%% Metrics for plotting - average median awake duration
%Finding average median awake duration for each parameter set. 
av_prey_awake_dur = squeeze(mean(all_data_table(1,transient:gennum,:),2));
av_pred_awake_dur = squeeze(mean(all_data_table(3,transient:gennum,:),2));

pred_a_iqr = squeeze(all_data_table(20,transient:gennum,:));
prey_a_iqr = squeeze(all_data_table(18,transient:gennum,:));
Ref_matrix=[3,9, -9, 0; 8,13, 10*60, 45*60; 10,17, -1,1; 22,21 (24/(0.99729*22))^2, (24/(0.99729*26))^2]; %prey median index from table of saved values

prey_awake_dur_matrix = reshape(av_prey_awake_dur, num_c_coefs,num_b_coefs,num_a_coefs);
pred_awake_dur_matrix = reshape(av_pred_awake_dur, num_c_coefs,num_b_coefs,num_a_coefs);

%% Awake dur graphs
 for j = 1:num_b_coefs
    
    figure(12)
    subplot(3,4,j) %2*j-1);
    imagesc(a_coef_str, c_coef_str, squeeze(prey_awake_dur_matrix(:,j,:)))
    zlabel('a coefs')
    ylabel('c coefs')
    zlabel('average awake dur')
    title(strcat(num2str(b_coef_str(j)),' = b coef prey awake dur'))
    caxis([0.2,0.8])
    shading flat
    view(2)
    colorbar
    
    figure(13)
    subplot(3,4,j) %2*j)
    imagesc(a_coef_str, c_coef_str, squeeze(pred_awake_dur_matrix(:,j,:)))
    xlabel('a coefs')
    ylabel('c coefs')
    zlabel('average awake dur')
    title(strcat(num2str(b_coef_str(j)),' = b coef pred awake dur'))
    caxis([0.2,0.8])
    shading flat
    view(2)
    colorbar
 end
 
 
 %% Residual from ideal sleep durations
%  close all
awakedurresid=[];
  for j = 1:num_b_coefs
 figure(2)
 subplot(3,4,j)%ceil(num_b_coefs/2),j)
 awakedurresid(:,:,j) = sqrt(squeeze(((prey_awake_dur_matrix(:,j,:)-0.6).^2))+squeeze(((pred_awake_dur_matrix(:,j,:)-0.4).^2)));
    imagesc(a_coef_str, c_coef_str, awakedurresid(:,:,j))
    shading flat
    xlabel('prey cross int a coef')
    ylabel('pred cross int c coef') 
    title(strcat(num2str(b_coef_str(j)),' = b coef residuals'))
     caxis([0,0.4])
    colorbar
    view(2)
 %    saveas(gcf, strcat(num2str(b_coef_str(j)),'= b coef residuals.fig'))
  end    
 
%% Scanning all paramsets using 'draw now'

% {'Median_prey_awake_dur'(1), 'IQR_prey_awake_dur'(2), 'Median_pred_awake_dur', 'IQR_pred_awake_dur', ...
 %'Median_prey_awake_bouts'(5), 'IQR_prey_awake_bouts', 'Median_pred_awake_bouts', 'IQR_pred_awake_bouts', ...
 %'PreyD0med'(9),'PreyD0iqr'(10),'PredD0med','PredD0iqr','PreyCHImed'(13),'PreyCHIiqr','PredCHImed'(15),'PredCHIiqr', ...
 %'PreyAmed'(17),'PreyAiqr'(18),'PredAmed','PredAiqr'(20),'PreyTAUmed'(21),'PreyTAUiqr'(22),'PredTAUmed'(23),
 %'PredTAUiqr'(24)};
medplot =1;
iqrplot =0;
preyparamind=21; %(1,5,9 - D0,13 - chi,17 - a,21 - tau)
preyiqrind =preyparamind+1;
predparamind =preyparamind+2;
prediqrind = preyparamind+3; 

n=0;

for k=122
    n = n+1;
    paramsets(k,:)
    
    if medplot == 1
    figure(preyparamind)
    drawnow
    
    plot(all_data_table(preyparamind,:,k),'k');
    hold on
    plot(all_data_table(predparamind,:,k),'r');
    title('prey black, pred red median')
    
    hold off
    else
    end
    if iqrplot==1
%         if preyiqrind
    figure(predparamind)
     drawnow
     plot(all_data_table(preyiqrind,:,k),'k');
     hold on
     plot(all_data_table(prediqrind,:,k),'r');
     title('prey black, pred red iqr')
     hold off
    else 
    end
    pause(1)
end

%% a gene distance

% 1) Define the A-distance for a simulation: mean(abs(PredAmed(21:end) - PreyAmed(21:end)))
% That is, on average, how far apart are the median A values for the predator vs. prey populations, 
% throwing away the first 20 gens.
% Try making surface plots of this -- intermediate values are probably of most interest. 
% Values close to 0 or >1 would suggest static solutions, whereas intermediate values may have dynamics.
% 2) Could attempt step 1) for the other 3 genes. Might be harder to interpret, but worth a try.
param2graphprey = 17;
param2graphpred = 19;

% all_data_table(metric,generations,paramsetind)
%scale = max(all_data_table(param2graphprey,:,1))-min(all_data_table(param2graphprey,:,1))

av_prey_med_gene = squeeze(mean(all_data_table(param2graphprey,transient:gennum,:),2));
av_pred_med_gene = squeeze(mean(all_data_table(param2graphpred,transient:gennum,:),2));
av_dist = abs(av_prey_med_gene - av_pred_med_gene);
av_dist = squeeze(reshape(av_dist, num_c_coefs,num_b_coefs,num_a_coefs));

% Creating the data for plotting
for j = 1:num_b_coefs
 figure(1)
 subplot(2,ceil(num_b_coefs/2),j)
    surf(a_coef_str, c_coef_str, squeeze(av_dist(:,j,:)))
    zlabel('a coefs')
    ylabel('c coefs')
    zlabel('average param dist')
    caxis([0,2])
    title(strcat(num2str(b_coef_str(j)),' = b coef, a gene mean dist'))
    shading flat
    view(2)
    h = colorbar
%     set(h, 'ylim', [0,1.4])
end
%% Relative metrics
% 3) Define the A-variability for a simulation: mean(PredAiqr(21:end)) +  mean(PreyAiqr(21:end))
% 4) Could attempt step 3) for the other 3 genes.
% {'Median_prey_awake_dur'(1), 'IQR_prey_awake_dur'(2), 'Median_pred_awake_dur', 'IQR_pred_awake_dur', ...
 %'Median_prey_awake_bouts'(5), 'IQR_prey_awake_bouts', 'Median_pred_awake_bouts', 'IQR_pred_awake_bouts', ...
 %'PreyD0med'(9),'PreyD0iqr'(10),'PredD0med','PredD0iqr','PreyCHImed'(13),'PreyCHIiqr','PredCHImed'(15),'PredCHIiqr', ...
 %'PreyAmed'(17),'PreyAiqr'(18),'PredAmed','PredAiqr'(20),'PreyTAUmed'(21),'PreyTAUiqr'(22),'PredTAUmed'(23),
 %'PredTAUiqr'(24)};
param2graphprey = 18;
param2graphpred = 20;
% all_data_table(metric,generations,paramsetind)
%scale = max(all_data_table(param2graphprey,:,1))-min(all_data_table(param2graphprey,:,1))
av_prey_iqr_gene = squeeze(mean(all_data_table(param2graphprey,transient:gennum,:),2));
av_pred_iqr_gene = squeeze(mean(all_data_table(param2graphpred,transient:gennum,:),2));
%av_var = av_prey_iqr_gene + av_pred_iqr_gene;
%av_var = reshape(av_var, num_c_coefs,num_b_coefs,num_a_coefs);

av_pred_iqr_gene = reshape(av_pred_iqr_gene, num_c_coefs,num_b_coefs,num_a_coefs);
av_prey_iqr_gene = reshape(av_prey_iqr_gene, num_c_coefs,num_b_coefs,num_a_coefs);

% Creating the data for plotting
for j = 1:num_a_coefs
 figure(10)
 subplot(4,ceil(num_a_coefs/4),j)
    surf(b_coef_str, c_coef_str, av_pred_iqr_gene(:,:,j))
    zlabel('b coefs')
    ylabel('c coefs')
    zlabel('average param dist')
    title(strcat(num2str(a_coef_str(j)),' = a coef, pred av iqr'))
    shading flat
    view(2)
    h = colorbar
    set(h, 'ylim', [0,0.15])
figure(11)
 subplot(4,ceil(num_a_coefs/4),j)
    surf(b_coef_str, c_coef_str, av_prey_iqr_gene(:,:,j))
    zlabel('b coefs')
    ylabel('c coefs')
    zlabel('average param dist')
    title(strcat(num2str(a_coef_str(j)),' = a coef, prey iqr av'))
    shading flat
    view(2)
    colorbar
    h = colorbar
    set(h, 'ylim', [0,0.15])
end

%%

[gen_pred_a,n_pred_a]=find(pred_a_iqr>0.40);
osc_psets = paramsets(unique(n_pred_a),:);
figure(40)
hold on
plot3(osc_psets(:,1),osc_psets(:,2),osc_psets(:,3),'*','MarkerSize',10);
xlabel('a coef'), ylabel('b coef')



awakedurresid = sqrt(squeeze(((prey_awake_dur_matrix(:,j,:)-0.5).^2))+squeeze(((pred_awake_dur_matrix(:,j,:)-0.4).^2)));
[gen_pred_a,n_pred_awakedurres]=find(awakedurresid<0.2);
res_psets = paramsets(unique(n_pred_awakedurres),:);
figure(40)
hold on
plot3(res_psets(:,1),res_psets(:,2),res_psets(:,3),'.','MarkerSize',15);
xlabel('a coef'), ylabel('b coef'), zlabel('c coef')
title('Interesting param sets')
legend('a gene osc','low awake resids')
