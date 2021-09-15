%% all gene freq figure

%load the awakestate data for correct parameter set
%load paramsetdata_all
% close all; 
k =52; %coefficient set

preyparamind=[9,13,17,21]; %(1,5,9 - D0,13 - chi,17 - a,21 - tau)
predparamind =preyparamind+2;

% figure position details
figrows = 6;
figcols = 1;
boxes = figrows*figcols;

%gen values - prey^&pred^, prey_&pred^, prey_&pred_, prey^&pred_, prey^&pred^
% start_gen = 339;
% end_gen = 275;
% gen_values = floor(linspace(start_gen,end_gen,figcols));
gen_values = [137,187, 285,367,387 ];
%
genes = [9,13,17,21];
gene_names = {'D0','chi','a','tau'};

%Ref_matrix = [par ind, table ind, minimum par value, maximum par value]
Ref_matrix=[3,9, -30, 13; ...
    8,13, 5*60, 45*60; ...
    10,17, -1,1; ... 
    22,21, (24/(0.99729*26))^2,(24/(0.99729*22))^2]; %prey median index from table of saved values
       
% mph=60;
% hpd=24;
% 
% 
% rang(1,1)=-30; rang(1,2) =13; %D0 range [-30,5] %%keep at -30
% rang(2,1) = log(5*mph); rang(2,2) = log(45*mph); %chi range
% rang(3,1) = -1.0; rang(3,2) =1.0; %a range
% rang(4,2) = (hpd/(0.99729*22))^2; rang(4,1) = (hpd/(0.99729*26))^2; %tau range
% 
% Ref_ind= [9,1;13,2;17,3;21,4];


gens2disp = [0,gennum];
% gene frequency plots

for t=1:length(preyparamind)
    [val,ind] = find(preyparamind(t) == genes);
    title_gene = char(gene_names(ind));
    figure(20)
    subplot(figrows,figcols,t)
    plot(all_data_table(preyparamind(t),:,k),'k');
    hold on
    plot(all_data_table(predparamind(t),:,k),'r');
    rectangle('Position', [gen_values(1) Ref_matrix(ind,3) (gen_values(2)-gen_values(1)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
    rectangle('Position', [gen_values(2) Ref_matrix(ind,3) (gen_values(3)-gen_values(2)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
    rectangle('Position', [gen_values(3) Ref_matrix(ind,3) (gen_values(4)-gen_values(3)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
    rectangle('Position', [gen_values(4) Ref_matrix(ind,3) (gen_values(5)-gen_values(4)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
    
    ylim([Ref_matrix(ind,3),Ref_matrix(ind,4)])
    title([title_gene ,' gene frequency' ])
    set(gca,'xtick',[])
    legend('prey','pred','Location','northwest')
    xlim([gens2disp])
end

%median awake duration
ylimmax = 1;
figure(20)
subplot(figrows,figcols,5)
plot(all_data_table(1,:,k),'k');hold on
plot(all_data_table(3,:,k),'r');hold on
rectangle('Position', [gen_values(1) 0 (gen_values(2)-gen_values(1)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(2) 0 (gen_values(3)-gen_values(2)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(3) 0 (gen_values(4)-gen_values(3)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(4) 0 (gen_values(5)-gen_values(4)) ylimmax], 'EdgeColor', 'b' )

ylim([0,ylimmax])
xlim([gens2disp])
set(gca,'xtick',[])
title(['median awake duration' ])
legend('prey','pred','Location','northwest')

%median awake bouts
ylimmax1 = max(all_data_table(5,:,k));
ylimmax2 = max(all_data_table(7,:,k));
ylimmax = max([ylimmax1,ylimmax2])+1;

figure(20)
subplot(figrows,figcols,6)
plot(all_data_table(5,:,k),'k');hold on
plot(all_data_table(7,:,k),'r');hold on
rectangle('Position', [gen_values(1) 0 (gen_values(2)-gen_values(1)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(2) 0 (gen_values(3)-gen_values(2)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(3) 0 (gen_values(4)-gen_values(3)) ylimmax], 'EdgeColor', 'b' )
rectangle('Position', [gen_values(4) 0 (gen_values(5)-gen_values(4)) ylimmax], 'EdgeColor', 'b' )
xlim([gens2disp])
ylim([0,ylimmax])
title(['median awake bouts' ])
legend('prey','pred','Location','northwest')
xlabel('generations')
% set(gca,'xtick',[])

% %%
% figure(20)
% subplot(figrows,figcols,(figcols*(7-1)+1):(figcols*7-1))
% for t = 1:gennum
% cor(t)=corr(avprey(:,t),avpred(:,t));
% end
% plot(cor,'k');hold on
% rectangle('Position', [gen_values(1) -1 (gen_values(end)-gen_values(1)) ylimmax], 'EdgeColor', 'b' )
% xlim([gens2disp])
% ylim([-1,1])
% xlabel('generations')
% title(['awake state correlation plot'])

%% subplot timeseries
figrows = 6;
figcols = 2;

linth = 0:0.01:3*24;
Imax = 1;%ones(size(linth));
dawn = 6;
dusk = 18;
I = zeros(size(linth));
I = Imax*(mod(linth,24)>dawn).*(mod(linth,24)<=dusk);

for c = 1:length(gen_values)
    hs = figure(21)
    subplot(figrows,figcols,2*c-1)
    area(linth,I,'FaceColor', 'y','FaceAlpha',0.5, 'LineStyle','none');hold on
    plot(linth, avprey(gen_values(c),:),'k'); hold on
    plot(linth, avpred(gen_values(c),:),'r'); hold on
    plot(linth, fitprey(gen_values(c),:)+0.005,'k'); hold on
    plot(linth, fitpred(gen_values(c),:)-0.005,'r'); hold on
        title(['generation', ' ', num2str(gen_values(c))])

    xlabel('time hrs')
    ylabel('awake state')
    
    ylim([-0.1,1.1])
    xlim([0,72])
%     title()
end

% raster plots

%Creating raster vector

raster = fitprey*2 + fitpred;
raster(raster==1)=0;
raster(raster==2)=1;
raster(raster==3)=2;

for c = 1:length(gen_values)
    hs = figure(21)
    subplot(figrows,figcols,2*c)
    imagesc(raster(gen_values(c),:)); hold on
    caxis([0,2])

    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

subplot(figrows,figcols,11:12)
axis(gca,'off')
colorbar('Ticks',[0,0.5,1],'TickLabels',{'fittest prey asleep','only fittest prey awake','fittest pred and fittest prey awake'}, ...
        'Location','southoutside')
colormap(parula(3))

% set(gcf, 'Position', get(0, 'Screensize'));
% filename = num2str(k)
% filename = strcat(filename,'_raster.png')
% saveas(gcf,filename)
%% Making data strucutres/ data workspaces 
clear all

for t = [25, 26, 52, 66]
paramsets = readtable('0_paramsetrec.xlsx');
paramsets = table2array(paramsets(:,1:3));  

genvec = readtable('0_gen_vec.xlsx');
genvec = table2array(genvec);
gennum = genvec(end);

% filenames string
filenames = {'_avpreyavawakestate.xlsx','_avpredavawakestate.xlsx', '_fittestpreyawakestate.xlsx','_fittestpredawakestate.xlsx'};
all_files = [];
    for f = 1:length(filenames)
    all_files = [all_files, strcat(num2str(t),filenames(f))];
    end
    
% creating tables
tic
avprey = table2array(readtable(char(all_files(1))));
fprintf('avprey done\n')
toc
avpred = table2array(readtable(char(all_files(2))));
fprintf('avpred done\n')
toc
fitprey = table2array(readtable(char(all_files(3))));
fprintf('fitprey done\n')
toc
fitpred = table2array(readtable(char(all_files(4))));
fprintf('fitpred done\n')
toc 

workspace_name = strcat(num2str(t),'_awakedata.mat')
save(workspace_name)
clear all

end

%% calculating fitness values
k = 37;
a = paramsets(k,1)
b = paramsets(k,2)
c = paramsets(k,3)
for i = 1:500;
 
    preyusefulawake(i) = length(find(fitprey(i,:)==1))/length(fitprey(i,:));
    preyusefulsleep(i) = length(find(fitprey(i,:)==0))/length(fitprey(i,:));

    predusefulawake(i) = length(find(fitpred(i,:)==1))/length(fitpred(i,:));
    predusefulsleep(i) = length(find(fitpred(i,:)==0))/length(fitpred(i,:));

    preyfit(i) = (-max(0,mean(fitprey(i,:).*avprey(i,:))-a*mean(fitprey(i,:).*avpred(i,:))))*(preyusefulsleep(i).^b*preyusefulawake(i));
    predfit(i) = -sqrt(mean(fitpred(i,:).*avpred(i,:)))*(predusefulsleep(i).^c)*(mean(fitpred(i,:).*avprey(i,:))); %this is to be minimised  
end 

%% one gene freq figure

%load the awakestate data for correct parameter set
%load paramsetdata_all

k = 4; %coefficient set

preyparamind=9 %(1,5,9 - D0,13 - chi,17 - a,21 - tau)
predparamind =preyparamind+2;

% figure position details
figrows = 3;
figcols = 3;
boxes = figrows*figcols;

%gen values 
start_gen = 244;
end_gen = 257;
gen_values = [244,250,257];%floor(linspace(start_gen,end_gen,figcols));
%

genes = [9,13,17,21];
gene_names = {'D0','chi','a','tau'};
[val,ind] = find(preyparamind == genes);
title_gene = char(gene_names(ind));

%Ref_matrix = [par ind, table ind, minimum par value, maximum par value]
Ref_matrix=[3,9, -9, 0; ...
    8,13, 10*60, 45*60; ...
    10,17, -1,1; ... 
    22,21, (24/(0.99729*26))^2,(24/(0.99729*22))^2]; %prey median index from table of saved values

%main plot
hs = figure(20)
subplot(figrows,figcols,[1:(figrows*(figcols-1))])
plot(all_data_table(preyparamind,:,k),'k');
hold on
plot(all_data_table(predparamind,:,k),'r');
rectangle('Position', [gen_values(1) Ref_matrix(ind,3) (gen_values(end)-gen_values(1)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
ylim([Ref_matrix(ind,3),Ref_matrix(ind,4)])
    title([title_gene ,' gene frequency' ])
xlabel('generations')
legend('prey','pred','Location','northwest')
