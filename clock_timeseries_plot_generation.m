%code that generates cheese_clock, gene_freq and pop_awake_states plots
%enter the k value before running code. 

clc;
clear all;
close all; 
k =91; %manually entercoefficient set



if k == 26
load('26_awakedata.mat')
load('paramset_data.mat')

gen_values = [203,210,215,227,230,235,238,243];
gens2disp = [100,450];

end

if k == 38
load('38_awakedata.mat')
load('paramset_data.mat')

gen_values = [220,230,270,340,430,445];
gens2disp = [100,500];

end

if k == 52
load('52_awakedata.mat')
load('paramset_data.mat')

gen_values = [136,140,160,200,215,350,375,385];
gens2disp = [100,450];

end


if k == 91
load('91_awakedata.mat')
load('paramset_data.mat')
gen_values = [250,255,260,265,270,275,277];
gens2disp = [200,500];

end




preyparamind=[9,13,17,21]; %(1,5,9 - D0,13 - chi,17 - a,21 - tau)
predparamind =preyparamind+2;
linth = 0:0.01:3*24;

%
genes = [9,13,17,21];
gene_names = {'D0','chi','a','intrinsic period'};


%% transforming tau gene - run this only once
all_data_table(preyparamind(end),:,:) = 24/(0.99729*sqrt(all_data_table(preyparamind(end),:,:)));
all_data_table(predparamind(end),:,:) = 24/(0.99729*sqrt(all_data_table(predparamind(end),:,:)));

%%
Ref_matrix=[3,9, -30, 13; ...
    8,13, 5*60, 45*60; ...
    10,17, -1,1; ... 
    22,21, 22,25.5]; %prey median index from table of saved values

figrows = length(gen_values);
figcols = 1;

%% day and night vectors

dawn = 6;
dusk = 18;

Imax = 1;
Imax_clock = 0.35;%ones(size(linth));

I_day = zeros(size(linth));

I_day = Imax*(mod(linth,24)>dawn).*(mod(linth,24)<=dusk);
x_day = I_day.*cos(linth*pi/12-pi/2).*Imax_clock; 
y_day = I_day.*sin(linth*pi/12-pi/2).*Imax_clock; 


I_night = I_day==0;
I_night = Imax.*I_night;
x_night = I_night.*cos(linth*pi/12-pi/2).*Imax_clock; 
y_night = I_night.*sin(linth*pi/12-pi/2).*Imax_clock; 


%% Time series plots

for c = 1:length(gen_values)
    h_ts = figure(24)
    subplot(figrows,figcols,c)
    area(linth,I_day,'FaceColor', 'y','FaceAlpha',0.2, 'LineStyle','none');hold on
    area(linth,I_night,'FaceColor', 'b','FaceAlpha',0.1, 'LineStyle','none');hold on

    plot(linth, avprey(gen_values(c),:),'k'); hold on
    plot(linth, avpred(gen_values(c),:),'r'); hold on

    ylabel(['gen', ' ', num2str(gen_values(c))],'FontWeight','bold')
    xticks([0:6:72])
    xticklabels([0:6:18,0:6:18,0:6:18,0:6:18,24])
    ylim([-0.1,1.1])
    xlim([0,72])
%     title()
end
    xlabel('time hrs')


saveas(h_ts, strcat(num2str(k),'_pop_awake_states'))
%% clock plots
% function for plotting vector

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),varargin{:});

%transformations for nicer figures
[prey_v_ang,prey_v_mag] = cart2pol(real(vecsumprey),imag(vecsumprey));
[pred_v_ang,pred_v_mag] = cart2pol(real(vecsumpred),imag(vecsumpred));

prey_v_ang_t = -prey_v_ang-pi/2;
pred_v_ang_t = -pred_v_ang-pi/2;


[vecsumprey_real,vecsumprey_imag] = pol2cart(prey_v_ang_t,prey_v_mag);
[vecsumpred_real,vecsumpred_imag] = pol2cart(pred_v_ang_t,pred_v_mag);

vecsumprey_compl = complex(vecsumprey_real,vecsumprey_imag);
vecsumpred_compl = complex(vecsumpred_real,vecsumpred_imag);

for g=1:figrows
    av_vecsum_prey = sum(vecsumprey_compl(:,gen_values(g))/720100);
    av_vecsum_pred = sum(vecsumpred_compl(:,gen_values(g))/720100);

    h_cc = figure(1)
    subplot(2,figrows,g)
    drawArrow([0,real(av_vecsum_prey)],[0,imag(av_vecsum_prey)],'MaxHeadSize',7,'linewidth',5,'color',[0.10,0.10,0.10]); hold on
    scatter(vecsumprey_real(:,gen_values(g))/length(linth),vecsumprey_imag(:,gen_values(g))/length(linth),6, 'filled',...
    'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none'); hold on;  
    patch(x_day,y_day,'yellow','FaceAlpha',0.2,'EdgeColor','none'); hold on;
    patch(x_night,y_night,'blue','FaceAlpha',0.2,'EdgeColor','none')
    text(0,0.38,'midday','HorizontalAlignment','center')
    h1 = text(0.4,0,'dusk','HorizontalAlignment','center')
    set(h1,'Rotation',270)
    text(0,-0.38,'midnight','HorizontalAlignment','center')
    h2 = text(-0.4,0,'dawn','HorizontalAlignment','center')
    set(h2,'Rotation',90)

    xlim([-0.4,0.4])
    ylim([-0.4,0.4])
    title(strcat('gen ', num2str(gen_values(g))))
    set(gca,'Xtick',[])
    set(gca,'ytick',[])
    %    set(gca,'Visible','off')
    set(gca,'XColor','none')
    set(gca,'YColor','none')

    h_cc = figure(1)
    subplot(2,figrows,g+figrows)
      drawArrow([0,real(av_vecsum_pred)],[0,imag(av_vecsum_pred)],'MaxHeadSize',3,'linewidth',5,'color',[1, 0.3250, 0]); hold on    
    scatter(vecsumpred_real(:,gen_values(g))/length(linth),vecsumpred_imag(:,gen_values(g))/length(linth),6, 'filled',...
    'MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none'); hold on;
    patch(x_day,y_day,'yellow','FaceAlpha',0.2,'EdgeColor','none')
    patch(x_night,y_night,'blue','FaceAlpha',0.2,'EdgeColor','none')
    text(0,0.38,'midday','HorizontalAlignment','center')
    h3 = text(0.4,0,'dusk','HorizontalAlignment','center')
    set(h3,'Rotation',270)
    text(0,-0.38,'midnight','HorizontalAlignment','center')
    h4 = text(-0.4,0,'dawn','HorizontalAlignment','center')
    set(h4,'Rotation',90)

    xlim([-0.4,0.4])
    ylim([-0.4,0.4])
    %     set(gca,'Visible','off')
    set(gca,'Xtick',[])
    set(gca,'ytick',[])
    set(gca,'XColor','none')
    set(gca,'YColor','none')
end

saveas(h_cc,strcat(num2str(k),'_cheese_clock.fig'))
%% gene freq plot


% figure position details
figrows = 6;
figcols = 1;
boxes = figrows*figcols;
% gens2disp = [0,gennum];
% gens2disp = [gen_values(1)-10,gen_values(end)+10];

for t=1:length(preyparamind)
    [val,ind] = find(preyparamind(t) == genes);
    title_gene = char(gene_names(ind));
    h_gf = figure(20)
    subplot(figrows,figcols,t)
    plot(all_data_table(preyparamind(t),:,k),'k');
    hold on
    plot(all_data_table(predparamind(t),:,k),'r');
    
    for r= 1:(length(gen_values)-1)
    rectangle('Position', [gen_values(r) Ref_matrix(ind,3) (gen_values(r+1)-gen_values(r)) (Ref_matrix(ind,4)-Ref_matrix(ind,3))], 'EdgeColor', 'b' )
    end
    ylim([Ref_matrix(ind,3),Ref_matrix(ind,4)])
    title([title_gene ,' gene frequency' ])
    set(gca,'xtick',[])
    legend('prey','pred','Location','northwest')
    xlim([gens2disp])
end
%
%median awake duration
ylimmax = 1;
h_gf = figure(20)
subplot(figrows,figcols,5)
plot(all_data_table(1,:,k),'k');hold on
plot(all_data_table(3,:,k),'r');hold on
for r=1:(length(gen_values)-1)
rectangle('Position', [gen_values(r) 0 (gen_values(r+1)-gen_values(r)) ylimmax], 'EdgeColor', 'b' )
end
ylim([0,ylimmax])
xlim([gens2disp])
set(gca,'xtick',[])
title(['median awake duration' ])
legend('prey','pred','Location','northwest')

%median awake bouts
ylimmax1 = max(all_data_table(5,:,k));
ylimmax2 = max(all_data_table(7,:,k));
ylimmax = max([ylimmax1,ylimmax2])+1;

h_gf = figure(20)
subplot(figrows,figcols,6)
plot(all_data_table(5,:,k),'k');hold on
plot(all_data_table(7,:,k),'r');hold on
for r = 1:(length(gen_values)-1)
rectangle('Position', [gen_values(r) 0 (gen_values(r+1)-gen_values(r)) ylimmax], 'EdgeColor', 'b' )
end
xlim([gens2disp])
ylim([0,ylimmax])
title(['median awake bouts' ])
legend('prey','pred','Location','northwest')
xlabel('generations')

saveas(h_gf,strcat(num2str(k),'_gene_freq.fig'))