clc; clear all; close all;

%read paramsetrec and find all paramsets that have b= 1.6
paramsets = readtable('0_paramsetrec.xlsx');
paramsets = table2array(paramsets(:,1:3));  
genvec = readtable('0_gen_vec.xlsx');
genvec = table2array(genvec);
gennum = length(genvec);
transient = 100; 
filename = {'_paramsetdata.xlsx'};

axis_names =  {'Med prey awake dur', 'IQR prey awake dur', 'Med pred awake dur', 'IQR pred awake dur', ...
 'Med prey awake bouts', 'IQR prey awake bouts', 'Med pred awake bouts', 'IQR pred awake bouts', ...
 'PreyD0med','PreyD0iqr','PredD0med','PredD0iqr','PreyCHImed','PreyCHIiqr','PredCHImed','PredCHIiqr', ...
 'PreyAmed','PreyAiqr','PredAmed','PredAiqr','PreyTAUmed','PreyTAUiqr','PredTAUmed',...
 'PredTAUiqr'};
%
no_params = length(paramsets);

num_a_coefs = length(unique(paramsets(:,1)));
num_b_coefs = length(unique(paramsets(:,2)));
num_c_coefs = length(unique(paramsets(:,3)));

a_coef_str = unique(paramsets(:,1));
b_coef_str = unique(paramsets(:,2));
c_coef_str = unique(paramsets(:,3));

%reading excel sheets into table

datasets = []; 
for t = 1:length(paramsets)
    
    filename1 = strcat(num2str(t),filename);
    if isfile(filename1)
    dataset = table2struct(readtable(string(filename1),'ReadVariableNames',true)); %read file for relevant columns
  
    datasets = [datasets,dataset]; 
    
    end
    
end
datasets = struct2cell(datasets); %Converting whole structure to table
datasets = cell2mat(datasets); %[gene, generations, paramset]
% datasets_plot = datasets;
datasets([13:16],:,:) = log(datasets([13:16],:,:));
%%
close all; 
varlist = {'a_coef_str','c_coef_str','gennum','num_a_coefs','paramsets',...    
'axis_names','datasets','genvec','num_b_coefs','transient',... 
'b_coef_str','filename','no_params','num_c_coefs','datasets'};
clearvars('-except',varlist{:});
%'PreyD0med'(9),'PreyD0iqr'(10),'PredD0med','PredD0iqr',
 %'PreyCHImed'(13),'PreyCHIiqr','PredCHImed'(15),'PredCHIiqr', ...
 %'PreyAmed'(17),'PreyAiqr'(18),'PredAmed','PredAiqr'(20),
 %'PreyTAUmed'(21),'PreyTAUiqr'(22),'PredTAUmed'(23),'PredTAUiqr'(24)};
 
x_col =17; %Choose col from axis names
y_col = x_col+2; %Choose col from axis names
lag = 0:100;
desc = strcat(axis_names(x_col),' vs', {' '}, axis_names(y_col));

%
cor_rec_detrend = [];
for n = 1:length(paramsets) 
   
        x_data = datasets(x_col,transient:end,n);
        y_data = datasets(y_col,transient:end,n);
         
        time_vec = [1:length(x_data)];
        p_x = polyfit(time_vec, x_data,1);
        p_y = polyfit(time_vec, y_data,1);

        detrend_x = polyval(p_x,time_vec);
        detrend_y = polyval(p_y,time_vec);

        x_data_0mean = x_data-detrend_x;
        y_data_0mean = y_data-detrend_y;
 
        x_data_std(n) = std(x_data_0mean);
        y_data_std(n) = std(y_data_0mean);
        
        for g = 1:length(lag)
        x_data_lag = x_data_0mean(1:(end-lag(g)));
        y_data_lag = y_data_0mean(1+lag(g):end);
        
        cor_rec_detrend(n,g) = corr(x_data_lag',y_data_lag');
 
        end
        
        cor_rec_detrend_smooth(n,:)= movmean(cor_rec_detrend(n,:),5); %smoothing correlation data - moving average window size of 5.

    clear dataset
    clear filename1
    clear x_data
    clear y_data

    
end

% Conds for finding oscs - conditions by file name index
[max_cor_val,max_cor_lag] = max(cor_rec_detrend_smooth'); %USING SMOOTHED COR DATA FOR EVERYTHING

% CONDITION 1: minimum max value ~0.4 
cond1= find(max_cor_val>0.4); %%%

% CONDITION 2: at least one peak anywhere (but not at 0 lag)
% Look for where the slope crosses from positive to negative between adjacent points. 
%(Hopefully there are only 0 or 1 instances per paramset)

local_max = [];
local_max_num = []; 
for t= 1:length(paramsets) 
        cor_der = diff(cor_rec_detrend_smooth(t,:));
        local_max{t} = (find([0 diff(sign(cor_der))]<0)); %finds local maximums of smoothed correlation curve
        local_max_num(t) = numel(cell2mat((local_max(t)))); %finding number of local maximums
end

cond2 = find(local_max_num~=0); %%%

% CONDITION 3: maximum value not at 0 lag - are there any at 0? 

cond3 = find(max_cor_lag >1); %%% %not 0 because this is actually the index of lag vector

% CONDITION 4: cor range is greater than ~0.6

cond4 = find(range(cor_rec_detrend_smooth')>0.55); %%%

% CONDITION 5: minimum variation/std deviation
% range of genes
mph=60;
hpd=24;

rang(1,1)=-30; rang(1,2) =13; %D0 range [-30,5] %%keep at -30
rang(2,1) = log(5*mph); rang(2,2) = log(45*mph); %chi range
rang(3,1) = -1.0; rang(3,2) =1.0; %a range
rang(4,2) = (hpd/(0.99729*22))^2; rang(4,1) = (hpd/(0.99729*26))^2; %tau range

Ref_ind= [9,1;13,2;17,3;21,4];

ind = find(x_col==Ref_ind(:,1));
range_rec = abs(rang(ind,2)-rang(ind,1));

x_std_perc = x_data_std*100/range_rec;
y_std_perc = y_data_std*100/range_rec;

std_thresh =12; %%%
cond5 = find(x_std_perc>std_thresh& y_std_perc>std_thresh);


% intersection of conditions (file_inds)

cond12 = intersect(cond1,cond2);
cond123 = intersect(cond12, cond3);
cond1234 = intersect(cond123, cond4);
cond12345_all = intersect(cond1234, cond5);

% text that goes on graphs
text_cond1 = 'max cor is ';
text_cond2 = ', #peaks = ';
text_cond3 = ', max cor @ lag ';
text_cond4 = ', cor range is '; 
text_cond5_pry = ' prey std ';
text_cond5_prd = ' pred std ';

cond1_vals = max_cor_val(cond12345_all);
cond2_vals = local_max_num(cond12345_all);
cond3_vals = max_cor_lag(cond12345_all); 
cond4_vals = range(cor_rec_detrend(cond12345_all,:),2);
cond5_val_pry = x_std_perc(cond12345_all);
cond5_val_prd = y_std_perc(cond12345_all);


failed_inds = 1:no_params;
failed_inds(cond12345_all)=0;

failed_inds = failed_inds(failed_inds~=0);

cond1_vals_fail = max_cor_val(failed_inds);
cond2_vals_fail = local_max_num(failed_inds);
cond3_vals_fail = max_cor_lag(failed_inds); 
cond4_vals_fail = range(cor_rec_detrend(failed_inds,:),2);
cond5_val_pry_fail = x_std_perc(failed_inds);
cond5_val_prd_fail = y_std_perc(failed_inds);


%% plotting
%switches
gene_osc_switch = 1;
gene_no_osc_switch = abs(gene_osc_switch-1);  
contour_plot_switch = 0; 



% PLOTS
close all

ind = find(x_col==Ref_ind(:,1));
ylims = [rang(ind,1),rang(ind,2)];

%oscilation plots
if gene_osc_switch==1 

     for m = 1:length(cond12345_all)
    dataset = datasets(:,:,cond12345_all(m));
        txt = strcat(text_cond1, num2str(cond1_vals(m)), text_cond2, num2str(cond2_vals(m)), text_cond3, num2str(cond3_vals(m)), text_cond4, num2str(cond4_vals(m)),...
             text_cond5_pry, num2str(cond5_val_pry(m)),text_cond5_prd, num2str(cond5_val_prd(m)));
        h = figure;
        h.WindowState = 'maximized';
        subplot(2,1,1)
        plot(dataset(x_col,:))
        hold on; plot(dataset(y_col,:)) 
         ylim(ylims)
        legend(axis_names{x_col},axis_names{y_col})
        title(strcat(num2str(cond12345_all(m)),' filename',desc))
        text(1, min([dataset(x_col,:),dataset(y_col,:)]), txt)

        subplot(2,1,2)
        plot(lag,cor_rec_detrend(cond12345_all(m),:))
        hold on; plot(lag,cor_rec_detrend_smooth(cond12345_all(m),:))
        ylim([-1,1])
        title(strcat(num2str(cond12345_all(m)),' filename cor plot ',desc))
      
    end
end


% failed criteria

%no oscilation plots
if gene_no_osc_switch==1 
    for k =40:50%length(failed_inds)

        dataset = datasets(:,:,failed_inds(k));
        txt = strcat(text_cond1, num2str(cond1_vals_fail(k)), text_cond2, num2str(cond2_vals_fail(k)), text_cond3, num2str(cond3_vals_fail(k)), text_cond4, num2str(cond4_vals_fail(k)), ...
        text_cond5_pry, num2str(cond5_val_pry_fail(k)),text_cond5_prd, num2str(cond5_val_prd_fail(k)));
        h = figure;
        h.WindowState = 'maximized';
        subplot(2,1,1)
        plot(dataset(x_col,:))
        hold on; plot(dataset(y_col,:)) 
        ylim(ylims)
        legend(axis_names{x_col},axis_names{y_col})
        title(strcat(num2str(failed_inds(k)),' filename',desc))
        text(1, min([dataset(x_col,:),dataset(y_col,:)]), txt)

        subplot(2,1,2)
        plot(lag,cor_rec_detrend(failed_inds(k),:))
        ylim([-1,1])
        hold on; plot(lag,cor_rec_detrend_smooth(failed_inds(k),:))
        title(strcat(num2str(failed_inds(k)),' cor rec',desc))
        
    end
end


% contour plot

if contour_plot_switch==1
    
    osc_map = zeros(1,no_params);
    osc_map([cond12345_all])=1; %
    oscs_shaped = squeeze(reshape(osc_map, num_c_coefs,num_a_coefs));

    figure
    imagesc(a_coef_str, c_coef_str, squeeze(oscs_shaped)); hold on
    colormap(parula); hold on    
    xlabel('a coefs')
    ylabel('c coefs')
    title('b = 1.6 coef osc. map')
end
