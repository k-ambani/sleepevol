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

%% filenames that have oscillations
oscs_D0 = [38, 51, 52, 53, 64, 79, 80]; % 7 paramsets
oscs_chi = [25, 26, 39]; %3 paramsets
oscs_a = [38, 51, 52, 53,...
            64, 65, 66, 79, 80,...
            82, 93, 94]; %14 paramsets
oscs_tau = [25, 26, 38, 39, 51,...
            52, 53, 64, 65, 66,...
            67, 77, 78, 79, 80,...
            91, 92, 93, 94, 95,...
            103, 105]; %24 paramsets
        
 
% gene ranges
       
mph=60;
hpd=24;


rang(1,1)=-30; rang(1,2) =13; %D0 range [-30,5] %%keep at -30
rang(2,1) = log(5*mph); rang(2,2) = log(45*mph); %chi range
rang(3,1) = -1.0; rang(3,2) =1.0; %a range
rang(4,2) = (hpd/(0.99729*22))^2; rang(4,1) = (hpd/(0.99729*26))^2; %tau range

Ref_ind= [9,1;13,2;17,3;21,4];

%%

%'PreyD0med'(9),'PreyD0iqr'(10),'PredD0med','PredD0iqr',
 %'PreyCHImed'(13),'PreyCHIiqr','PredCHImed'(15),'PredCHIiqr', ...
 %'PreyAmed'(17),'PreyAiqr'(18),'PredAmed','PredAiqr'(20),
 %'PreyTAUmed'(21),'PreyTAUiqr'(22),'PredTAUmed'(23),'PredTAUiqr'(24)};
 
x_col =17; %Choose col from axis names
y_col = x_col+2; %Choose col from axis names
desc = strcat(axis_names(x_col),' vs', {' '}, axis_names(y_col));
ind = find(x_col==Ref_ind(:,1));
ylims = [rang(ind,1),rang(ind,2)];

%%
close all; clc;

cond12345_all = oscs_a;


for m = 1:length(cond12345_all)
    dataset = datasets(:,:,cond12345_all(m));
        
        h = figure;
        h.WindowState = 'maximized';
        plot(dataset(x_col,:))
        hold on; plot(dataset(y_col,:)) 
         ylim(ylims)
        legend(axis_names{x_col},axis_names{y_col})
        title(strcat(num2str(cond12345_all(m)),' filename',desc))
        
 end
 