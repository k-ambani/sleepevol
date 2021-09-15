function X = evol_fun_para_ts(n)
  % Simulate evolution as a function, taking i as an input.
  % For n = 0, parameter set are saved. 
  if ischar(n)==1 % Check if it's a string
    n=str2num(n); % If so, convert to a number
else
end

% ----- Parameters used by sleepde.m are stored in a vector ----
%Graph options available
%Ref parameters give min and max values for alleles in gene pool.
format long
tic
X = [];
%% Specify and select fitness function parameters

P = []; % Initialize matrix of parameter values
as = [0.2:0.2:2.4];
bs = 1.6;
cs = [0.5:0.5:6];

for i1=1:length(as)
    for i2=1:length(bs)
        for i3=1:length(cs)
            P = [P; as(i1),bs(i2),cs(i3)];
        end
    end
end

%Matrix of parameter sets has been formed. Need to break into different
%arrays. 
if n == 0
else
a = P(n,1)
b = P(n,2)
c = P(n,3)
end
%% Run settings (Duration, cont, merge)
% ----- Run settings -----
dur    = 10*24;            % run duration (hours)
ts     = 10;                % Initial # data points per time unit (min^{-1})
mph=60;
hpd=24;
% This parameter is only used by the 'rgk4' solver
dur = dur*mph;
started = datestr(now)
%Simulation parameters - enter these numbers
gennum =500; %number of generations
fittest = 50; %number of fittest members KEEP EVEN
rep = 2; %popsize/fittest, number of replicants per fittest mem (offspring)
preypopsize =fittest*rep;  %number of prey
predpopsize = preypopsize; %number of predators
avmutationpergen = preypopsize*1.0;
buffermins = 10;
mutatingparams = [3,8,10,22]; %define the parameters you want to mutate

if n==0
filename1 = strcat(num2str(n),'_gen_vec.xlsx')
T3 = table([1:gennum]','VariableNames',{'Gen_vec'})
writetable(T3, filename1)
filename2 = strcat(num2str(n),'_paramsetrec.xlsx')
						     headings = { 'a','b','c'}
						     T2 = table(P(:,1),P(:,2),P(:,3),'VariableNames', headings)
						     writetable(T2, filename2)
 else
    

%% Sleep/wake switch parameters
Ref = zeros(29,3);  %[min, max, increment]

%%%%%connection strength/tau (mV/min) to VLPO from MA NOT USED
Ref(1,1)  = -2.1/(mph);       Ref(1,2) = -2.0/(mph);        Ref(1,3) = -2.2/(mph) ;
%%%%% connection strength/tau (mV/min) to VLPO from DMH
  Ref(2,1) = -0.17/(mph);       Ref(2,2) =  -0.20/(mph);      Ref(2,3) = -0.13/(mph);
%%%%% D0 (-8,0) gives awakedur from ~0.45 to ~0.66 - controls total amount of sleep
Ref(3,1) = -4.8;              Ref(3,2) = -30;                Ref(3,3) =13;
%%%%% connection strength to VLPO from Homeostatic drive mu_vh (mV nM^{-1} min^{-1})
  Ref(4,1)  = 1;                Ref(4,2) = 0.7;               Ref(4,3) = 1.3;
%%%%% connection strength/tau (mV/m) to MA from VLPO NOT USED
  Ref(5,1)  = -1.8/(mph);       Ref(5,2)  = -1.8/(mph);       Ref(5,3)  = -1.8/(mph);
%%%%% 'A0' average effect from cholinergic and other sources to MA(mV/min)
Ref(6,1) = 1.3;               Ref(6,2) = 1.0;               Ref(6,3) = 1.6;

%Parameter values for (chi)dH/dt + H = mu*Qm
			     % mu (nM min)  %%%%4.4 or 7? in original code
			     Ref(7,1) = 4.4/mph;    Ref(7,2) = 4.4/mph;    Ref(7,3) = 4.4/mph;
%%%%%chi (min) 'characteristic time for somnogen acculumation and
% clearance' (0.01-45*mph) have default at 1 because otherwise may
% cause negative chi values when adding noise
Ref(8,1) = 20*mph;     Ref(8,2) = 5*mph;     Ref(8,3) = 45*mph;
%Modulation by SPZ, etc
					       % mu_md (mV s)
					       Ref(9,1) = .01/(mph);         Ref(9,2) = .01/(mph);         Ref(9,3) = .01/(mph);
% a
Ref(10,1) = 1;                Ref(10,2) = -1.0;                Ref(10,3)=1.0;
% k (min^{-1})
Ref(11,1) = 17*mph;           Ref(11,2) = 17*mph;           Ref(11,3) = 17*mph;
% delta
Ref(12,1) = 2.8;              Ref(12,2) = 2.8;              Ref(12,3) = 2.8;
% b (min^{-1})
Ref(13,1) = 4.8*mph;          Ref(13,2) = 4.8*mph;          Ref(13,3) = 4.8*mph;

%Masking parameter
% muvb (masking effect of light on VLPO)
  Ref(14,1) =0*-880/(mph);      Ref(14,2)= 0*-880/(mph);      Ref(14,3) =0*-880/(mph);
%Circadian parameter values
% Parameters used by sleepde.m are stored in a vector. Parameter values
% are from St. Hilaire et al. Addition of a non-photic component to a
% light-based mathematical model the human circadian pacemaker, J.
								% Theor. Biol. (2007).

% oscillator stiffness
  Ref(15,1) = 0.13;      Ref(15,2) = 0.13;      Ref(15,3) = 0.13;
% c1
Ref(16,1) = 1/3;              Ref(16,2) = 1/3;              Ref(16,3) = 1/3;
% c2
Ref(17,1) = 4/3;              Ref(17,2) = 4/3;              Ref(17,3) = 4/3;
%c3
Ref(18,1) = 0;                Ref(18,2) = 0;                Ref(18,3) = 0;
%c4
Ref(19,1) = -256/105;         Ref(19,2) = -256/105;         Ref(19,3) = -256/105;
%q
Ref(20,1) = 1/3;              Ref(20,2) = 1/3;              Ref(20,3) = 1/3;
% h photic drive strength
Ref(21,1) = -0.2;             Ref(21,2) = -0.2;             Ref(21,3) = -0.2;
%%%%%intrinsic period for DD conditions, period of x (0.9,1.1)
							      Ref(22,1) = (hpd/(0.99729*24.15))^2;Ref(22,2) = (hpd/(0.99729*22))^2; Ref(22,3) = (hpd/(0.99729*26))^2;

%Idrive parameter values ('photic drive')
  % beta (min^{-1}) rate of change of activated to ready photoreceptors
  Ref(23,1) = 0.007;     Ref(23,2) = 0.007;     Ref(23,3) = 0.007;
% rho used in NS NOT USED
Ref(24,1) = 1*0.032;         Ref(24,2) = 1*0.032;          Ref(24,3) = 1*0.032;
%'C' constant in alpha, alpha0/(I0)^p
							     Ref(25,1) = 0.1/9500^0.5;    Ref(25,2) = 0.1/9500^0.5;      Ref(25,3) = 0.1/9500^0.5;
%%%%%'G' constant in Bhat
Ref(26,1) = 37;              Ref(26,2) = 37;                 Ref(26,3) = 74;
%'p'
Ref(27,1) = 0.5;      Ref(27,2) = 0.5;       Ref(27,3) =0.5;
% I_1 (lux)
Ref(28,1) = 100;      Ref(28,2) = 100;       Ref(28,3) = 100;
%modulation 'r' in paper
Ref(29,1) = 0.4;            Ref(29,2) = 0.4;            Ref(29,3)= 0.4;

pars = Ref(:,1);
%%
% Creating allele pool
allele3 = linspace(Ref(3,2),Ref(3,3),preypopsize);
allele8 = linspace(Ref(8,2),Ref(8,3),preypopsize);
allele10 = linspace(Ref(10,2),Ref(10,3),preypopsize);
allele22 = linspace(Ref(22,3),Ref(22,2),preypopsize);

%% Initialising
gencoords = repmat((1:preypopsize),gennum,1);
%preypopgenome = repmat(pars,1,preypopsize);
preyno = [1:preypopsize]';
%predpopgenome = repmat(pars,1,predpopsize);
predno = [1:predpopsize]';
preyfit=zeros(preypopsize,gennum);
predfit=zeros(predpopsize,gennum);


%Initialising results vectors
avpreyawakestate = zeros(gennum,7201);
avpredawakestate = zeros(gennum,7201);
preyconstate = zeros(preypopsize,7201);
preyusefulawake = zeros(preypopsize,gennum);
preyusefulsleep = zeros(preypopsize,gennum);
preyvecsum = zeros(preypopsize,gennum);
preyDm = zeros(preypopsize,7201);
preyDv = zeros(preypopsize,7201);
%preyH = zeros(preypopsize,7201);
preyawakestate = zeros(preypopsize,7201);
preyawakebouts = zeros(preypopsize,gennum);
preyawakedur = zeros(preypopsize,gennum);
preydiurnality = zeros(preypopsize,gennum);
preycob = zeros(preypopsize,gennum);
preymag = zeros(preypopsize,gennum);
preydistratio = zeros(preypopsize,gennum);
preydistx = zeros(preypopsize,gennum);
preyperiod = zeros(preypopsize,gennum);

predvecsum = zeros(predpopsize, gennum);
predconstate = zeros(predpopsize,7201);
predusefulawake = zeros(predpopsize,gennum);
predusefulsleep = zeros(predpopsize,gennum);
predDm = zeros(preypopsize,7201);
predDv = zeros(preypopsize,7201);

predawakestate = zeros(predpopsize,7201);
predawakebouts = zeros(predpopsize,gennum);
predawakedur = zeros(predpopsize,gennum);
preddiurnality = zeros(predpopsize,gennum);
predcob = zeros(predpopsize,gennum);
predmag = zeros(predpopsize,gennum);
preddistratio = zeros(predpopsize,gennum);
preddistx = zeros(predpopsize,gennum);
predperiod = zeros(predpopsize,gennum);

preytaufreqrec = zeros(numel(allele22),gennum); %(#unique alleles, frequency, #generations + 1)
predtaufreqrec = zeros(numel(allele22),gennum);
preychifreqrec = zeros(numel(allele8),gennum);
predchifreqrec = zeros(numel(allele8),gennum);
preyafreqrec = zeros(numel(allele3),gennum);
predafreqrec = zeros(numel(allele3),gennum);
preyD0freqrec = zeros(numel(allele10),gennum);
predD0freqrec = zeros(numel(allele10),gennum);

%% Genetic Algorithm


for j = 1:gennum
    generation = j
	  predpopgenome(:,:,j) = repmat(pars,1,predpopsize);
preypopgenome(:,:,j) = repmat(pars,1,preypopsize);
    % Mutations!
    if j == 1
  for m = 1:length(mutatingparams)
            preypopgenome(3,:,j) =allele3(randperm(numel(allele3)));
preypopgenome(8,:,j) =allele8(randperm(numel(allele8)));
preypopgenome(10,:,j) =allele10(randperm(numel(allele10)));
preypopgenome(22,:,j) =allele22(randperm(numel(allele22)));
            
predpopgenome(3,:,j) =allele3(randperm(numel(allele3)));
predpopgenome(8,:,j) =allele8(randperm(numel(allele8)));
predpopgenome(10,:,j) =allele10(randperm(numel(allele10)));
predpopgenome(22,:,j) =allele22(randperm(numel(allele22)));
            
        end
 else
   newallele3prey = [];
newallele8prey = [];
newallele10prey = [];
newallele22prey = [];
newallele3pred = [];
newallele8pred = [];
newallele10pred = [];
newallele22pred = [];
for k = 1:fittest/2
            
	  preyparents = [preypopgenome(:,2*k-1,j-1),preypopgenome(:,2*k,j-1)];
allele3indprey(1) = find(preyparents(3,1) == allele3);
allele3indprey(2) = find(preyparents(3,2) == allele3);
allele3kidsprey = allele3(randi([min(allele3indprey),max(allele3indprey)],[1,4]));
newallele3prey = [allele3kidsprey,newallele3prey];
            
            
allele8indprey(1) = find(preyparents(8,1) == allele8);
allele8indprey(2) = find(preyparents(8,2) == allele8);
allele8kidsprey = allele8(randi([min(allele8indprey),max(allele8indprey)],[1,4]));
newallele8prey = [allele8kidsprey,newallele8prey];
            
allele10indprey(1) = find(preyparents(10,1) == allele10);
allele10indprey(2) = find(preyparents(10,2) == allele10);
allele10kidsprey = allele10(randi([min(allele10indprey),max(allele10indprey)],[1,4]));
newallele10prey = [allele10kidsprey,newallele10prey];
            
allele22indprey(1) = find(preyparents(22,1) == allele22);
allele22indprey(2) = find(preyparents(22,2) == allele22);
allele22kidsprey = allele22(randi([min(allele22indprey),max(allele22indprey)],[1,4]));
newallele22prey = [allele22kidsprey,newallele22prey];
            
predparents = [predpopgenome(:,2*k-1,j-1),predpopgenome(:,2*k,j-1)];
            
allele3indpred(1) = find(predparents(3,1) == allele3);
allele3indpred(2) = find(predparents(3,2) == allele3);
allele3kidspred = allele3(randi([min(allele3indpred),max(allele3indpred)],[1,4]));
newallele3pred = [allele3kidspred,newallele3pred];
            
allele8indpred(1) = find(predparents(8,1) == allele8);
allele8indpred(2) = find(predparents(8,2) == allele8);
allele8kidspred = allele8(randi([min(allele8indpred),max(allele8indpred)],[1,4]));
newallele8pred = [allele8kidspred,newallele8pred];
            
allele10indpred(1) = find(predparents(10,1) == allele10);
allele10indpred(2) = find(predparents(10,2) == allele10);
allele10kidspred = allele10(randi([min(allele10indpred),max(allele10indpred)],[1,4]));
newallele10pred = [allele10kidspred,newallele10pred];
            
allele22indpred(1) = find(predparents(22,1) == allele22);
allele22indpred(2) = find(predparents(22,2) == allele22);
allele22kidspred = allele22(randi([min(allele22indpred),max(allele22indpred)],[1,4]));
newallele22pred = [allele22kidspred,newallele22pred];
            
            
        end
        preypopgenome(3,:,j) = newallele3prey;
preypopgenome(8,:,j) = newallele8prey;
preypopgenome(10,:,j) = newallele10prey;
preypopgenome(22,:,j) = newallele22prey;
predpopgenome(3,:,j) = newallele3pred;
predpopgenome(8,:,j) = newallele8pred;
predpopgenome(10,:,j) = newallele10pred;
predpopgenome(22,:,j) = newallele22pred;
        
        %random mutations
        pois = poissrnd(avmutationpergen);
        if pois==0
	  else
            for m = 1:pois
		      mutalleleIND = [mutatingparams(randi(numel(mutatingparams))), randi([1,preypopsize])];
mutallelevec = allele3*(mutalleleIND(1,1)==3)+allele8*(mutalleleIND(1,1)==8)+allele10*(mutalleleIND(1,1)==10)+allele22*(mutalleleIND(1,1)==22);
mutallele =  mutallelevec(randi(numel(mutallelevec)));
preypopgenome(mutalleleIND(1,1), mutalleleIND(1,2),j) = mutallele;
                
mutalleleIND = [mutatingparams(randi(numel(mutatingparams))), randi([1,predpopsize])];
mutallelevec = allele3*(mutalleleIND(1,1)==3)+allele8*(mutalleleIND(1,1)==8)+allele10*(mutalleleIND(1,1)==10)+allele22*(mutalleleIND(1,1)==22);
mutallele =  mutallelevec(randi(numel(mutallelevec)));
predpopgenome(mutalleleIND(1,1), mutalleleIND(1,2),j) = mutallele;
                
            end
        end
    end
    %Calculating metrics for ea member
    parfor i=1:predpopsize
        [preyconstate(i,:), preyusefulawake(i,j), preyusefulsleep(i,j), preyvecsum(i,j), preyDm(i,:), preyDv(i,:),~,preyawakebouts(i,j), preyawakedur(i,j),~, preycob(i,j), preymag(i,j), preydistratio(i,j), preydistx(i,j), preyperiod(i,j)]= prey_pred_hardswtichmodel_bufferincluded(preypopgenome(:,i,j),buffermins, 0);
        [predconstate(i,:), predusefulawake(i,j), predusefulsleep(i,j), predvecsum(i,j), predDm(i,:), predDv(i,:),~,predawakebouts(i,j), predawakedur(i,j),~, predcob(i,j), predmag(i,j), preddistratio(i,j), preddistx(i,j),predperiod(i,j)]= prey_pred_hardswtichmodel_bufferincluded(predpopgenome(:,i,j),buffermins, 0);
        
    end
    linth = 0:0.01:3*24; %3 days
    % Cost Function
    avpreyawakestate(j,:) = mean(preyconstate);
avpredawakestate(j,:) = mean(predconstate);
    
for i = 1:predpopsize
	
    
    preyfit(i,j) = (-max(0,mean(preyconstate(i,:).*avpreyawakestate(j,:))-a*mean(preyconstate(i,:).*avpredawakestate(j,:))))*(preyusefulsleep(i,j).^b*preyusefulawake(i,j));
predfit(i,j) = -sqrt(mean(predconstate(i,:).*avpredawakestate(j,:)))*(predusefulsleep(i,j).^c)*(mean(predconstate(i,:).*avpreyawakestate(j,:))); %this is to be minimised  
    end
    %Sorting prey and pred according to fitness
[preyfit(:,j), rankingprey] = sort(preyfit(:,j));
[predfit(:,j), rankingpred] = sort(predfit(:,j));
    %Sorting genomes from fittest to least fittest
    preypopgenome(:,:,j) = preypopgenome(:,rankingprey,j);
predpopgenome(:,:,j) = predpopgenome(:,rankingpred,j);

  rankedpreyconstate = preyconstate(rankingprey,:);
    rankedpredconstate = predconstate(rankingpred,:);
    
    fittestpreyconstate(j,:) = rankedpreyconstate(1,:);
    fittestpredconstate(j,:) = rankedpredconstate(1,:);
   %Calculating awake dur and awake bouts stats
   PREYawakedurmed(:,j) = median(preyawakedur(:,j));
PREDawakedurmed(:,j) = median(predawakedur(:,j));
        
PREYawakeduriqr(:,j) = iqr(preyawakedur(:,j));
PREDawakeduriqr(:,j) = iqr(predawakedur(:,j));
        
PREYawakeboutsmed(:,j) = median(preyawakebouts(:,j));
PREDawakeboutsmed(:,j) = median(predawakebouts(:,j));
    
PREYawakeboutsiqr(:,j) = iqr(preyawakebouts(:,j));
PREDawakeboutsiqr(:,j) = iqr(predawakebouts(:,j));

%Prey D0 gene, param 3
  PREYD0med(j) = median(preypopgenome(3,:,j)); 
PREYD0iqr(j) = iqr(preypopgenome(3,:,j));
%Pred D0 gene, param 3
  PREDD0med(j) = median(predpopgenome(3,:,j)); 
PREDD0iqr(j) = iqr(predpopgenome(3,:,j));
%Prey chi gene, param 8
  PREYCHImed(j) = median(preypopgenome(8,:,j));
PREYCHIiqr(j) = iqr(preypopgenome(8,:,j));
%Pred chi gene, param 8
  PREDCHImed(j) = median(predpopgenome(8,:,j));
PREDCHIiqr(j) = iqr(predpopgenome(8,:,j));
%Prey a gene, param 10
  PREYAmed(j) = median(preypopgenome(10,:,j));
PREYAiqr(j) = iqr(preypopgenome(10,:,j));
%Pred a gene, param 10
  PREDAmed(j) = median(predpopgenome(10,:,j));
PREDAiqr(j) = iqr(predpopgenome(10,:,j));
        %Prey tau param/gene
        PREYTAUmed(j) = median(preypopgenome(22,:,j));
PREYTAUiqr(j) = iqr(preypopgenome(22,:,j)); 
        %Pred tau param/gene
        PREDTAUmed(j) = median(predpopgenome(22,:,j));
PREDTAUiqr(j) = iqr(predpopgenome(22,:,j)) ;

end


% FOR EVERY n (parameter set), AN EXCEL FILE IS CREATED. 
  filename = strcat(num2str(n),'_paramsetdata.xlsx')
  varnames = { 'Median_prey_awake_dur', 'IQR_prey_awake_dur', 'Median_pred_awake_dur', 'IQR_pred_awake_dur', ...
	       'Median_prey_awake_bouts', 'IQR_prey_awake_bouts', 'Median_pred_awake_bouts', 'IQR_pred_awake_bouts', ...
	       'PreyD0med','PreyD0iqr','PredD0med','PredD0iqr','PreyCHImed','PreyCHIiqr','PredCHImed','PredCHIiqr', ...
	       'PreyAmed','PreyAiqr','PredAmed','PredAiqr','PreyTAUmed','PreyTAUiqr','PredTAUmed','PredTAUiqr'};
    
T1 = table( PREYawakedurmed',  PREYawakeduriqr',PREDawakedurmed',  PREDawakeduriqr',...
	    PREYawakeboutsmed',  PREYawakeboutsiqr',PREDawakeboutsmed',  PREDawakeboutsiqr',...
	    PREYD0med',PREYD0iqr',PREDD0med',PREDD0iqr', PREYCHImed',PREYCHIiqr',PREDCHImed',PREDCHIiqr',...
	    PREYAmed',PREYAiqr',PREDAmed',PREDAiqr', PREYTAUmed',PREYTAUiqr',PREDTAUmed',PREDTAUiqr',...
	    'VariableNames',varnames);
writetable(T1,filename)
%% awakestates
filename1 = strcat(num2str(n),'_avpreyavawakestate.xlsx');
allavpreyawakestate = table(avpreyawakestate);
writetable(allavpreyawakestate,filename1);

filename1 = strcat(num2str(n),'_avpredavawakestate.xlsx');
allavpredawakestate = table(avpredawakestate);
writetable(allavpredawakestate,filename1);

filename3 = strcat(num2str(n),'_fittestpreyawakestate.xlsx');
T3 = table(fittestpreyconstate);     
writetable(T3,filename3);

filename4 = strcat(num2str(n),'_fittestpredawakestate.xlsx');
T4 = table(fittestpredconstate);    
writetable(T4,filename4);
ended = datestr(now)
toc
end
end

%% All functions used in evol. 
  function [constate, usefulawakedur, usefulsleepdur, vecsum, Dm, Dv, awakestate,awakebouts, awakedur, diurnality, cob, mag, distratio, distx, period] = prey_pred_hardswtichmodel_bufferincluded(pars,buffermins, plotting)
  %function []=hardswtichmodel(pars,plotting)
  %% Code for simulating the sleep/wake model in Phillips et al. (2013)

%% Settings for simulation
% Conversion factors
  hpd = 24; % hours per day
  mph=60; % minutes per hour


  %% Run settings (Duration, cont, merge)
    days = 3;
dur    = days*24;            % run duration (hours)
  dur = dur*mph;

%% Initial Conditions (Declare V and C)
% Initialise everything
  V = zeros(1,4);

% Set up the base initial conditions

V(1) = 16;          % H (excluded later if ramp)
V(2) = -.7;           % x
V(3) = 0;           % xc
V(4) = .03;        % n

%% Solve the differential equations

p=pars;

%    Code to limit transients
[t,V]=ode23s(@hardswitchde,[0;24*10*mph],V',odeset,p);

[t,V]=ode23s(@hardswitchde,[0;dur],V(end,:)',odeset,p);
% Define the model variables
H = V(:,1); %homeostatic sleep drive (concentration of substance)
x = V(:,2); %SCN activity (dimensionless)
xc = V(:,3); %y (complementary variable) (dimensionless)
n = V(:,4); %fraction of activated photoreceptors (dimensionless)

th = t/60;% Time in hours

%% Interping
linth = min(th):0.01:max(th);
H = interp1(th,H,linth);
x = interp1(th,x,linth);
xc = interp1(th,xc,linth);
n = interp1(th,n,linth);
%% awake state
awakestate = sign(diff(H)); %+1 means awake, -1 means asleep
awakestate(end+1)=awakestate(end);
constate = awakestate;
awakebouts = sum(diff(awakestate)>0)/days; %awakebouts counted for each time model wakes

%% useful sleep and useful awake time for fitness function
sleep2wake = find(diff(constate)==2);
wake2sleep = find(diff(constate)==-2);
buffertimepnts = floor(100*buffermins/mph);
%first 10 minutes of wake is 0
for i = 1:length(sleep2wake)
    if sleep2wake(i)>=(7201-buffertimepnts)
    else
        constate(sleep2wake(i)+1: sleep2wake(i)+buffertimepnts) = zeros(1,buffertimepnts);
    end
end
%first 10 minutes of sleep is 0
for j = 1:length(wake2sleep)
    if wake2sleep(j)>=(7201-buffertimepnts)
    else
        constate(wake2sleep(j)+1:wake2sleep(j)+buffertimepnts) = zeros(1,buffertimepnts);
    end
end

constate=0.5*(constate+1);
%%
%%Calculating Dm and Dv
Dm = pars(6) + pars(9)*(pars(10)*pars(11)*(x + pars(12)) + pars(13));
Dv=pars(4)*H + pars(3) + pars(2)*(pars(10)*pars(11)*(x + pars(12)) + pars(13));


%% Dist ratios
distx = mean (abs(   x(intersect(find(linth>=0*hpd),find(linth<=1*hpd))) - x(intersect(find(linth>=hpd*2),find(linth<=3*hpd)))  ));

distx12 = mean (abs(   x(intersect(find(linth>=0*hpd),find(linth<=1*hpd))) - x(intersect(find(linth>=hpd*1),find(linth<=2*hpd)))  ));

distx23 = mean (abs(   x(intersect(find(linth>=1*hpd),find(linth<=2*hpd))) - x(intersect(find(linth>=hpd*2),find(linth<=3*hpd)))  ));

distratio=distx12/distx23;


%% Calculating duration model is awake
awakedur = length(find(awakestate==1))/length(awakestate); %1 means awake whole time, 0 means sleep whole time

usefulawakedur = length(find(constate==1))/length(constate);
usefulsleepdur = length(find(constate==0))/length(constate);
%% Calculating diurnality
I = elight(linth*60)/1000;
diurnality = (sum(awakestate.*I)+0.00001)/(sum(awakestate) + 0.00002)  ;

%% Average expressed period
phi = unwrap(-atan2(xc,x));

cycles = (phi(end)-phi(1))/(2*pi);
period = (linth(end)-linth(1))/cycles; %Average expressed period (how well model entrains)


%% Average wake clock time

activeangle = mod(linth, 24)*pi/12;
vecsum = sum(awakestate.*exp(sqrt(-1)*activeangle));
cob = 12/pi*angle(vecsum);
mag = abs(vecsum)/length(linth);

end

function dVdt = hardswitchde(t,V,pars)
persistent state
if isempty(state)
    state=1;
end

hpd = 24; % hours per day
mph=60;
kappainv = 2*pi/(hpd*mph);       % 1/kappa (min^{-1})

%% Determining drive values
%[gamma, B, ] = ldrive(V(2),V(3),V(4),elight(t)*(1-0.97*(state==1)),pars);
dawn = 6;
dusk = 18;
Imax = 1000;
I = (Imax*(mod(t/60,24)>dawn).*(mod(t/60,24)<=dusk))*(state==1);
%alpha = alpha_0.* ((I/I_0).^pars(27)) .* I ./ (I+I_1);
gamma = pars(25)* I^(pars(27)+1)./(I+pars(28));

%Bhat: activation of photoreceptors results in a photic drive to the pacemaker
%that is assumed to be proportional to the rate of activation. G is
%proportionality constant.
Bhat = (pars(26).* (1-V(4)) .* gamma);

%Response of circadian pacemaker to light.
B = Bhat .* (1-pars(29).*V(2)) .* (1-pars(29).*V(3));

Dm = pars(6) + pars(9)*(pars(10)*pars(11)*(V(2) + pars(12)) + pars(13));
Dv=(pars(4)*V(1) + pars(3) + pars(2)*(pars(10)*pars(11)*(V(2) + pars(12)) + pars(13)) + pars(14)*Bhat);

highDv = 0.4758 + 0.8526*Dm + 0.4720*Dm^2 + 0.0394*Dm^3;
lowDv = 0.5520 + 0.8430*Dm - 0.1342*Dm^2 + 0.0113*Dm^3;

cstate = -1*(Dv>highDv) + 1*(lowDv>Dv) ;

if cstate==0;
 else
   state = cstate;
end


dVdt = [(pars(7)*4.8*60*(state==1) - V(1))/pars(8); % H
					  kappainv*(V(3) + B +  pars(15)*(pars(16)*V(2) + pars(17)*V(2)^3 + pars(18)*V(2)^5) + pars(19)*V(2)^7); % x
					  kappainv*(B*(pars(20)*V(3) - pars(21)*V(2)) - V(2)*pars(22));  % xc
					  (gamma*(1-V(4)) - pars(23)*V(4))]; % n

end
% Function for describing environmental light patterns
%t is time of day in minutes
function I = elight(t)

  dawn = 6;
dusk = 18;
Imax = 1000;



I = Imax*(mod(t/60,24)>dawn).*(mod(t/60,24)<=dusk);

end
