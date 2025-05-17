%% V1.2: KaiABC simulation with binding correlation, period estimation, and TTFL

% Description: Simulates clock fidelity based off bulk experimental data &
% parameters. This version incoporates TTFL into the model and modifies the KaiB
% distribution methology.

% Info: Exports Model_Results.mat with final results and creates figure.

%% Begin Input Parameters %%
clear

%% Input Parameters: Encapsulation (Gamma distrbution)
corr_ABC    = 0.13;          % (For KaiAC only) Default: 0.13 (exp) %Correlation ratio for KaiA binding to KaiC (monomer ratio)
CV_avg      = 0.25;          % Default: 0.25 %Average CV value from literature
conc_stioch = [1.2*0.9,3.6*1.2,2.64*1.1]; % KaiA : KaiB : KaiC concentration stiochiometry (µM) @ 1X % Default saved: [1.2*0.9,3.6*1.2,2.64*1.1];
conc_rel    = [0.5 ;0.75; 1.0; 1.5; 2.5]; % List of Concentrations to run

SasA_CikA_Supp = 1;         %1 for SasA/CikA support, 0 for off

ratio_ABC = [1,1,1];        % Ratio of KaiA:KaiB:KaiC binding during loading (set to 1:1:1)

%% Input Parameters: Clock reaction simulation
sim_osc_all = [];           % Set Output Variable Name

diam = [2];                 % List of vesicle diameters to run
radi = diam/2;              % Calculate Radius
n    = 5000;                % # of simulated vesicles (default = 5000)

crit_conc_rel = [0.5,0.5,0.5]; % Default: 0.5 %critical concentration (relative) (from bulk experiments) for KaiA & KaiB & KaiC
crit_ratio    = 0.5;        % Default: 0.5 %Critical ratio relative 
max_ratio     = 3.0;        % Default: 3.0 %Max relative ratios
b             = 650;        % Default: 650 %KaiB Molecules bound per µm^2 surface area 
KaiB_free_per = 0.5;        % percent bound KaiB to membrane

title1   = strcat('b = ', num2str(b)); %, ',  BCcobind = ', num2str(kaic_cobind));
plotfigs = 1;            % Plot simulation fidelity data?

vol_um3 = (4/3)*pi*(radi).^3; % Calculate volume of vesicle (µm^3)
vol_L   = vol_um3 * 1e-15;    % Volume in L
SA      = 4*pi*(radi).^2;     % Surface area in µm^2
N_A     = 6.022 * 10^23;      % Avogadro's number

KaiC_mem_assoc = 0.0;

%% Inital Calculations and Parameter Loading %%
%% Critical Concentration  & Ratio calculation
% Calculates critical concentration in µM

if SasA_CikA_Supp == 0
    % Critical Concentration
    min_kaiA = 1.2*crit_conc_rel(1);   % KaiA limiting concentration µM % Set to 0.3 for SasA/CikA support
    min_kaiB = 3.5*crit_conc_rel(2);   % KaiB limiting concentration µM % Set to 0.9 for SasA/CikA support
    min_kaiC = 3.5*crit_conc_rel(3);   % KaiC limiting concentration µM % Set to 1.75 for SasA/CikA support
    
    % Limiting ratios
    min_R_kaiAC = crit_ratio*(1.2/3.5); % Limiting KaiA:KaiC ratio %Set to 0.09 for SasA/CikA support
    min_R_kaiBC = crit_ratio;           % Limiting KaiB:KaiC ratio %Set to 0.25 for SasA/CikA support

    % Maximum ratios
    max_R_kaiAC = max_ratio*(1.2/3.5);  % Max KaiA:KaiC ratio: 1.03

else 
    min_kaiA = 0.3;                    % KaiA limiting concentration µM % Set to 0.3 for SasA/CikA support
    min_kaiB = 0.9;                    % KaiB limiting concentration µM % Set to 0.9 for SasA/CikA support
    min_kaiC = 1.75;                   % KaiC limiting concentration µM % Set to 1.75 for SasA/CikA support
    
        % Limiting ratios
    min_R_kaiAC = 0.09;                % Limiting KaiA:KaiC ratio %Set to 0.09 for SasA/CikA support
    min_R_kaiBC = 0.25;                % Limiting KaiB:KaiC ratio %Set to 0.25 for SasA/CikA support

    % Maximum ratios
    max_R_kaiAC = max_ratio*(1.2/3.5);  % Max KaiA:KaiC ratio: 1.03

end

%% Other calculations


%% Period Simulation Setup using period offset lookup tables
%Loads period lookup tables:
load('Period Simulation\Period_lookup_tables.mat') 

%Loads experimental period data for comparison:
load('Period Simulation\Exp_periods.mat') 

% Loads number of oscilating vesicles by size:
load('Period Simulation\N_oscillating_vesbysize.mat')


%% Amplitude Simulation Setup - Based on bulk data and KaiB concentration

load('Amplitude Simulation\Amp_Lookup.mat')
amp_bulk_calcurve = @(x) (0.074/3.5)*x + 0.09; 

%% Initalize variables
sim_period_all = [];
sim_period_SD = [];
sim_osc_all = [];
adj_amp_mean = [];
adj_amp_SD = [];
sim_osc_all = [];


%% Begin Simulation Calculations %%


%% Simulate encapsulated protein concentrations using gamma distribution

conc_relstring = 'c'+strrep(string(num2str(conc_rel,'%.2f')),'.','_'); %Concentration converted into strings

KaiC1x = conc_stioch(3)*(1/6); % The multiplers signify hexamer (KaiC), tetramer (KaiB), and dimer (KaiA).
KaiB1x = conc_stioch(2)*(1/4); 
KaiA1x = conc_stioch(1)*(1/2); 

% Initatize variables
ABC_mean  = []; C_mean  = []; B_mean  = []; A_mean  = [];
ABC_gam   = []; C_gam   = []; B_gam   = []; A_gam   = [];
ABC_k     = []; C_k     = []; B_k     = []; A_k     = [];
ABC_theta = []; C_theta = []; B_theta = []; A_theta = [];

A_mean_list = [];

%Covert to Approriate Unit for SI unit liters (liter = dm^3)
b_2 = b * 10^(5)^(2); %from molecules/um^2 to molecules/dm^2 | (L = dm^3)& (um = 10^-5 dm)
radi_2 = radi * 10^(-5); %from um to dm (um = 10^-5) 

% Loop to run all selected concentrations
for i = 1:length(conc_rel)
    % Note: Mean KaiA & KaiB concentrations is reduced because of binding to KaiC
    % Equation for mean from gamma parameters: mean = k (shape) * theta (scale)

    %% KaiABC Encapsulation
    ABC_mean(i) = KaiC1x*conc_rel(i)*corr_ABC*6; % Mean # of KaiABC formed in monomers
    ABC_k(i)    = 1/CV_avg^2;
    ABC_theta(i)= ABC_mean(i)*CV_avg^2;

    %% KaiC Encapsulation
    C_mean(i)   = KaiC1x*conc_rel(i)- (KaiC1x*conc_rel(i)*KaiC_mem_assoc) - (ABC_mean(i)*ratio_ABC(3)*(1/6)); %Mean reduced due to forming KaiABC (by correlation factor)
    C_k(i)      = 1/CV_avg^2;
    C_theta(i)  = C_mean(i)*CV_avg^2;
        
    %% KaiB Encapsulation
    B_mean(i)   = KaiB1x*conc_rel(i);
    B_k(i)      = 1/CV_avg^2;
    B_theta(i)  = B_mean(i)*CV_avg^2;

    
    %% KaiA Encapsulation
    A_mean(i)   = KaiA1x*conc_rel(i)-(ABC_mean(i)*ratio_ABC(1)*(1/2)); %Mean reduced due to forming KaiABC (by correlation factor)
    A_k(i) = 1/CV_avg^2;
    A_theta(i) = A_mean(i)*CV_avg^2;
    A_mean_list(i) = A_mean(i);
    
    %% Random Assignment
    for j = 1:length(diam)
        ABC_gam{i}(:,j)  = gamrnd(ABC_k(i),ABC_theta(i),5000,1);
        C_gam{i}(:,j)    = gamrnd(C_k(i),C_theta(i),5000,1);
        C_gam{i}(:,j)    = (C_gam{i}(:,j) + (KaiC1x*conc_rel(i)*KaiC_mem_assoc))*6 + (ABC_gam{i}(:,j).*ratio_ABC(3));
        B_gam{i}(:,j)    = gamrnd(B_k(i),B_theta(i),5000,1)*4;

        % Check if subtraction reduces mean KaiA below 0
        if A_mean(i) > 0
            A_gam{i}(:,j) = gamrnd(A_k(i),A_theta(i),5000,1);
            A_gam{i}(:,j) = A_gam{i}(:,j)*2 + (ABC_gam{i}(:,j)*ratio_ABC(1));
        else
            A_mean(i)    = 0;
            A_gam{i}(:,j)     = zeros(5000,1); %If all KaiA is bound to KaiABC complexes there is no KaiA encapsulation
            limA_mean    = KaiA1x*conc_rel(i);
            % Special KaiABC distribution if KaiA is limited (so no extra KaiA appears)
            A_gam{i}(:,j)     = A_gam{i}(:,j)*2 + ABC_gam{i}(:,j)*limA_mean/ABC_mean(i); 
        end
    end
    % Validation generated means
    C_gam_mean = cellfun(@mean,C_gam(i),'UniformOutput',false);
    B_gam_mean = cellfun(@mean,B_gam(i),'UniformOutput',false);
    A_gam_mean = cellfun(@mean,A_gam(i),'UniformOutput',false);

end

KaiC = table(C_gam{:},'VariableNames',conc_relstring);
KaiB = table(B_gam{:},'VariableNames',conc_relstring);
KaiA = table(A_gam{:},'VariableNames',conc_relstring);


%% Simulate Clock Vesicle Fidelity (How many oscillate)

%% Combined encapsulation in vesicles.

for conc = 1:length(conc_rel)
    %% Load data and determine whether sim vesicles oscillate.
    % Note: Column 1 = KaiA Conc., Column 2 = KaiB Conc., Column 3 = KaiC Conc.(µM)

    i_KaiA{conc} = KaiA{:,conc};%conc_KaiABC_sim{conc}(:,1);
    i_KaiB{conc} = KaiB{:,conc};%conc_KaiABC_sim{conc}(:,2);
    i_KaiC{conc} = KaiC{:,conc};%conc_KaiABC_sim{conc}(:,3);

    %% Remove membrane associated proteins from free protein concentration

    KaiB_free{conc} = i_KaiB{conc} .* ones(1,length(vol_L)) .* KaiB_free_per;
    KaiA_free{conc} = i_KaiA{conc} .* ones(1,length(vol_L)); % no binding, multipled by arrays of 1s to create column for each size
    KaiC_free{conc} = i_KaiC{conc} .* ones(1,length(vol_L)); % no binding, multipled by arrays of 1s to create column for each size
   

    %% Check and replace negative values

    KaiB_free{conc}(KaiB_free{conc}<0) = 0;

    %% Calculate KaiA:KaiC and KaiB:KaiC ratios

    KaiAC_ratio{conc} = (KaiA_free{conc} ./ KaiC_free{conc});
    KaiBC_ratio{conc} = (KaiB_free{conc} ./ KaiC_free{conc});

    %% Check if conditions are met, and determining oscillating (osc) vesicles

    checkC_kaiA{conc} = KaiA_free{conc} >= min_kaiA;
    checkC_kaiB{conc} = KaiB_free{conc} >= min_kaiB;
    checkC_kaiC{conc} = KaiC_free{conc}   >= min_kaiC;
    checkR_kaiAC{conc} = KaiAC_ratio{conc} >= min_R_kaiAC;
    checkR_kaiBC{conc} = KaiBC_ratio{conc} >= min_R_kaiBC;
    checkR_maxkaiAC{conc} = KaiAC_ratio{conc} <= max_R_kaiAC;
    osc{conc} = checkC_kaiA{conc} & checkC_kaiB{conc} & checkC_kaiC{conc} & checkR_kaiAC{conc} & checkR_kaiBC{conc} & checkR_maxkaiAC{conc};


    %% Period Simulation using period offset lookup tables

    mean_period = interp1(period_bulk.conc,period_bulk.period,KaiC_free{conc}/3.5,'linear','extrap');
    
    period_conc{conc} = interp1(period_KaiAC.ratio,period_KaiAC.p_offset,KaiAC_ratio{conc},'linear','extrap')...
                    + interp1(period_KaiBC.ratio,period_KaiBC.p_offset,KaiBC_ratio{conc}, 'linear','extrap')+ mean_period(conc); %last number is mean period

    %Remove all periods for non-oscillating vesicles and replace with a NaN
    % Columns = size, Rows = individual vesicles
    period_conc{conc}(osc{conc} == 0) = NaN;
    period_conc_mat(:,conc) = period_conc{conc}(:);
    
    period_nonan{conc} = period_conc{conc}(~isnan(period_conc{conc}));
    
    % Note that (rows = concentration , columns = size)
    period_mean{conc} = mean(period_conc{conc},'omitnan');
    period_SD{conc} =   std(period_conc{conc},'omitnan');
    sim_period_all = vertcat(sim_period_all, period_mean{conc});
    sim_period_SD = vertcat(sim_period_SD, period_SD{conc});

    %% Amplitude Simulation - Based on bulk data and KaiB concentration

    % Note that (rows = concentration , columns = size)
    KaiC_free_osc{conc} = KaiC_free{conc};
    KaiC_free_osc{conc}(osc{conc}==0) = NaN;
    KaiB_free_osc{conc} = KaiB_free{conc};
    KaiB_free_osc{conc}(osc{conc}==0) = NaN;
    KaiA_free_osc{conc} = KaiA_free{conc};
    KaiA_free_osc{conc}(osc{conc}==0) = NaN;
    
    inital_amp{conc} = interp1(amp_base.conc,amp_base.amp,KaiC_free_osc{conc}./3.5,'linear','extrap');
    
    KaiCB_ratio{conc} = (KaiC_free{conc} ./ KaiB_free{conc});
    KaiCB_ratio_osc{conc} = KaiCB_ratio{conc};
    KaiCB_ratio_osc{conc}(osc{conc}==0) = NaN;
    
    KaiAC_ratio{conc} = (KaiA_free{conc} ./ KaiC_free{conc});
    KaiAC_ratio_osc{conc} = KaiAC_ratio{conc};
    KaiAC_ratio_osc{conc}(osc{conc}==0) = NaN;
    
    amp_adj{conc} = (interp1(amp_KaiCB.kaicb_ratio,amp_KaiCB.amp_delta,KaiCB_ratio_osc{conc},'linear','extrap')...
            + interp1(amp_KaiAC.kaiac_ratio,amp_KaiAC.amp_delta,KaiAC_ratio_osc{conc},'linear','extrap')+ inital_amp{conc});
    
    amp_min_check = amp_adj{conc}<0;
    amp_adj{conc}(amp_min_check) = 0;
    
    
    amp_nonan{conc} = amp_adj{conc}(~isnan(amp_adj{conc}));
    
    amp_adj_mat(:,conc)= amp_adj{conc}(:);

    %% Fraction oscillating 2
    % Fraction oscillating (rows = concentration , columns = size)
    sim_osc_all = vertcat(sim_osc_all,sum(osc{conc})/n);

end

%% Figure Plotting

% Concentration values
sim_conc = conc_rel; %Relative values (0.5x, 0.75x, etc...)
sim_size = diam;

% Plot Figure with simulation results
if plotfigs == 1
    h1 = figure;
    hold on; plot(sim_conc, sim_osc_all,'o-','LineWidth',3,'MarkerSize',10)
    colororder(flip(['#123624';'#00634F';'#178D97';'#64A7CE';'#A1BAD9';'#D0D1E6';'#123624';'#00634F';'#178D97';'#64A7CE';'#A1BAD9';'#D0D1E6']));
    set(gca,'fontsize',18,'LineWidth',3); grid off; box on;
    xticks([0.5:0.5:2.5])
    title(strcat('C:',num2str(crit_conc_rel(1)),', R:',num2str(crit_ratio),', b=', num2str(b),', Corr-BC:',num2str(corr_ABC),'/AC:',num2str(corr_ABC,2)));
    xlabel('Concentration (×)');
    ylabel('Fidelity');
    xlim([0.4 2.6])
    set(gca,'XMinorTick','on','YMinorTick','on')
    hA = gca; hA.XAxis.MinorTickValues = 0:0.25:3;
    legend(append(string(diam),' µm'))
end

%% Save Variables to File
conc_name = append(string(sim_conc'),'x');
diam_list = append(string(diam),' µm');

fidelity_table = array2table(sim_osc_all','VariableNames',conc_name,'RowNames',diam_list);
period_table = array2table(period_conc_mat,'VariableNames',conc_name); 
period_table_osc = cell2table(amp_nonan,'VariableNames',conc_name);
amplitude_table = array2table(amp_adj_mat,'VariableNames',conc_name);
amplitude_table_osc = cell2table(period_nonan,'VariableNames',conc_name); 

a=pwd;
outputFileNameMAT1 = 'Model_Results.mat';
save(outputFileNameMAT1,'fidelity_table','fidelity_table','period_table','period_table_osc','amplitude_table','amplitude_table_osc'); 