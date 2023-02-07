% Initial settings
close all;
clear;
clc;

%%% Initial settings: Set your path to dynasim and working directory.
% dynasim_path = '';                    % Add your DynaSim path here¡
% addpath(genpath(dynasim_path));
% root_path = '';                       % Add your working directory path here¡
% addpath(genpath(root_path));

%%% MODEL DEFINITION %%%
npop = 4;                            % Number of different populations: D1 and D2
nSPNs = 150;                         % Number of SPN cells in the pool
multiplicity = ones(1,npop)*nSPNs;   % Create a neuron matrix

%%% Define equations of cell model
eqns = {
  'dV/dt = @current'
  'V(0) = -65+5*randn(1,Npop)'       % mV, Initial conditions
  'monitor @current'
};

%%% Mechanisim parameters
g_l_D1 = 0.075;                      % spn_leak input current
g_l_D2 = 0.09;

g_cat_D1 = 0.018;                    % spn_Cat input current
g_cat_D2 = 0.025;

%%% Background Async mean input
g_extInp = 0.00035;                  % in mS/cm^2 milisiemens
baseline_extInp = 30e3;              % in spks/s

%%% Pfc input to SPN (excitatory beta input)
g_pfcInp = g_extInp;
test = 1.25;                         % new variable to increase the magnitude of the PFC input

DC_spont = 0;
DC = 44e3*test;
AC_spont = 0;
AC =8e3*test;

nofreq = 0;
freq_alphaInp = 10;                 % Alpha
freq_betaInp = 25;                  % High beta

latency_Inp = 10;                   % in ms
dutyPercentage_pfcInp = 50;
thrDUTY = cos(pi*dutyPercentage_pfcInp/100);

%%% related time variables
transient = 1000;
stim_duration = 1000;

%%% PFC input types: uncomment to simulute the different types of PFC inputs in case that you wanna see how the modeled neural circuit will respond to them.
%pfc_type = 'Spontaneous';
%pfc_type = 'DC';
%pfc_type = 'Alpha';
%pfc_type = 'Beta';

%%% Trial Types: uncomment the type of trial that you will simulate.
%pfc_type = 'Color';
pfc_type = 'Orientation';

if strcmp(pfc_type, 'Spontaneous')
    % Simulation time definition
    n_seq_stim = 1;
    tOff = transient + n_seq_stim*stim_duration;

    % Stimulus timing
    tOn_Inp = transient;
    tOff_Inp = transient + stim_duration;
    tOn_Inp_O = tOn_Inp;
    tOn_Inp_C = tOn_Inp;
    tOff_Inp_O = tOff_Inp;
    tOff_Inp_C = tOff_Inp;

%%%%%%%%%%%% %%% Spont %%% %%%%%%%%%%%%%%%%
    DC_Inp = DC_spont;

    AC_ref_O = AC_spont;
    AC_ref_C = AC_spont;

    freq_ref_C = nofreq;
    freq_ref_O = nofreq;

elseif strcmp(pfc_type, 'DC')
    % Simulation time definition
    n_seq_stim = 1;
    tOff = transient + n_seq_stim*stim_duration;

    % Stimulus timing
    tOn_Inp = transient;
    tOff_Inp = transient + stim_duration;
    tOn_Inp_O = tOn_Inp;
    tOn_Inp_C = tOn_Inp;
    tOff_Inp_O = tOff_Inp;
    tOff_Inp_C = tOff_Inp;

%%%%%%%%%%%%% %%% DC %%% %%%%%%%%%%%%%%%%
    DC_Inp = DC;

    AC_ref_O = AC_spont;
    AC_ref_C = AC_spont;

    freq_ref_C = nofreq;
    freq_ref_O = nofreq;

elseif strcmp(pfc_type, 'Alpha')
    % Simulation time definition
    n_seq_stim = 1;
    tOff = transient + n_seq_stim*stim_duration;

    % Stimulus timing
    tOn_Inp = transient;
    tOff_Inp = transient + stim_duration;
    tOn_Inp_O = tOn_Inp;
    tOn_Inp_C = tOn_Inp;
    tOff_Inp_O = tOff_Inp;
    tOff_Inp_C = tOff_Inp;

%%%%%%%%%%%%% %%% AC alpha%%% %%%%%%%%%%%%%%%
    DC_Inp = DC;

    AC_ref_O = AC;
    AC_ref_C = AC;

    freq_ref_C = freq_alphaInp;
    freq_ref_O = freq_alphaInp;

elseif strcmp(pfc_type, 'Beta')
    % Simulation time definition
    n_seq_stim = 1;
    tOff = transient + n_seq_stim*stim_duration;

    % Stimulus timing
    tOn_Inp = transient;
    tOff_Inp = transient + stim_duration;
    tOn_Inp_O = tOn_Inp;
    tOn_Inp_C = tOn_Inp;
    tOff_Inp_O = tOff_Inp;
    tOff_Inp_C = tOff_Inp;

%%%%%%%%%%%%% %%% AC beta%%% %%%%%%%%%%%%%%%
    DC_Inp = DC;

    AC_ref_O = AC;
    AC_ref_C = AC;

    % %%% AC Beta %%%
    freq_ref_C = freq_betaInp;
    freq_ref_O = freq_betaInp;

elseif strcmp(pfc_type, 'Orientation')
    % Simulation time definition
    n_seq_stim = 2;
    tOff = transient + n_seq_stim*stim_duration;

    % Stimulus timing
    tOn_Inp = transient + stim_duration;    % ms
    tOff_Inp = transient + n_seq_stim*stim_duration;  % ms
    tOn_Inp_O = tOn_Inp;
    tOn_Inp_C = tOn_Inp;
    tOff_Inp_O = tOff_Inp;
    tOff_Inp_C = tOff_Inp;

%%%%%%%%%%%%% %%% PFC input %%% %%%%%%%%%%%%%%%
    DC_Inp = DC;

    % Orientation trial specification: 'No Alpha's effect'.
    AC_ref_O = AC;
    colorReducedAmplFactor = 1.6; % high beta amplitude of color wrt orientation
    AC_ref_C= AC/colorReducedAmplFactor; % this factor only applies in colour populations during orientation trials.
    freq_ref_C = freq_betaInp;
    freq_ref_O = freq_betaInp;
else
    % Simulation time definition
    n_seq_stim = 2;
    tOff = transient + n_seq_stim*stim_duration;

%%%%%%%%%%%%% %%% PFC input %%% %%%%%%%%%%%%%%%
    DC_Inp = DC;

    AC_pfc_alphaInp = AC;
    colorTrialReducedAmplFactor = 1.3;
    AC_pfc_betaInp = AC/colorTrialReducedAmplFactor;  % High beta equal amplitude reduction factor for both kind of pops when color trial it's on going.

    % Orientation population
    % Stimulus timing
    tOn_alphaInp = transient;
    tOff_alphaInp = transient + stim_duration;   % ms
    tOn_betaInp = transient + stim_duration;     % ms
    tOff_betaInp = transient + n_seq_stim*stim_duration;  % ms

    tOn_Inp_O= [tOn_alphaInp; tOn_betaInp];
    tOff_Inp_O = [tOff_alphaInp; tOff_betaInp];

    freq_ref_O = [freq_alphaInp; freq_betaInp];
    AC_ref_O = [AC_pfc_alphaInp; AC_pfc_betaInp];

    % Color population
    % Stimulus timing
    tOn_Inp_C = transient + stim_duration;    % ms
    tOff_Inp_C = transient + n_seq_stim*stim_duration;  % ms
    freq_ref_C = freq_betaInp;
    AC_ref_C = AC_pfc_betaInp;
end

% Simulation time definition
tOn = 0;
tspan = [tOn,tOff];

%% Connectivity parameters
conn_prob_gaba_D1toD1 = 0.26;
conn_prob_gaba_D1toD2 = 0.06;
conn_prob_gaba_D2toD1 = 0.27;
conn_prob_gaba_D2toD2 = 0.36;

% Calculate initial Gaba conductance as full model
cross_inhibition_factor = 1.5;
connec_Fac_long = 0.55;
connec_Fac_short = 0.45;

g_gaba_D1_lowDA = 0.5/multiplicity(1);
g_gaba_D2_lowDA = 1.65*g_gaba_D1_lowDA*multiplicity(1)/multiplicity(2);
f_gaba_D1_highDA = 1.3;
f_gaba_D2_highDA = 0.5;

g_gaba_D1toD1 = g_gaba_D1_lowDA*f_gaba_D1_highDA;
g_gaba_D1toD2 = g_gaba_D1_lowDA*f_gaba_D1_highDA*cross_inhibition_factor;
g_gaba_D2toD1 = g_gaba_D2_lowDA*f_gaba_D2_highDA*cross_inhibition_factor;
g_gaba_D2toD2 = g_gaba_D2_lowDA*f_gaba_D2_highDA;

g_s_D1toD1 = connec_Fac_short*g_gaba_D1toD1;
g_s_D1toD2 = connec_Fac_short*g_gaba_D1toD2;
g_s_D2toD1 = connec_Fac_short*g_gaba_D2toD1;
g_s_D2toD2 = connec_Fac_short*g_gaba_D2toD2;

g_ns_D1toD1 = connec_Fac_long*g_gaba_D1toD1;
g_ns_D1toD2 = connec_Fac_long*g_gaba_D1toD2;
g_ns_D2toD1 = connec_Fac_long*g_gaba_D2toD1;
g_ns_D2toD2 = connec_Fac_long*g_gaba_D2toD2;

Tnames = {'D1toD1';'D1toD2';'D2toD1';'D2toD2'};
g_S = [g_s_D1toD1; g_s_D1toD2; g_s_D2toD1; g_s_D2toD2];
g_nS = [g_ns_D1toD1; g_ns_D1toD2; g_ns_D2toD1; g_ns_D2toD2];
Tconnec = table(g_S,g_nS,'RowName',Tnames);

% Calculate Tau Gaba for each connection
median_tau_gaba = 30.4;      % ms
sig_tau_gaba = 8.2;          % ms
tau_gaba_D1 = median_tau_gaba + sig_tau_gaba*randn(1,multiplicity(1));
tau_gaba_D2 = median_tau_gaba + sig_tau_gaba*randn(1,multiplicity(2));

% Tau STD
tauD_D1_highDA = 1030;    % from 275 % ms
tauD_D2_highDA = 210;     % from 415

deltaD_gabaD1 = 0.35;
deltaD_gabaD2 = 0.35;

tau_riseD = 10; % 200; % 500;ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model definition and simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% %%% Population definitions %%% %%%%%%%%%%%%%
s = [];

s.populations(1).name = 'OSPN_D1';    % ORIENTATION pops
s.populations(1).size = multiplicity(1);
s.populations(1).equations = eqns;
s.populations(1).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
s.populations(1).parameters = {'cm',1,'g_l',g_l_D1,'g_cat',g_cat_D1,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_Inp_O,'offset_pfc_poisson',tOff_Inp_O,'f_pfc_poisson',freq_ref_O,'thresholddutycycle_pfc_poisson',thrDUTY,'latency_pfc_poisson',latency_Inp,'DC_pfc_poisson',DC_Inp,'AC_pfc_poisson',AC_ref_O};
s.populations(2).name = 'OSPN_D2';
s.populations(2).size = multiplicity(2);
s.populations(2).equations = eqns;
s.populations(2).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
s.populations(2).parameters = {'cm',1,'g_l',g_l_D2,'g_cat',g_cat_D2,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_Inp_O,'offset_pfc_poisson',tOff_Inp_O,'f_pfc_poisson',freq_ref_O,'thresholddutycycle_pfc_poisson',thrDUTY,'latency_pfc_poisson',latency_Inp,'DC_pfc_poisson',DC_Inp,'AC_pfc_poisson',AC_ref_O};

s.populations(3).name = 'CSPN_D1';    % COLOR pops
s.populations(3).size = multiplicity(3);
s.populations(3).equations = eqns;
s.populations(3).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
s.populations(3).parameters = {'cm',1,'g_l',g_l_D1,'g_cat',g_cat_D1,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_Inp_C,'offset_pfc_poisson',tOff_Inp_C,'f_pfc_poisson',freq_ref_C,'thresholddutycycle_pfc_poisson',thrDUTY,'latency_pfc_poisson',latency_Inp,'DC_pfc_poisson',DC_Inp,'AC_pfc_poisson',AC_ref_C};
s.populations(4).name = 'CSPN_D2';
s.populations(4).size = multiplicity(4);
s.populations(4).equations = eqns;
s.populations(4).mechanism_list = {'spn_iNa','spn_iK','spn_iLeak','spn_iM','spn_iCa','spn_CaBuffer','spn_iKca','spn_iPoisson','spn_ipfcPoisson'};
s.populations(4).parameters = {'cm',1,'g_l',g_l_D2,'g_cat',g_cat_D2,'g_poisson',g_extInp,'baseline_poisson',baseline_extInp,'g_pfc_poisson',g_pfcInp,'onset_pfc_poisson',tOn_Inp_C,'offset_pfc_poisson',tOff_Inp_C,'f_pfc_poisson',freq_ref_C,'thresholddutycycle_pfc_poisson',thrDUTY,'latency_pfc_poisson',latency_Inp,'DC_pfc_poisson',DC_Inp,'AC_pfc_poisson',AC_ref_C};

%%%%%%%%%%% %%% Connectivity definition %%% %%%%%%%%%%%%%
s.connections(1).direction = 'OSPN_D1->OSPN_D1';   % Model INTRACONNECTIONS: ORIENTATION
s.connections(1).mechanism_list = {'spn_iGABA'};
s.connections(1).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD1,'g_gaba',g_s_D1toD1,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(2).direction = 'OSPN_D1->OSPN_D2';
s.connections(2).mechanism_list = {'spn_iGABA'};
s.connections(2).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD2,'g_gaba',g_s_D1toD2,'tau_gaba', tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(3).direction = 'OSPN_D2->OSPN_D1';
s.connections(3).mechanism_list = {'spn_iGABA'};
s.connections(3).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD1,'g_gaba',g_s_D2toD1,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(4).direction = 'OSPN_D2->OSPN_D2';
s.connections(4).mechanism_list = {'spn_iGABA'};
s.connections(4).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD2,'g_gaba',g_s_D2toD2,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(5).direction = 'CSPN_D1->CSPN_D1';    % Model INTRACONNECTIONS: COLOUR
s.connections(5).mechanism_list = {'spn_iGABA'};
s.connections(5).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD1,'g_gaba',g_s_D1toD1,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(6).direction = 'CSPN_D1->CSPN_D2';
s.connections(6).mechanism_list = {'spn_iGABA'};
s.connections(6).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD2,'g_gaba',g_s_D1toD2,'tau_gaba', tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(7).direction = 'CSPN_D2->CSPN_D1';
s.connections(7).mechanism_list = {'spn_iGABA'};
s.connections(7).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD1,'g_gaba',g_s_D2toD1,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(8).direction = 'CSPN_D2->CSPN_D2';
s.connections(8).mechanism_list = {'spn_iGABA'};
s.connections(8).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD2,'g_gaba',g_s_D2toD2,'tau_gaba', tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(9).direction='OSPN_D1->CSPN_D1';   % Model INTERCONNECTIONS: ORIENTATION
s.connections(9).mechanism_list='spn_iGABA';
s.connections(9).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD1,'g_gaba',g_ns_D1toD1,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(10).direction='OSPN_D1->CSPN_D2';
s.connections(10).mechanism_list='spn_iGABA';
s.connections(10).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD2,'g_gaba',g_ns_D1toD2,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(11).direction='OSPN_D2->CSPN_D1';
s.connections(11).mechanism_list='spn_iGABA';
s.connections(11).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD1,'g_gaba',g_ns_D2toD1,'tau_gaba',tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(12).direction='OSPN_D2->CSPN_D2';
s.connections(12).mechanism_list='spn_iGABA';
s.connections(12).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD2,'g_gaba',g_ns_D2toD2,'tau_gaba',tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(13).direction='CSPN_D1->OSPN_D1';   % Model INTERCONNECTIONS: COLOUR
s.connections(13).mechanism_list='spn_iGABA';
s.connections(13).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD1,'g_gaba',g_ns_D1toD1,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(14).direction='CSPN_D1->OSPN_D2';
s.connections(14).mechanism_list='spn_iGABA';
s.connections(14).parameters = {'conn_prob_gaba',conn_prob_gaba_D1toD2,'g_gaba',g_ns_D1toD2,'tau_gaba',tau_gaba_D1,'tauD_gaba',tauD_D1_highDA,'deltaD_gaba',deltaD_gabaD1,'tau_riseD_gaba',tau_riseD};

s.connections(15).direction='CSPN_D2->OSPN_D1';
s.connections(15).mechanism_list='spn_iGABA';
s.connections(15).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD1,'g_gaba',g_ns_D2toD1,'tau_gaba',tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

s.connections(16).direction='CSPN_D2->OSPN_D2';
s.connections(16).mechanism_list='spn_iGABA';
s.connections(16).parameters = {'conn_prob_gaba',conn_prob_gaba_D2toD2,'g_gaba',g_ns_D2toD2,'tau_gaba',tau_gaba_D2,'tauD_gaba',tauD_D2_highDA,'deltaD_gaba',deltaD_gabaD2,'tau_riseD_gaba',tau_riseD};

%%%%%%%%%%% %%% Simulation definition %%% %%%%%%%%%%%%%
sim_log_flag = 1;       % whether to display simulation progress
data = dsSimulate(s,'time_limits',tspan,'dt',0.05,'solver','rk4','sim_log_flag',sim_log_flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting Resuls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors_OC = [
            254/255    153/255   82/255     % naranja
            47/255     191/255   246/255];  % azul clarito
xlabel_Ds = {'','0','0.5','1'};
xlabel_OC = {'-1.5','-1','-0.5','0','0.5','1'};

if strcmp(pfc_type, 'Spontaneous')
%%%%%%%%%%%% %%% Spont %%% %%%%%%%%%%%%%%%%
    % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp tOff_Inp];
    voltage_spont{1} = data.OSPN_D1_V;
    voltage_spont{2} = data.OSPN_D2_V;
    voltage_spont{3} = data.CSPN_D1_V;
    voltage_spont{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D2O, and then also, D1C Vs. D2C
    plotInstFR(voltage_spont, time, transient, tl, xlabel_Ds);

elseif strcmp(pfc_type, 'DC')
%%%%%%%%%%%%% %%% DC %%% %%%%%%%%%%%%%%%%
     % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp tOff_Inp];
    voltage_DC{1} = data.OSPN_D1_V;
    voltage_DC{2} = data.OSPN_D2_V;
    voltage_DC{3} = data.CSPN_D1_V;
    voltage_DC{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D2O, and then also, D1C Vs. D2C
    plotInstFR(voltage_DC, time, transient, tl, xlabel_Ds);

elseif strcmp(pfc_type, 'Alpha')
%%%%%%%%%%%%% %%% AC alpha %%% %%%%%%%%%%%%%%%
    % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp tOff_Inp];
    voltage_alpha{1} = data.OSPN_D1_V;
    voltage_alpha{2} = data.OSPN_D2_V;
    voltage_alpha{3} = data.CSPN_D1_V;
    voltage_alpha{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D2O, and then also, D1C Vs. D2C
    plotInstFR(voltage_alpha, time, transient, tl, xlabel_Ds);

elseif strcmp(pfc_type, 'Beta')
%%%%%%%%%%%%% %%% AC beta %%% %%%%%%%%%%%%%%%
    % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp tOff_Inp];
    voltage_beta{1} = data.OSPN_D1_V;
    voltage_beta{2} = data.OSPN_D2_V;
    voltage_beta{3} = data.CSPN_D1_V;
    voltage_beta{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D2O, and then also, D1C Vs. D2C
    plotInstFR(voltage_beta, time, transient, tl, xlabel_Ds);

elseif strcmp(pfc_type, 'Orientation')
    % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp tOff_Inp];
    Ovoltage_btw{1} = data.OSPN_D1_V;
    Ovoltage_btw{2} = data.CSPN_D1_V;
    Ovoltage_btw{3} = data.OSPN_D2_V;
    Ovoltage_btw{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D1C, and then also, D2O Vs. D2C
    Ovoltage{1} = data.OSPN_D1_V;
    Ovoltage{2} = data.OSPN_D2_V;
    Ovoltage{3} = data.CSPN_D1_V;
    Ovoltage{4} = data.CSPN_D2_V;

    plotInstFR(Ovoltage_btw, time, transient, tl, xlabel_OC, colors_OC);
    plotInstFR(Ovoltage, time, transient, tl, xlabel_OC);
else
    % keep on eye on the order of voltage variable 'cause depending on the order it will show different comparations.
    time = data.time;
    tl = [tOn_Inp_C tOff_Inp_C];
    Cvoltage_btw{1} = data.OSPN_D1_V;
    Cvoltage_btw{2} = data.CSPN_D1_V;
    Cvoltage_btw{3} = data.OSPN_D2_V;
    Cvoltage_btw{4} = data.CSPN_D2_V;
    % On this case, it will show compared results for D1O Vs. D1C, and then also, D2O Vs. D2C
    Cvoltage{1} = data.OSPN_D1_V;
    Cvoltage{2} = data.OSPN_D2_V;
    Cvoltage{3} = data.CSPN_D1_V;
    Cvoltage{4} = data.CSPN_D2_V;

    plotInstFR(Cvoltage_btw, time, transient, tl, xlabel_OC, colors_OC);
    plotInstFR(Cvoltage, time, transient, tl, xlabel_OC);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Saving Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment it to save the previous figures.

% for iFig = 1:length(findobj('type','figure'))
%     saveas(figure(iFig),sprintf('%s_trial_0%d.svg',pfc_type,iFig),'svg')
% end
