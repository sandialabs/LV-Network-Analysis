%% main_run_lv_network_hca
% About
% This script is for running PV hosting capactiy analysis (HCA) on meshed LV
% distribution networks using one of two different methods, then plots the
% results:
%
% 1) The detailed method determines that maximum PV size that can be installed at any 3-phase LV
%    bus without causing reverse current to exceed a predefined threshold (i.e., revThresh) at any network
%    protector relay location under a given loading level (i.e., loadLevel).
% 2) The streamlined method does not require power flow analysis and
%    determines the hosting capacity results based on characteristics of the network model

%% Start
clear;
currentDirectory = cd;

dirSave = currentDirectory;

%% Check that GridPV Ver 2.2 is installed
infoMATLAB = ver;
filterGridPV = strcmpi('GridPV',{infoMATLAB.Name}');
if ~any(filterGridPV)
    error(sprintf('\n===This script requires the GridPV Toolbox (Ver 2.2), which can be downloaded for free at:\n https://pvpmc.sandia.gov/tools/gridpv-toolbox/ \n\n'));
elseif ~strcmpi(infoMATLAB(filterGridPV).Version,'2.2')
    error(sprintf('\n===This script requires the GridPV Toolbox (Ver 2.2), which can be downloaded for free at:\n https://pvpmc.sandia.gov/tools/gridpv-toolbox/ \n\n'));
end

%% Define the OpenDSS model directory
% A version of the IEEE 342-node LV Network Test System is included with
% OpenDSS, and is used as an example with this code. 

cd('C:\Program Files\OpenDSS\IEEETestCases\LVTestCaseNorthAmerican');
dirDSS = cd;

%% Initialize OpenDSS
cd(dirDSS);
fprintf('=====Starting OpenDSS...(%s)\n\n',datestr(now))
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircObj.AllowForms = false;
DSSCircuit = DSSCircObj.ActiveCircuit;

%% Compile Reduced IEEE 342-node Model
cd(dirDSS);
DSSCircObj.DataPath = dirDSS;
fprintf('==Compiling Circuit...(%s)\n',datestr(now()));
DSSText.command = 'ClearAll';
DSSText.command = 'Compile "Master.dss"';
fprintf('==Finished...(%s)\n\n',datestr(now()));
cd(currentDirectory);

DSSText.Command = 'set mode=direct';
DSSText.Command = 'solve';
DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
DSSText.Command = 'solve';

% Adding an energymeter to the example circuit (the next three lines can be
% commented out if using a different circuit model)
DSSText.Command = 'new energymeter.IEEE342 element=line.1';
DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
DSSText.Command = 'solve';


%% Run the Detailed LV HCA
Results = lv_network_hca_detailed(DSSCircObj,'revThresh',[0]','loadLevel',100,'Vpu_UB',1.05,'Vpu_LB',0.95,'ShowPlots',true,'CalcVoltage',true,'CalcProtection',true,'CalcThermal',true,'PV_limit',2.25);

% Run Streamlined HCA

%% Save Results
cd(dirSave);
filename = sprintf('detailedHCA_%s%s.mat',datestr(now,30),DSSCircuit.Name);
save(filename,'Results','-v7.3');
cd(currentDirectory);

%% End of script
fprintf('\n\n===Simulation Finished...(%s)===\n\n\n',datestr(now));




























