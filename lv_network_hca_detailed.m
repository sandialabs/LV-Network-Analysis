%% lv_network_hca_detailed
% This script runs a detailed locational PV hosting capactiy analysis on
% a meshed LV distribution network. The script determines the maximum PV
% size that can be installed at any 3-phase LV bus without causing reverse
% current to exceed a predefined threshold (i.e., revThresh) at any network
% protector relay location under a given loading level (i.e., loadLevel).
% The hosting capacity results are then plotted for each possible location.
% 
% This function assumes there are network protector (NP) relays on the LV side of each network transformer.
% Typically, the NPs are set to trip for any amount of reverse current,
% which significantly limits the DER installation capacity since reverse
% current caused by DER would result in false trips. 
% The NP can be desensitized by setting revThresh above 0. In general,
% DER hosting capacity will increase as this variable increases. 
%
%% Syntax
%  Results = LVnetworkHCA_detailed(DSSCircObj);
%  Results = LVnetworkHCA_detailed(DSSCircObj,'revThresh',100);
%
%% Description
% Function to run a detailed locational hosting capacity analysis on a meshed LV network, then plots the results. 
% The user can set the reverse current threshold and loading level. By default, all types are plotted.
%
%% Inputs
% * *|DSSCircObj|* - link to OpenDSS active circuit and command text (from DSSStartup)
% * *|revThresh|* - parameter variable (numeric) for setting for NP reverse current trip threshold (default = 0). 
% * *|loadLevel|* - parameter variable (numeric) to set the loading level on the LV network, scales all loads between 0-125% of original values.
% * *|ShowPlots|* - parameter variable (logical) to toggle the plotting on (default) or off

%% Outputs
% *|Results|* is a structure with the simulation settings and HC results for each LV bus.

%% Copyright 2018
% Sandia Corporation. Under the terms of Contract XXXX with Sandia Corporation, the U.S. Government retains certain rights in this software.
% See the license agreement for full terms and conditions.
%
% Please acknowledge any contributions of the GridPV Toolbox by citing:
% M. J. Reno and K. Coogan, "Grid Integrated Distributed PV (GridPV) Version 2," Sandia National Laboratories SAND2013-20141, 2014.
%
%% Example
% % Runs the HCA on the LV network, plots the results, and returns them in the Results struct

%%
% [DSSCircObj, DSSText, gridpvPath] = DSSStartup;
% DSSText.command = ['Compile "' gridpvPath 'ExampleCircuit\master_Ckt24.dss"'];
% DSSText.command = 'solve';
% Results =
% LVnetworkHCA_detailed(DSSCircObj,'revThresh',1000,'loadLevel',50); % run HCA under 50% loading and allowing 1000A of reverse current at all network protector relays.


%% Start of Function
function Results = lv_network_hca_detailed(DSSCircObj,varargin)

%% Parse inputs
p = inputParser; %setup parse structure
p.addRequired('DSSCircObj', @isinterfaceOpenDSS);

%p.addParameter('revThresh',0, @(x)isnumeric(x) && ((x>=0) && (x<=10000)) && size(x,1)==1);
p.addParameter('revThresh', [0]', @(x)isvector(x) && (all(x>=0) && all(x<=20000)));
p.addParameter('loadLevel',100, @(x)isnumeric(x) && ((x>=0) && (x<=125)));
p.addParameter('Vpu_UB',1.05, @(x)isnumeric(x) && ((x>=1) && (x<=1.5)));
p.addParameter('Vpu_LB',0.95, @(x)isnumeric(x) && ((x>=0.5) && (x<=1)));
p.addParameter('ShowPlots',true, @(x)islogical(x));
p.addParameter('HC_bus','file', @(x)ischar(x) && isfile(x))
p.addParameter('Contingency','file', @(x)ischar(x) && isfile(x))
p.addParameter('CalcProtection',true, @(x)islogical(x));
p.addParameter('CalcThermal',true, @(x)islogical(x));
p.addParameter('CalcVoltage',true, @(x)islogical(x));
p.addParameter('PV_limit',1.5, @(x)isnumeric(x));
%%add a parameter for LV bus names to add Matt's bus list, cell array
%%strings of bus names

%%provide options to select between thermal,voltage,protection


p.parse(DSSCircObj, varargin{:}); %parse inputs
nargoutchk(0,1);

allFields = fieldnames(p.Results); %set all parsed inputs to workspace variables
for ii=1:length(allFields)
    eval([allFields{ii}, ' = ', 'p.Results.',allFields{ii},';']);
end

%% Check for GridPV
infoMATLAB = ver;
filterGridPV = strcmpi('GridPV',{infoMATLAB.Name}');
if ~any(filterGridPV)
    error(sprintf('\n===This script requires the GridPV Toolbox (Ver 2.2), which can be downloaded for free at:\n https://pvpmc.sandia.gov/tools/gridpv-toolbox/ \n\n'));
elseif ~strcmpi(infoMATLAB(filterGridPV).Version,'2.2')
    error(sprintf('\n===This script requires the GridPV Toolbox (Ver 2.2), which can be downloaded for free at:\n https://pvpmc.sandia.gov/tools/gridpv-toolbox/ \n\n'));
end

%% Get Circuit Info
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText = DSSCircObj.Text;

Lines = getLineInfo(DSSCircObj);
Load_Original = getLoadInfo(DSSCircObj);
load_Original_All =sum([Load_Original.powerReal]);

Bus_old = getBusInfo(DSSCircObj);
busNames={Bus_old.name}';
Bus= getBusInfo(DSSCircObj,busNames,1);


Xfmr = getTransformerInfo(DSSCircObj);
NetXfmr = Xfmr([Xfmr.bus2kVBase]<1);
SubXfmr = Xfmr([Xfmr.bus2kVBase]>=1);
NetworkProtectors=NetXfmr;
namesNP = {NetworkProtectors.name}';



%%Contingency_scenario with input files
if ~strcmp(Contingency,'file') 
    if ~isfile(Contingency)
        error('The file %s does not exist. Please provide the right file.', filename); %%already checking this while providing inputs, not required
    else
        
        option=detectImportOptions(Contingency,'Delimiter',',', ...            % comma-delimited file
                                        'NumHeaderLines',0, ...         % skip the header lines 
                                        'ExtraColumnsRule','ignore'); 
        Scenarios=readtable(Contingency, option);
    end
else
    Scenarios=0;
    
end

%%HC_bus_list with input files
if ~strcmp(HC_bus,'file') 
    if ~isfile(HC_bus)
        error('The file %s does not exist. Please provide the right file.', filename); %%already checking this while providing inputs, not required
    else
        option=detectImportOptions(HC_bus,'Delimiter',',', ...            % comma-delimited file
                                        'NumHeaderLines',0, ...         % skip the header lines 
                                        'ExtraColumnsRule','ignore'); 
        HC_bus_list=readtable(HC_bus, option);
    end
else
    HC_bus_list=0;
    
end



%%Using HC_bus_list instead of Bus_LV if available
if istable(HC_bus_list)
    [tf, idx] = ismember(lower(HC_bus_list.name), {Bus.name}');
    HC_bus_list.kVBase=nan(height(HC_bus_list),1);
    HC_bus_list.coordinates=nan(height(HC_bus_list),2);
    Bus_kV=[Bus.kVBase]';
    Bus_coordinates={Bus.coordinates}';
    Bus_coordinates=cat(1,Bus_coordinates{:});
    HC_bus_list.kVBase(tf) = Bus_kV(idx(tf))';
    HC_bus_list.coordinates(tf)= Bus_coordinates(idx(tf));
    HC_bus_list.coordinates(tf,2)=Bus_coordinates(idx(tf),2);
    Bus_LV=HC_bus_list;
    Bus_LV=table2struct(Bus_LV);

else
    Bus_LV = Bus([Bus.kVBase]'<1);% use it if the user does not pass anything
    NP_Bus1 = {NetworkProtectors.bus1}';
    NP_Bus1 = extractBefore(NP_Bus1,'.');
    
    [tf, idx] = ismember(lower({Bus_LV.name}'),NP_Bus1);
    Bus_LV = Bus_LV(~tf); % not considering installations directly at NP buses
    Bus_LV = Bus_LV([Bus_LV.numPhases]==3); % only consider 3-phase buses
end

%%use when HC bus list has at least 1 bus that is above 1 kV
if istable(HC_bus_list) && any ([Bus_LV.kVBase]'>=1)
    additional_tx=Bus_LV([Bus_LV.kVBase]'>=1); % find the bus name
    [tf, idx] = ismember({SubXfmr.bus2}',lower(additional_tx.name)); % check which transformers are connected to that bus
    indices = find(idx == 1);
    SubXfmr_np=SubXfmr(indices);
    NetXfmr = [NetXfmr' SubXfmr_np']'; %% add them to the network transformer list
    NetworkProtectors=NetXfmr;
    namesNP = {NetworkProtectors.name}';
end


DSSText.Command = 'set mode=direct';
DSSText.Command = 'solve';
DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
DSSText.Command = 'solve';

%% Adjust Loading Scenario
DSSText.Command = 'batchedit load..* vmaxpu=2 vlowpu=0.5 vminpu=0.5 model=1';

if (loadLevel>=0) && (loadLevel<=125)
        for ll = 1:length(Load_Original)
            DSSText.Command = sprintf('edit load.%s kw=%1.5f kvar=%1.5f',Load_Original(ll).name,(loadLevel/100)*(Load_Original(ll).kW),...
                (loadLevel/100)*(Load_Original(ll).kvar));
        end
        DSSText.Command = 'set mode=direct';
        DSSText.Command = 'solve';
        DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
        DSSText.Command = 'solve';
        Load_Now = getLoadInfo(DSSCircObj);
else
    warning('===Load level selection is outside of bounds. Loads will remain at original values.===\n\n');
    DSSText.Command = 'set mode=direct';
    DSSText.Command = 'solve';
    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
    DSSText.Command = 'solve';
    Load_Now = getLoadInfo(DSSCircObj);
end

%% Add PV system to a bus and refresh data

busPV = Bus_LV(1).name;
busPVkV=Bus_LV(1).kVBase;
load_kW_All = sum([Load_Now.powerReal]);
load_kVAR_All = sum([Load_Now.powerReactive]);

penetrationPCT_start = 0;
pen_resolution = 0.1;
penetrationPCT = penetrationPCT_start;
penetration = penetrationPCT/100;

kW_min = 1;

DSSText.Command = sprintf('New PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true %%pmpp=200 kva=%1.5f pmpp=%1.5f varfollowinverter=true vmaxpu=2 vminpu=0.5 irradiance=1',...
    busPV,busPVkV,kW_min,kW_min);

DSSText.Command = 'set mode=direct';
DSSText.Command = 'solve';
DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 ';
DSSText.Command = 'solve';
PV = getPVInfo(DSSCircObj);

%% Pre-Allocate the results variables
if istable(Scenarios)

    HC_kW_revCurrent = nan(length(Bus_LV),width(Scenarios));
    HC_penetration_revCurrent = nan(length(Bus_LV),width(Scenarios));
    
    I1_final = nan(length(namesNP),length(Bus_LV));
    bus_issue_rc = nan(length(Bus_LV),width(Scenarios));
else
    HC_kW_revCurrent = nan(length(Bus_LV),length(revThresh));
    HC_penetration_revCurrent = nan(length(Bus_LV),length(revThresh));
    bus_issue_rc = nan(length(Bus_LV),length(revThresh));
    I1_final = nan(length(namesNP),length(Bus_LV));
end
%% Solve the HCA - Reverse Current Threshold

convergenceThresh = 1;
if CalcProtection
% Conduct HCA for each reverse current threshold
    for rr = 1:length(revThresh)
        h = waitbar(0,sprintf('Solving HCA (%i A reverse current limit)...',revThresh(rr)));
    
        % Find the max PV size at each LV bus without causing reverse current
        % above the revThresh(rr) at any Network Protector location
        
        for bb = 1:length(Bus_LV)
    
            waitbar(bb/length(Bus_LV),h);
            if istable(Scenarios)
                no_of_scenarios=width(Scenarios);
                Scenarios_new=table2struct(Scenarios);
                fields=fieldnames (Scenarios_new);
    
                
                for cc=1:numel(fields)
                
    
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=false',Scenarios_new(dd).(fields{cc}));
                        
                    end
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                    
                    %DSSText.Command = sprintf('batchedit PVsystem..* bus1=%s enabled=false',...
                        %Bus_LV(bb).name);%%%add kv on every case %%%change batchedit to edit
            
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
            
                    
                    Xfmr = getTransformerInfo(DSSCircObj);
                    NetXfmr = Xfmr([Xfmr.bus2kVBase]<1);
                    SubXfmr = Xfmr([Xfmr.bus2kVBase]>=1);
                    if istable(HC_bus_list) && any ([Bus_LV.kVBase]'>=1)

                        additional_tx=Bus_LV([Bus_LV.kVBase]'>=1);
                        [tf, idx] = ismember({SubXfmr.bus2}',lower(additional_tx.name));
                        indices = find(idx == 1);
                        SubXfmr_np=SubXfmr(indices);
                        NetXfmr = [NetXfmr' SubXfmr_np']';
                        
                    end
                    NetworkProtectors=NetXfmr;
                    namesNP = {NetworkProtectors.name}';
                    NP_now = getTransformerInfo(DSSCircObj,namesNP);
                    vseq = {NP_now.bus2CplxSeqVoltages}';
                    vseq = cat(1,vseq{:},[]);
                    v1 = complex(vseq(:,3),vseq(:,4));
                    v1_mag = abs(v1);
                    v1_ang = rad2deg(angle(v1));
            
                    Iseq = {NP_now.bus2CplxSeqCurrents}';
                    Iseq = cat(1,Iseq{:},[]);
                    I1 = complex(Iseq(:,3),Iseq(:,4));
                    I1_mag = abs(I1);
                    I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                    %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                    hasReverseCurrent = cosd(I1_ang)<0;
                    I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                    NP_Imin = min(I1_mag);
                    I_headroom = NP_Imin + revThresh(rr);
            
                    if I_headroom < convergenceThresh
                        HC_kW_revCurrent(bb,cc) = 0;
                        HC_penetration_revCurrent(bb,cc) = 0;
                        %I1_final(:,bb) = I1_mag;
                    else

                        %DSSText.Command = 'batchedit PVsystem..* enabled=true';
                        DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                        Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                        DSSText.Command = 'set mode=direct';
                        DSSText.Command = 'solve';
            
                        DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                        DSSText.Command = 'solve';
            
                        
                        NP_now = getTransformerInfo(DSSCircObj,namesNP);
                        vseq = {NP_now.bus2CplxSeqVoltages}';
                        vseq = cat(1,vseq{:},[]);
                        v1 = complex(vseq(:,3),vseq(:,4));
                        v1_mag = abs(v1);
                        v1_ang = rad2deg(angle(v1));
            
                        Iseq = {NP_now.bus2CplxSeqCurrents}';
                        Iseq = cat(1,Iseq{:},[]);
                        I1 = complex(Iseq(:,3),Iseq(:,4));
                        I1_mag = abs(I1);
                        I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                        %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                        hasReverseCurrent = cosd(I1_ang)<0;
                        I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                        NP_Imin = min(I1_mag);
                        I_headroom = NP_Imin + revThresh(rr);
            
                        iter=1;
                        maxIter=10000;
                        iterSovle = true;
            
                        while  (abs(I_headroom)>convergenceThresh && (iter<=maxIter) && (PV.kVA>0)) && iterSovle %%adding a new condition so simulation stops when PV size is 200% over the total load
                            DSSCircuit.PVSystems.Name = 'PV1';
                            PVold = DSSCircuit.PVSystems.kVArated;
            
                            if revThresh(rr)<1500
                                penetrationPCT = penetrationPCT+pen_resolution;
                                P_pv_new = (penetrationPCT/100)*load_kW_All ;
                            else
                                k_converge = 0.5;
            
                                deltaPV = k_converge*((3*I_headroom*PV.voltage) / 1e3);
            
                                P_pv_new = max((PVold + deltaPV),0);
                            end
                            DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                P_pv_new,P_pv_new);

            
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
            
                            DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                            DSSText.Command = 'solve';
            
                            if ~DSSCircuit.Solution.Converged
                                warnSt = circuitCheck(DSSCircObj);
                                bus_issue_rc(bb,cc)=1;
                                DSSText.Command = 'set mode=direct';
                                DSSText.Command = 'solve';
                                %error('Circuit did not converge...\n\n');
                            end
                            NP_now = getTransformerInfo(DSSCircObj,namesNP);
                            %NP_now = getLineInfo(DSSCircObj,namesNP);
                            vseq = {NP_now.bus2CplxSeqVoltages}';
                            vseq = cat(1,vseq{:},[]);
                            v1 = complex(vseq(:,3),vseq(:,4));
                            v1_mag = abs(v1);
                            v1_ang = rad2deg(angle(v1));
            
                            Iseq = {NP_now.bus2CplxSeqCurrents}';
                            Iseq = cat(1,Iseq{:},[]);
                            I1 = complex(Iseq(:,3),Iseq(:,4));
                            I1_mag = abs(I1);
                            I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                            %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                            hasReverseCurrent = cosd(I1_ang)<0;
                            I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                            NP_Imin = min(I1_mag);
                            I_headroom = NP_Imin + revThresh(rr);
            
                            if ((I_headroom<0) && (revThresh(rr)<1500)) || (PVold>PV_limit*load_kW_All)
                                iterSovle=false;
                                DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                    PVold,PVold);
                                % DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
                                %     PVold,PVold);
                                DSSText.Command = 'set mode=direct';
                                DSSText.Command = 'solve';
            
                                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                                DSSText.Command = 'solve';
                            end
                            iter=iter+1;
                        end
            
                        DSSCircuit.PVSystems.Name = 'PV1';
                        PVold = DSSCircuit.PVSystems.kVArated;
                        HC_kW_revCurrent(bb,cc) = PVold;
                        HC_penetration_revCurrent(bb,cc) = (PVold/load_kW_All)*100;
                        % HC_kW_revCurrent(bb,rr) = PVold;
                        % HC_penetration_revCurrent(bb,rr) = (PVold/load_kW_All)*100;
            
                        %I1_final(:,bb) = I1_mag;
            
                        if iter > maxIter
                            warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                        end
                    end
            
                    % reset the pv system size and simulation parameters
                    penetrationPCT = penetrationPCT_start;
                    penetration = penetrationPCT/100;
                    DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                        kW_min,kW_min);
                    % DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
                    %     kW_min,kW_min);
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=true',Scenarios_new(dd).(fields{cc}));
                    end
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                    DSSText.Command = 'solve';
                end
            else
                
                DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                    
                    %DSSText.Command = sprintf('batchedit PVsystem..* bus1=%s enabled=false',...
                        %Bus_LV(bb).name);%%%add kv on every case %%%change batchedit to edit
            
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                DSSText.Command = 'solve';
            
                %NP_now = getLineInfo(DSSCircObj,namesNP);
                Xfmr = getTransformerInfo(DSSCircObj);
                NetXfmr = Xfmr([Xfmr.bus2kVBase]<1);
                SubXfmr = Xfmr([Xfmr.bus2kVBase]>=1);
                if istable(HC_bus_list) && any ([Bus_LV.kVBase]'>=1)

                    additional_tx=Bus_LV([Bus_LV.kVBase]'>=1);
                    [tf, idx] = ismember({SubXfmr.bus2}',lower(additional_tx.name));
                    indices = find(idx == 1);
                    SubXfmr_np=SubXfmr(indices);
                    NetXfmr = [NetXfmr' SubXfmr_np']';
                    
                end
                NetworkProtectors=NetXfmr;
                namesNP = {NetworkProtectors.name}';
                NP_now = getTransformerInfo(DSSCircObj,namesNP);
                vseq = {NP_now.bus2CplxSeqVoltages}';
                vseq = cat(1,vseq{:},[]);
                v1 = complex(vseq(:,3),vseq(:,4));
                v1_mag = abs(v1);
                v1_ang = rad2deg(angle(v1));
        
                Iseq = {NP_now.bus2CplxSeqCurrents}';
                Iseq = cat(1,Iseq{:},[]);
                I1 = complex(Iseq(:,3),Iseq(:,4));
                I1_mag = abs(I1);
                I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                hasReverseCurrent = cosd(I1_ang)<0;
                I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                NP_Imin = min(I1_mag);
                I_headroom = NP_Imin + revThresh(rr);
            
                if I_headroom < convergenceThresh
                    HC_kW_revCurrent(bb,rr) = 0;
                    HC_penetration_revCurrent(bb,rr) = 0;
                    I1_final(:,bb) = I1_mag;
                else
                    % for out=1:length(lines_out)
                    %     DSSText.Command ='Disable Line.%s')
                    %DSSText.Command = 'batchedit PVsystem..* enabled=true';
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
        
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
        
                    %NP_now = getLineInfo(DSSCircObj,namesNP);
                    NP_now = getTransformerInfo(DSSCircObj,namesNP);
                    vseq = {NP_now.bus2CplxSeqVoltages}';
                    vseq = cat(1,vseq{:},[]);
                    v1 = complex(vseq(:,3),vseq(:,4));
                    v1_mag = abs(v1);
                    v1_ang = rad2deg(angle(v1));
        
                    Iseq = {NP_now.bus2CplxSeqCurrents}';
                    Iseq = cat(1,Iseq{:},[]);
                    I1 = complex(Iseq(:,3),Iseq(:,4));
                    I1_mag = abs(I1);
                    I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                    %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                    hasReverseCurrent = cosd(I1_ang)<0;
                    I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                    NP_Imin = min(I1_mag);
                    I_headroom = NP_Imin + revThresh(rr);
        
                    iter=1;
                    maxIter=10000;
                    iterSovle = true;
        
                    while  (abs(I_headroom)>convergenceThresh && (iter<=maxIter) && (PV.kVA>0)) && iterSovle
                        DSSCircuit.PVSystems.Name = 'PV1';
                        PVold = DSSCircuit.PVSystems.kVArated;
        
                        if revThresh(rr)<1500
                            penetrationPCT = penetrationPCT+pen_resolution;
                            P_pv_new = (penetrationPCT/100)*load_kW_All ;
                        else
                            k_converge = 0.5;
        
                            deltaPV = k_converge*((3*I_headroom*PV.voltage) / 1e3);
        
                            P_pv_new = max((PVold + deltaPV),0);
                        end
                        DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                            P_pv_new,P_pv_new);
                        % DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
                        %     P_pv_new,P_pv_new);
        
                        DSSText.Command = 'set mode=direct';
                        DSSText.Command = 'solve';
        
                        DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                        DSSText.Command = 'solve';
        
                        if ~DSSCircuit.Solution.Converged
                            warnSt = circuitCheck(DSSCircObj);
                            bus_issue_rc(bb,rr)=1;
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
                            %error('Circuit did not converge...\n\n');
                        end
                        NP_now = getTransformerInfo(DSSCircObj,namesNP);
                        %NP_now = getLineInfo(DSSCircObj,namesNP);
                        vseq = {NP_now.bus2CplxSeqVoltages}';
                        vseq = cat(1,vseq{:},[]);
                        v1 = complex(vseq(:,3),vseq(:,4));
                        v1_mag = abs(v1);
                        v1_ang = rad2deg(angle(v1));
        
                        Iseq = {NP_now.bus2CplxSeqCurrents}';
                        Iseq = cat(1,Iseq{:},[]);
                        I1 = complex(Iseq(:,3),Iseq(:,4));
                        I1_mag = abs(I1);
                        I1_ang = -180 + rad2deg(angle(I1)) - v1_ang; %%my line to use bus 2
                        %I1_ang = rad2deg(angle(I1)) - v1_ang;     %degrees (negative sign = lagging)
                        hasReverseCurrent = cosd(I1_ang)<0;
                        I1_mag(hasReverseCurrent) = I1_mag(hasReverseCurrent)*-1;
                        NP_Imin = min(I1_mag);
                        I_headroom = NP_Imin + revThresh(rr);
        
                        if ((I_headroom<0) && (revThresh(rr)<1500)) || (PVold>PV_limit*load_kW_All)
                            iterSovle=false;
                            DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                PVold,PVold);
                            % DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
                            %     PVold,PVold);
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
        
                            DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                            DSSText.Command = 'solve';
                        end
                        iter=iter+1;
                    end
        
                    DSSCircuit.PVSystems.Name = 'PV1';
                    PVold = DSSCircuit.PVSystems.kVArated;
                    HC_kW_revCurrent(bb,rr) = PVold;
                    HC_penetration_revCurrent(bb,rr) = (PVold/load_kW_All)*100;
                    % HC_kW_revCurrent(bb,rr) = PVold;
                    % HC_penetration_revCurrent(bb,rr) = (PVold/load_kW_All)*100;
            
                        %I1_final(:,bb) = I1_mag;
            
                    if iter > maxIter
                        warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                    end
                end
        
                % reset the pv system size and simulation parameters
                penetrationPCT = penetrationPCT_start;
                penetration = penetrationPCT/100;
                DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                    kW_min,kW_min);
                % DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
                %     kW_min,kW_min);
    
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                DSSText.Command = 'solve';
            end
    
            
        end
        close(h);
    end
end

%% Solve the HCA - Thermal (DER Gen Limit)
if istable(Scenarios)

    HC_kW_Thermal_DERgen = nan(length(Bus_LV),width(Scenarios));
    HC_penetration_Thermal_DERgen = nan(length(Bus_LV),width(Scenarios));
    bus_issue_th = nan(length(Bus_LV),width(Scenarios));
    
else
    HC_kW_Thermal_DERgen = nan(length(Bus_LV),length(revThresh));
    HC_penetration_Thermal_DERgen = nan(length(Bus_LV),length(revThresh));
    bus_issue_th = nan(length(Bus_LV),length(revThresh));
end


PV_kW_min = 100;
stepSize_kW = 100;

if CalcThermal
    for rr = 1:length(revThresh)%%%we don't need it
        h = waitbar(0,'Solving HCA (DER Gen. Thermal Limit)...');
    
        % Find the max PV size at each LV bus without causing thermal overload
        % at any network transformer
        for bb = 1:length(Bus_LV)
            waitbar(bb/length(Bus_LV),h);
            if istable(Scenarios)
    
                no_of_scenarios=width(Scenarios);
                Scenarios_new=table2struct(Scenarios);
                fields=fieldnames (Scenarios_new);
                for cc=1:numel(fields)
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=false',Scenarios_new(dd).(fields{cc}));                    
                    end
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f kva=%1.5f pmpp=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase,PV_kW_min,PV_kW_min);
    
    
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
                    
                    Xfmr = getTransformerInfo(DSSCircObj);
                    NetXfmr = Xfmr([Xfmr.bus2kVBase]<1);
                    SubXfmr = Xfmr([Xfmr.bus2kVBase]>=1);
                    if istable(HC_bus_list) && any ([Bus_LV.kVBase]'>=1)

                        additional_tx=Bus_LV([Bus_LV.kVBase]'>=1);
                        [tf, idx] = ismember({SubXfmr.bus2}',lower(additional_tx.name));
                        indices = find(idx == 1);
                        SubXfmr_np=SubXfmr(indices);
                        NetXfmr = [NetXfmr' SubXfmr_np']';
                        
                    end
                    xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
            
                    wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                    wdg2_kva_rating = [xfmr_now.wdg2kVA]';
            
                    bus1Q = [xfmr_now.bus1PowerReactive]';
                    bus1P = [xfmr_now.bus1PowerReal]';
                    bus1S = sqrt(bus1P.^2+bus1Q.^2);
            
                    bus2Q = [xfmr_now.bus2PowerReactive]';
                    bus2P = [xfmr_now.bus2PowerReal]';
                    bus2S = sqrt(bus2P.^2+bus2Q.^2);
    
                    if any((bus1S>wdg1_kva_rating) | (bus2S>wdg2_kva_rating))
                        HC_kW_Thermal_DERgen(bb,cc) = 0;
            
                    else
                        DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                        Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                        %DSSText.Command = sprintf('batchedit PVsystem..* enabled=true');
                        DSSText.Command = 'set mode=direct';
                        DSSText.Command = 'solve';
            
                        DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                        DSSText.Command = 'solve';
            
            
                        xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
            
                        wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                        wdg2_kva_rating = [xfmr_now.wdg2kVA]';
            
                        bus1Q = [xfmr_now.bus1PowerReactive]';
                        bus1P = [xfmr_now.bus1PowerReal]';
                        bus1S = sqrt(bus1P.^2+bus1Q.^2);
            
                        bus2Q = [xfmr_now.bus2PowerReactive]';
                        bus2P = [xfmr_now.bus2PowerReal]';
                        bus2S = sqrt(bus2P.^2+bus2Q.^2);
            
            
            
                        if any((bus1S>wdg1_kva_rating) | (bus2S>wdg2_kva_rating))
                            HC_kW_Thermal_DERgen(bb,cc) = 0;
                        else
            
                            pv_old = PV_kW_min;
            
                            iter=1;
                            maxIter=1000;
                            foundHC = false;
            
                            while  all(bus1S<=wdg1_kva_rating) && all(bus2S<=wdg2_kva_rating) && (iter<=maxIter) && ~foundHC %%&& (pv_old<PV_limit*load_kW_All)
            
                                pvnew = pv_old + stepSize_kW;
            
                                DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                    pvnew,pvnew);
            
                                DSSText.Command = 'set mode=direct';
                                DSSText.Command = 'solve';
            
                                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                                DSSText.Command = 'solve';
            
                                if ~DSSCircuit.Solution.Converged
                                    warnSt = circuitCheck(DSSCircObj);
                                    bus_issue_th(bb,cc)=1;
                                    DSSText.Command = 'set mode=direct';
                                    DSSText.Command = 'solve';
                                    %error('Circuit did not converge...\n\n');
                                end
            
                                xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
            
                                wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                                wdg2_kva_rating = [xfmr_now.wdg2kVA]';
            
                                bus1Q = [xfmr_now.bus1PowerReactive]';
                                bus1P = [xfmr_now.bus1PowerReal]';
                                bus1S = sqrt(bus1P.^2+bus1Q.^2);
            
                                bus2Q = [xfmr_now.bus2PowerReactive]';
                                bus2P = [xfmr_now.bus2PowerReal]';
                                bus2S = sqrt(bus2P.^2+bus2Q.^2);
            
                                if all(bus1S<=wdg1_kva_rating) && all(bus2S<=wdg2_kva_rating)
                                    pv_old = pvnew;
                                    if (pv_old>PV_limit*load_kW_All)
                                        HC_kW_Thermal_DERgen(bb,cc) = pv_old;
                                        foundHC = true;
                                    end
                                else
                                    HC_kW_Thermal_DERgen(bb,cc) = pv_old;
                                    foundHC = true;
                                end
            
                                iter=iter+1;
            
                            end
            
                            if iter > maxIter
                                warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                            end
            
                        end
                    end
            
                    % reset the pv system size and simulation parameters
                    DSSText.Command = sprintf('Edit PVsystem.PV1 enabled=false kva=%1.5f pmpp=%1.5f',...
                        PV_kW_min,PV_kW_min);
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=true',Scenarios_new(dd).(fields{cc}));
                    end
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                    DSSText.Command = 'solve';
                end
            else
                
                DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f kva=%1.5f pmpp=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase,PV_kW_min,PV_kW_min);
    
    
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                DSSText.Command = 'solve';
                
                Xfmr = getTransformerInfo(DSSCircObj);
                NetXfmr = Xfmr([Xfmr.bus2kVBase]<1);
                SubXfmr = Xfmr([Xfmr.bus2kVBase]>=1);
                if istable(HC_bus_list) && any ([Bus_LV.kVBase]'>=1)
                    additional_tx=Bus_LV([Bus_LV.kVBase]'>=1);
                    [tf, idx] = ismember({SubXfmr.bus2}',lower(additional_tx.name));
                    indices = find(idx == 1);
                    SubXfmr_np=SubXfmr(indices);
                    NetXfmr = [NetXfmr' SubXfmr_np']';
                    
                end
                xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
        
                wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                wdg2_kva_rating = [xfmr_now.wdg2kVA]';
        
                bus1Q = [xfmr_now.bus1PowerReactive]';
                bus1P = [xfmr_now.bus1PowerReal]';
                bus1S = sqrt(bus1P.^2+bus1Q.^2);
        
                bus2Q = [xfmr_now.bus2PowerReactive]';
                bus2P = [xfmr_now.bus2PowerReal]';
                bus2S = sqrt(bus2P.^2+bus2Q.^2);
    
                if any((bus1S>wdg1_kva_rating) | (bus2S>wdg2_kva_rating))
                    HC_kW_Thermal_DERgen(bb,rr) = 0;
        
                else
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                    %DSSText.Command = sprintf('batchedit PVsystem..* enabled=true');
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
        
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
        
        
                    xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
        
                    wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                    wdg2_kva_rating = [xfmr_now.wdg2kVA]';
        
                    bus1Q = [xfmr_now.bus1PowerReactive]';
                    bus1P = [xfmr_now.bus1PowerReal]';
                    bus1S = sqrt(bus1P.^2+bus1Q.^2);
        
                    bus2Q = [xfmr_now.bus2PowerReactive]';
                    bus2P = [xfmr_now.bus2PowerReal]';
                    bus2S = sqrt(bus2P.^2+bus2Q.^2);
        
        
        
                    if any((bus1S>wdg1_kva_rating) | (bus2S>wdg2_kva_rating))
                        HC_kW_Thermal_DERgen(bb,rr) = 0;
                    else
        
                        pv_old = PV_kW_min;
        
                        iter=1;
                        maxIter=1000;
                        foundHC = false;
        
                        while  all(bus1S<=wdg1_kva_rating) && all(bus2S<=wdg2_kva_rating) && (iter<=maxIter) && ~foundHC %%&& (pv_old<PV_limit*load_kW_All)
        
                            pvnew = pv_old + stepSize_kW;
        
                            DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                pvnew,pvnew);
        
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
        
                            DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                            DSSText.Command = 'solve';
        
                            if ~DSSCircuit.Solution.Converged
                                warnSt = circuitCheck(DSSCircObj);
                                bus_issue_th(bb,rr)=1;
                                DSSText.Command = 'set mode=direct';
                                DSSText.Command = 'solve';
                                %error('Circuit did not converge...\n\n');
                            end
        
                            xfmr_now = getTransformerInfo(DSSCircObj,{NetXfmr.name}');
        
                            wdg1_kva_rating = [xfmr_now.wdg1kVA]';
                            wdg2_kva_rating = [xfmr_now.wdg2kVA]';
        
                            bus1Q = [xfmr_now.bus1PowerReactive]';
                            bus1P = [xfmr_now.bus1PowerReal]';
                            bus1S = sqrt(bus1P.^2+bus1Q.^2);
        
                            bus2Q = [xfmr_now.bus2PowerReactive]';
                            bus2P = [xfmr_now.bus2PowerReal]';
                            bus2S = sqrt(bus2P.^2+bus2Q.^2);
        
                            if all(bus1S<=wdg1_kva_rating) && all(bus2S<=wdg2_kva_rating)
                                pv_old = pvnew;
                                if (pv_old>PV_limit*load_kW_All)
                                    HC_kW_Thermal_DERgen(bb,rr) = pv_old;
                                    foundHC = true;
                                end
                            else
                                HC_kW_Thermal_DERgen(bb,rr) = pv_old;
                                foundHC = true;
                            end
        
                            iter=iter+1;
        
                        end
        
                        if iter > maxIter
                            warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                        end
                    end
                end
    
            
                % reset the pv system size and simulation parameters
                DSSText.Command = sprintf('Edit PVsystem.PV1 enabled=false kva=%1.5f pmpp=%1.5f',...
                PV_kW_min,PV_kW_min);
    
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                DSSText.Command = 'solve';
            end
        end
    
        close(h);
    end
end


%% Solve the HCA - Voltage Upper Limit
if istable(Scenarios)

    HC_kW_Vlimit_DERgen = nan(length(Bus_LV),width(Scenarios));
    HC_penetration_Vlimit_DERgen = nan(length(Bus_LV),width(Scenarios));
    bus_issue = nan(length(Bus_LV),width(Scenarios));
    
    
else
    HC_kW_Vlimit_DERgen = nan(length(Bus_LV),length(revThresh));
    HC_penetration_Vlimit_DERgen = nan(length(Bus_LV),length(revThresh));
    bus_issue = nan(length(Bus_LV),length(revThresh));
    
end



PV_kW_min = 100;
if CalcVoltage
    for rr = 1:length(revThresh)
        h = waitbar(0,'Solving HCA (Voltage Upper Limit)...');
    
        % Find the max PV size at each LV bus without causing voltage
        % violations (upper bound) at any LV network node
        for bb = 1:length(Bus_LV)
            waitbar(bb/length(Bus_LV),h);
            if istable(Scenarios)
                no_of_scenarios=width(Scenarios);
                Scenarios_new=table2struct(Scenarios);
                fields=fieldnames (Scenarios_new);
    
                for cc=1:numel(fields)
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=false',Scenarios_new(dd).(fields{cc}));
                        
                    end
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f kva=%1.5f pmpp=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase,PV_kW_min,PV_kW_min);
                    % DSSText.Command = sprintf('batchedit PVsystem..* bus1=%s enabled=false kva=%1.5f pmpp=%1.5f',...
                    % Bus_LV(bb).name,PV_kW_min,PV_kW_min);
    
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
            
                    Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                    vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                    vpu_max_noPV = max(vpu_temp);
            
                    DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                    Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                    %DSSText.Command = sprintf('batchedit PVsystem..* enabled=true');
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                    DSSText.Command = 'solve';
            
                    Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                    vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                    vpu_max = max(vpu_temp);
            
                    if any(vpu_temp>Vpu_UB)
                        HC_kW_Vlimit_DERgen(bb,cc) = 0;
                    else
            
                        pv_old = PV_kW_min;
                        vpu_max_old = vpu_max;
            
                        iter=1;
                        maxIter=1000;
                        foundHC = false;
            
                        while  (vpu_max<=Vpu_UB) && (iter<=maxIter) && ~foundHC %%&& (pv_old<PV_limit*load_kW_All)%&& bb~=12%&& (pv_old<HC_kW_revCurrent(bb,rr))
            
                            %pvnew = pv_old + 100; %Joe's line
                            pvnew = pv_old + 100;
            
                            DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                                pvnew,pvnew);
            
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
            
                            DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 maxcontroliter=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton, tried algorithm=normal and NCIM, did not work and max iteration reached, with direct solution max iteration reached too
                            DSSText.Command = 'solve';
            
                            if ~DSSCircuit.Solution.Converged
                                
                                warnSt = circuitCheck(DSSCircObj);
                                bus_issue(bb,cc)=1;
                                DSSText.Command = 'set mode=direct';
                                DSSText.Command = 'solve';
            
                                %error('Circuit did not converge...\n\n');
                            end
            
                            Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                            vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                            vpu_max = max(vpu_temp);
            
                            if (vpu_max<=Vpu_UB)
                                pv_old = pvnew;
                                vpu_max_old = vpu_max;
                                if (pv_old>PV_limit*load_kW_All)
                                    HC_kW_Vlimit_DERgen(bb,cc) = pv_old;
                                    foundHC = true;
                                end
                            else
                                HC_kW_Vlimit_DERgen(bb,cc) = pv_old;
                                foundHC = true;
                            end
            
                            iter=iter+1;
            
                        end
            
                        if iter > maxIter
                            warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                            HC_kW_Vlimit_DERgen(bb,cc) = pv_old;
                        end
            
                    end
            
                    % reset the pv system size and simulation parameters
                    DSSText.Command = sprintf('Edit PVsystem.PV1 enabled=false kva=%1.5f pmpp=%1.5f',...
                        PV_kW_min,PV_kW_min);
                    for dd=1:length(Scenarios_new)
                        DSSText.Command =sprintf('Edit %s enabled=true',Scenarios_new(dd).(fields{cc}));
                    end
                    DSSText.Command = 'set mode=direct';
                    DSSText.Command = 'solve';
            
                    DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                    DSSText.Command = 'solve';
                end
            else
                DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f kva=%1.5f pmpp=%1.5f pf=1 pfpriority=true balanced=true enabled=false',...
                Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase,PV_kW_min,PV_kW_min);
                % DSSText.Command = sprintf('batchedit PVsystem..* bus1=%s enabled=false kva=%1.5f pmpp=%1.5f',...
                % Bus_LV(bb).name,PV_kW_min,PV_kW_min);
    
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                DSSText.Command = 'solve';
        
                Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                vpu_max_noPV = max(vpu_temp);
        
                DSSText.Command =sprintf('Edit PVsystem.PV1 phases=3 bus1=%s kv=%1.5f pf=1 pfpriority=true balanced=true enabled=true',...
                Bus_LV(bb).name,sqrt(3)*Bus_LV(bb).kVBase);
                %DSSText.Command = sprintf('batchedit PVsystem..* enabled=true');
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton
                DSSText.Command = 'solve';
        
                Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                vpu_max = max(vpu_temp);
        
                if any(vpu_temp>Vpu_UB)
                    HC_kW_Vlimit_DERgen(bb,rr) = 0;
                else
        
                    pv_old = PV_kW_min;
                    vpu_max_old = vpu_max;
        
                    iter=1;
                    maxIter=1000;
                    foundHC = false;
        
                    while  (vpu_max<=Vpu_UB) && (iter<=maxIter) && ~foundHC %&& (pv_old<PV_limit*load_kW_All)%&& bb~=12%&& (pv_old<HC_kW_revCurrent(bb,rr))
        
                        %pvnew = pv_old + 100; %Joe's line
                        pvnew = pv_old + 100;
        
                        DSSText.Command = sprintf('Edit PVsystem.PV1 kva=%1.5f pmpp=%1.5f',...
                            pvnew,pvnew);
        
                        DSSText.Command = 'set mode=direct';
                        DSSText.Command = 'solve';
        
                        DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 maxcontroliter=1000 miniterations=1 tolerance=0.0001'; %algorithm=newton, tried algorithm=normal and NCIM, did not work and max iteration reached, with direct solution max iteration reached too
                        DSSText.Command = 'solve';
        
                        if ~DSSCircuit.Solution.Converged
                            
                            warnSt = circuitCheck(DSSCircObj);
                            bus_issue(bb,rr)=1;
                            DSSText.Command = 'set mode=direct';
                            DSSText.Command = 'solve';
        
                            %error('Circuit did not converge...\n\n');
                        end
        
                        Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
                        vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
                        vpu_max = max(vpu_temp);
        
                        if (vpu_max<=Vpu_UB)
                            pv_old = pvnew;
                            vpu_max_old = vpu_max;
                            if (pv_old>PV_limit*load_kW_All)
                                HC_kW_Vlimit_DERgen(bb,rr) = pv_old;
                                foundHC = true;
                            end
                        else
                            HC_kW_Vlimit_DERgen(bb,rr) = pv_old;
                            foundHC = true;
                        end
        
                        iter=iter+1;
        
                    end
        
                    if iter > maxIter
                        warning('==Max iterations reached rr=%i bb=%i...\n\n',rr,bb);
                        HC_kW_Vlimit_DERgen(bb,rr) = pv_old;
                    end
        
                end
        
                % reset the pv system size and simulation parameters
                DSSText.Command = sprintf('Edit PVsystem.PV1 enabled=false kva=%1.5f pmpp=%1.5f',...
                    PV_kW_min,PV_kW_min);
    
                DSSText.Command = 'set mode=direct';
                DSSText.Command = 'solve';
        
                DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';
                DSSText.Command = 'solve';
            end
    
        end
        close(h);
    end
end


            %                 vpu_max_old = vpu_max_noPV;
            %                 pv_old1 = 0;
            %                 pv_old2 = PV_kW_min;


%                 while  (abs(vpu_max-Vpu_UB)>convergenceThresh) && (iter<=maxIter) && iterSovle
% 
%                     pvnew = interp1([vpu_max_old,vpu_max],[pv_old1,pv_old2],Vpu_UB,'linear','extrap');
% 
%                     DSSText.Command = sprintf('batchedit PVsystem..* kva=%1.5f pmpp=%1.5f',...
%                         pvnew,pvnew);
% 
%                     DSSText.Command = 'set mode=direct';
%                     DSSText.Command = 'solve';
% 
%                     DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 maxcontroliter=100 miniterations=1 tolerance=0.0001'; %algorithm=newton
%                     DSSText.Command = 'solve';
% 
%                     if ~DSSCircuit.Solution.Converged
%                         error('Circuit did not converge...\n\n');
%                     end
% 
%                     vpu_max_old = vpu_max;
%                     pv_old1 = pv_old2;
%                     pv_old2 = pvnew;
% 
%                     Bus_LV_Vpu_now = getBusInfo(DSSCircObj,{Bus_LV.name}');
%                     vpu_temp = [Bus_LV_Vpu_now.voltagePU]';
%                     vpu_max = max(vpu_temp);
% 
%                     iter=iter+1;
% 
% %                            pvNow = getPVInfo(DSSCircObj);
% 
%                 end


%% Disable the PV system and set all loads back to their original values
DSSText.Command = 'batchedit PVsystem..* enabled=false';
DSSText.Command = 'set mode=direct';
DSSText.Command = 'solve';


DSSText.Command = 'set mode=snap algorithm=newton controlmode=static maxiterations=1000 miniterations=1 tolerance=0.0001';%algorithm=newton
DSSText.Command = 'solve';

for ll = 1:length(Load_Original)
    DSSText.Command = sprintf('edit load.%s kw=%1.5f kvar=%1.5f',Load_Original(ll).name, Load_Original(ll).kW, Load_Original(ll).kvar);
end
DSSText.Command = 'set mode=snap controlmode=static';
DSSText.Command = 'solve';

%% plot HCA results on circuit plot (% of load)

cmap = colormap(jet(1000));
cmapCumPct = ((1:1:length(cmap))./length(cmap))';

VHC_all = HC_penetration_revCurrent;
VHC_all = VHC_all(:);
vhcRange_upper = max(VHC_all);
vhcRange_lower = 0;
maxVal = max(abs([vhcRange_upper,vhcRange_lower]));
vhcRange = maxVal;


%% create plot
if ShowPlots && CalcProtection
figure;
if length(revThresh)>1
    for rr=1:length(revThresh)
        subplot(2,ceil(length(revThresh)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA',revThresh(rr)));
        hold on;

        VHC = HC_penetration_revCurrent(:,rr);
        idxC = (VHC)./vhcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (% of Peak Load)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    for rr=1:length(fields)
        subplot(2,ceil(length(fields)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA Scenario %i',revThresh, rr));
        hold on;

        VHC = HC_penetration_revCurrent(:,rr);
        idxC = (VHC)./vhcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (% of Peak Load)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
else
    Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
    ax = gca;
    xticks('');
    yticks('');
    hold on;

    VHC = HC_penetration_revCurrent(:,rr);
    idxC = (VHC)./vhcRange;
    cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

    cbh = colorbar();
    colormap(cmap);
    caxis manual;
    caxis([0,maxVal]);
    ylabel(cbh,'HC (% of Load)');
    for ii = 1:length(Bus_LV)
        plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
    end
    xticklabels('');
    yticklabels('');
    legend();
    set(gca,'FontSize',13,'FontWeight','Bold');
end
sgtitle(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('circuit_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));
%% plot HCA results on circuit plot (kW)
VHC_all = HC_kW_revCurrent/1e3;
VHC_all = VHC_all(:);
vhcRange_upper = max(VHC_all);
vhcRange_lower = 0;
maxVal = max(abs([vhcRange_upper,vhcRange_lower]));
vhcRange = maxVal;

%% Create plot

figure;
if length(revThresh)>1
    for rr=1:length(revThresh)
        subplot(2,ceil(length(revThresh)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA',revThresh(rr)));
        hold on;

        VHC = HC_kW_revCurrent(:,rr)/1e3;
        idxC = (VHC)./vhcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (MW)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    for rr=1:length(fields)
        subplot(2,ceil(length(fields)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA Scenario %i',revThresh, rr));
        hold on;

        VHC = HC_kW_revCurrent(:,rr)/1e3;
        idxC = (VHC)./vhcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (MW)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
else
    Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
    ax = gca;
    xticks('');
    yticks('');
    hold on;

    VHC = HC_kW_revCurrent(:,rr)/1e3;
    idxC = (VHC)./vhcRange;
    cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

    cbh = colorbar();
    colormap(cmap);
    caxis manual;
    caxis([0,maxVal]);
    ylabel(cbh,'HC (MW)');
    for ii = 1:length(Bus_LV)
        plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'-ko','MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
    end
    xticklabels('');
    yticklabels('');
    legend();
    set(gca,'FontSize',13,'FontWeight','Bold');
end
sgtitle(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('circuit_with_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));
% 
%% boxplots of HC results (MW)


if length(revThresh)>=1 && ~istable(Scenarios)
    % labels = arrayfun(@(x) sprintf('%d', x), 1:length(revThresh), 'UniformOutput', false);
    labels=revThresh;
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    labels = arrayfun(@(x) sprintf('Scenario %d', x), 1:length(fields), 'UniformOutput', false);
end

figure;

boxplot(HC_kW_revCurrent/1e3,'Labels', labels);
ylabel('HC (MW)');
%xticklabels(revThresh);
xlabel('Reverse Current Threshold (A)');
set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
ytickformat('%1.1f');

title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('boxplot_with_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

%% boxplots of HC results (%)
figure;
boxplot(HC_penetration_revCurrent,'Labels', labels);
ylabel('HC (% of Load)');
%xticklabels(revThresh);
xlabel('Reverse Current Desensitization (A)');
set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
ytickformat('%1.1f');

title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('boxplot_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

end
%% Solve the HCA - min of all
if istable(Scenarios)

    HC_kW_total = nan(length(Bus_LV),width(Scenarios));
    HC_penetration_pct_total = nan(length(Bus_LV),width(Scenarios));

    HC_kW_minall = nan(length(Bus_LV),1);
    HC_penetration_pct_minall = nan(length(Bus_LV),1);
       
    
else
    HC_kW_total = nan(length(Bus_LV),length(revThresh));
    HC_penetration_pct_total = nan(length(Bus_LV),length(revThresh));
   
    
end
% end

%% Create Results Structure

Results.SimSettings.revThresh = revThresh;
Results.SimSettings.loadLevel = loadLevel;
Results.BusInfo = Bus_LV;
Results.HC_kW_protection = HC_kW_revCurrent;
Results.Buserror = bus_issue;
Results.HC_kW_voltage = HC_kW_Vlimit_DERgen;
Results.HC_kW_thermal = HC_kW_Thermal_DERgen;
if istable(Scenarios)
    no_of_scenarios=width(Scenarios);
    HC_kW_minall=min([min(Results.HC_kW_protection,[],2),min(Results.HC_kW_voltage,[],2),min(Results.HC_kW_thermal,[],2)],[],2);
    HC_kW_minbytype=[min(Results.HC_kW_protection,[],2),min(Results.HC_kW_voltage,[],2),min(Results.HC_kW_thermal,[],2)];
    Results.HC_kW_minall=HC_kW_minall;
    Results.HC_penetration_pct_minall=(Results.HC_kW_minall/load_kW_All)*100;
    HC_penetration_pct_minall=Results.HC_penetration_pct_minall;
    
    for i=1:no_of_scenarios

        Results.HC_kW_total(:,i) = min([HC_kW_revCurrent(:,i),HC_kW_Vlimit_DERgen(:,i),HC_kW_Thermal_DERgen(:,i)],[],2);
        HC_kW_total(:,i)=Results.HC_kW_total(:,i);
        
    end
elseif length(revThresh)>1
    for i=1:length(revThresh)
        Results.HC_kW_total(:,i) = min([HC_kW_revCurrent(:,i),HC_kW_Vlimit_DERgen(:,i),HC_kW_Thermal_DERgen(:,i)],[],2);
        HC_kW_total(:,i)=Results.HC_kW_total(:,i);
    end
else
    Results.HC_kW_total = min([HC_kW_revCurrent,HC_kW_Vlimit_DERgen,HC_kW_Thermal_DERgen],[],2);
    HC_kW_total=Results.HC_kW_total;
end
Results.HC_penetration_pct_total_load = (Results.HC_kW_total/load_kW_All)*100;

HC_penetration_pct_total=Results.HC_penetration_pct_total_load;
Results.CircuitName = DSSCircuit.Name;



%%%plot for min of all
cmap = colormap(jet(1000));
cmapCumPct = ((1:1:length(cmap))./length(cmap))';

HC_all = HC_penetration_pct_total;
HC_all = HC_all(:);
hcRange_upper = max(HC_all);
hcRange_lower = 0;
maxVal = max(abs([hcRange_upper,hcRange_lower]));
hcRange = maxVal;


%% create plot
if ShowPlots && CalcVoltage && CalcThermal && CalcProtection
figure;
if length(revThresh)>1
    for rr=1:length(revThresh)
        subplot(2,ceil(length(revThresh)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA',revThresh(rr)));
        hold on;
        data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
        markers = {'o', 's', 'd'};
        [min_values, min_indices] = min(data_all, [], 2);
        
        HC = HC_penetration_pct_total(:,rr);
        idxC = (HC)./hcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (% of Peak Load)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    for rr=1:length(fields)
        subplot(2,ceil(length(fields)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA Scenario %i',revThresh, rr));
        hold on;

        data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
        markers = {'o', 's', 'd'};
        [min_values, min_indices] = min(data_all, [], 2);

        HC = HC_penetration_pct_total(:,rr);
        idxC = (HC)./hcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (% of Peak Load)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
else
    Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
    ax = gca;
    xticks('');
    yticks('');
    hold on;
    data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
    markers = {'o', 's', 'd'};
    [min_values, min_indices] = min(data_all, [], 2);
    HC = HC_penetration_pct_total(:,rr);
    idxC = (HC)./hcRange;
    cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

    cbh = colorbar();
    colormap(cmap);
    caxis manual;
    caxis([0,maxVal]);
    ylabel(cbh,'HC (% of Load)');
    for ii = 1:length(Bus_LV)
        plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
    end
    xticklabels('');
    yticklabels('');
    legend();
    set(gca,'FontSize',13,'FontWeight','Bold');
end
sgtitle(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('mincircuit_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));
%% plot HCA results on circuit plot (kW)
HC_all = HC_kW_total/1e3;
HC_all = HC_all(:);
hcRange_upper = max(HC_all);
hcRange_lower = 0;
maxVal = max(abs([hcRange_upper,hcRange_lower]));
hcRange = maxVal;

%% Create plot

figure;
if length(revThresh)>1
    for rr=1:length(revThresh)
        subplot(2,ceil(length(revThresh)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA',revThresh(rr)));
        hold on;
        data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
        markers = {'o', 's', 'd'};
        [min_values, min_indices] = min(data_all, [], 2);

        HC = HC_kW_total(:,rr)/1e3;
        idxC = (HC)./hcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (MW)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    for rr=1:length(fields)
        subplot(2,ceil(length(fields)/2),rr);
        Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
        ax = gca;
        xticks('');
        yticks('');
        title(sprintf('Desensitized %iA Scenario %i',revThresh, rr));
        hold on;
        data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
        markers = {'o', 's', 'd'};
        [min_values, min_indices] = min(data_all, [], 2);

        HC = HC_kW_total(:,rr)/1e3;
        idxC = (HC)./hcRange;
        cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

        cbh = colorbar();
        colormap(cmap);
        caxis manual;
        caxis([0,maxVal]);
        ylabel(cbh,'HC (MW)');
        for ii = 1:length(Bus_LV)
            plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
        end
        xticklabels('');
        yticklabels('');
        legend();
        set(gca,'FontSize',13,'FontWeight','Bold');
    end
else
    Handles2 = plotCircuitLines(DSSCircObj,'SwitchMarker','on','ServiceTransformerMarker','on','MVTransformerMarker','on'); %,'Coloring','voltagePU');
    ax = gca;
    xticks('');
    yticks('');
    hold on;
    data_all=[HC_kW_revCurrent(:,rr),HC_kW_Vlimit_DERgen(:,rr),HC_kW_Thermal_DERgen(:,rr)];
    markers = {'o', 's', 'd'};
    [min_values, min_indices] = min(data_all, [], 2);

    HC = HC_kW_total(:,rr)/1e3;
    idxC = (HC)./hcRange;
    cDataPlot = interp1(cmapCumPct,cmap,idxC','linear','extrap');

    cbh = colorbar();
    colormap(cmap);
    caxis manual;
    caxis([0,maxVal]);
    ylabel(cbh,'HC (MW)');
    for ii = 1:length(Bus_LV)
        plot(ax,Bus_LV(ii).coordinates(2),Bus_LV(ii).coordinates(1),'Marker', markers{min_indices (ii)},'MarkerSize',7.5,'MarkerFaceColor',cDataPlot(ii,:),'LineStyle','none','HandleVisibility','off');
    end
    xticklabels('');
    yticklabels('');
    legend();
    set(gca,'FontSize',13,'FontWeight','Bold');
end
sgtitle(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('mincircuit_with_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));
% 
%% boxplots of HC results (MW)
if length(revThresh)>=1 && ~istable(Scenarios)
    %labels = arrayfun(@(x) sprintf('%d', x), 1:length(revThresh), 'UniformOutput', false);
    labels=revThresh;
elseif istable(Scenarios)
    Scenarios_new=table2struct(Scenarios);
    fields=fieldnames (Scenarios_new);
    labels = arrayfun(@(x) sprintf('Scenario %d', x), 1:length(fields), 'UniformOutput', false);
end
figure;

boxplot(HC_kW_total/1e3,'Labels', labels);
ylabel('HC (MW)');
%xticklabels(revThresh);
xlabel('Reverse Current Threshold (A)');
set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
ytickformat('%1.1f');

title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('minboxplot_with_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

figure;
boxplot(HC_penetration_pct_total,'Labels', labels);
ylabel('HC (% of Load)');
%xticklabels(revThresh);
xlabel('Reverse Current Desensitization (A)');
set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
ytickformat('%1.1f');

title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
saveas(gcf, sprintf('minboxplot_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

if istable(Scenarios)
    figure;

    boxplot(HC_kW_minall/1e3,'Labels', revThresh);
    ylabel('HC (MW)');
    %xticklabels(revThresh);
    xlabel('Min HC of all contingencies');
    set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
    ytickformat('%1.1f');

    title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
    saveas(gcf, sprintf('minallboxplot_with_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

    figure;
    boxplot(HC_penetration_pct_minall,'Labels', revThresh);
    ylabel('HC (% of Load)');
    %xticklabels(revThresh);
    xlabel('Min HC of all contingencies');
    set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
    ytickformat('%1.1f');

    title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
    saveas(gcf, sprintf('minallboxplot_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));

    figure;
    boxplot(HC_kW_minbytype,'Labels',{'protection','voltage','thermal'} );
    ylabel('HC (MW)');
    %xticklabels(revThresh);
    xlabel('Min HC of all contingencies');
    set(gca,'YGrid','on','GridAlpha',0.2,'FontSize',12,'FontWeight','Bold');
    ytickformat('%1.1f');

    title(sprintf('Loading Scenario: %1.1f%%',loadLevel));
    saveas(gcf, sprintf('mintypeboxplot_withpercent_%s%s.jpeg',num2str(loadLevel),DSSCircuit.Name));
end
end


































