% Data Analysis 
clc; clear; 

load("Data/ExpData.mat");


%% Figure 1D: Creation 

Fig1Scenarios = [ "IKK_RelA_negativeControl";
                  "IKK_RelA_Single_30sec_pulse";
                  "IKK_RelA_Single_2min_pulse";
                  "IKK_RelA_Single_6min_pulse";
                  "IKK_RelA_Single_15min_pulse";
                  "IKK_RelA_Single_30min_pulse"];

Fig1Labels = ["Control","30 Second Pulse","2 Minute Pulse", ...
              "6 Minute Pulse","15 Minute Pulse","30 Minute Pulse"];

lengthOfPulse = {0 0.5 2 6 15 30}; 

plotMultiPanelPulseResponse(Fig1Scenarios, Fig1Labels, lengthOfPulse)
%exportPulseDataToExcel(Fig1Scenarios, Fig1Labels,"Figure1D_Data.xlsx")

%% Figure 1E,1F 
Fig1Labels = ["0/0","1/30 Sec","1/2Min","1/6min","1/15min","1/30min"]
plotAUC_Multipanel(Fig1Scenarios, Fig1Labels, lengthOfPulse,"Figure1E_F_Data.xlsx")


%% Figure 2A
Fig2Scenarios = ["IKK_RelA_Single_6min_pulse";
                "IKK_RelA_Two_3min_pulses_05min_gap"; 
                "IKK_RelA_Three_2min_pulses_05min_gap";
                "IKK_RelA_Four_90sec_pulses_05min_gap";
                "IKK_RelA_Single_15min_pulse"];

Fig2Labels = ["1X6 min","2X3 min","3X2 min","4X1.5 min","1X15 min"];

lengthOfPulse2 = {6; 
                 struct('nPulses', 2, 'pulseDuration', 3, 'gap', 5);
                 struct('nPulses', 3, 'pulseDuration', 2, 'gap', 5);
                 struct('nPulses', 4, 'pulseDuration', 1.5, 'gap', 5);
                 15}; 

plotMultiPanelPulseResponse(Fig2Scenarios, Fig2Labels, lengthOfPulse2)
%exportPulseDataToExcel(Fig2Scenarios, Fig2Labels,"Figure2A_Data.xlsx")

%% Figure 2B,2C
Fig2Labels = ["1/6min","2/6min","3/6min","4/6min","1/15min"];
plotAUC_Multipanel(Fig2Scenarios, Fig2Labels, lengthOfPulse2,"Figure2B_C_Data.xlsx")





%%
function exportPulseDataToExcel(FigScenarios, FigLabels, filename)
    % exportPulseDataToExcel - Exports multi-cell response data to Excel
    %
    % Inputs:
    %   FigScenarios - string array of field names in ExpData
    %   FigLabels    - labels for column headers in Excel
    %   filename     - filename of the Excel file to write (e.g., 'output.xlsx')

    allData = evalin('base', 'ExpData');
    output = []; 
    labels = []; 
    scenarioLabel = [""];
    %% First, determine the maximum number of timepoints
    for i = 1:length(FigScenarios)
        ithExp = FigScenarios{i}; 
        ithLabel = FigLabels{i}; 
        ithIKK = allData.(ithExp).IKK_spotNumber; 
        ithRelAFC = allData.(ithExp).RelA_FC; 
        
        nCells = size(ithIKK,2);
        nTimepoints = size(ithIKK,1);
        for j = 1:nCells
            output = [output ithIKK(:,j)];
            output = [output ithRelAFC(:,j)];
            labels = [labels "NEMO Spots (#)" "Nuclear RelA (Fold Change)"];
            jthLabel = ithLabel+" Cell #"+j;
            scenarioLabel = [scenarioLabel  jthLabel jthLabel ]; 
        end 
    end 
    time = 0:nTimepoints-1; 
    time = time*4; 
    output = [time' output]; 
    labels = ["Time [Minutes]" labels];
    %%
    output=num2cell(output);
    output = [num2cell(labels); output];
    output = [num2cell(scenarioLabel); output];
    % Write to Excel
%writecell(output, filename);
end



%% Functions 
function plotAUC_Multipanel(FigScenarios, FigLabels, lengthOfPulse,filename)
% plotAUC_Multipanel - Generates a 1x3 layout for AUC NEMO, RelA, and their correlation
% Inputs:
%   FigScenarios   - string array of scenario names (e.g., data field names)
%   FigLabels      - base labels for each group (e.g., '2 min', etc)
%   lengthOfPulse  - cell array; scalar for single pulse, struct for multi-pulse

% Load data
allData = evalin('base', 'ExpData');
nGroups = length(FigScenarios);
AUC_NEMO = [];
AUC_RelA = [];
labels = [];
displayLabels = []; 
cellNumber = [];
pulseInfo = struct('isMultiPulse', false(1,nGroups));
updatedLabels = strings(nGroups,1);

% Store scenario-specific AUC for correlation filtering
allAUC = struct('nemo', [], 'rela', [], 'group', []);

% Compute AUC and format labels
for i = 1:nGroups
    data = allData.(FigScenarios(i));
    dt = 4;
    if isfield(data, 'IKK_spotNumber')
        auc_nemo = trapz(data.IKK_spotNumber, 1);
        AUC_NEMO = [AUC_NEMO, auc_nemo];
    end
    if isfield(data, 'RelA_FC')
        ithRelA = data.RelA_FC-1;
        ithRelA(ithRelA < 0) = 0;
        auc_rela = trapz(ithRelA, 1) ;
        AUC_RelA = [AUC_RelA, auc_rela];
    end
    groupSize = size(data.IKK_spotNumber, 2);
    labels = [labels, repmat(i, 1, groupSize)];
    displayLabels = [displayLabels, repmat(FigLabels(i), 1, groupSize)];
    cellNumber = [cellNumber, 1:groupSize];

    pulseInput = lengthOfPulse{i};
    if isnumeric(pulseInput)
        updatedLabels(i) = FigLabels(i);
    elseif isstruct(pulseInput)
        totalDur = pulseInput.nPulses * pulseInput.pulseDuration;
        updatedLabels(i) = sprintf('%d/%d', pulseInput.nPulses, round(totalDur));
        pulseInfo(i).isMultiPulse = true;
    end

    % Save AUCs for correlation panel with group index
    allAUC(i).nemo = auc_nemo;
    allAUC(i).rela = auc_rela;
    allAUC(i).group = FigScenarios(i);
end

% Create tiled layout
figure('Position', [100 100 1000 400])
tiledlayout(1,3,'TileSpacing','compact')

% Panel 1: AUC NEMO
nexttile(1)
boxplot(AUC_NEMO, labels, 'Colors', [0.2 0.6 0.2])
set(gca, 'XTickLabel', updatedLabels)
xlabel({'Pulses (#)', '/ Total IL-1', 'duration (min)'})
ylabel('AUC NEMO spots (#*min)')
title('B')
xlim padded

% Panel 2: AUC RelA
nexttile(2)
boxplot(AUC_RelA, labels, 'Colors', [0.9 0.4 0.3])
set(gca, 'XTickLabel', updatedLabels)
xlabel({'Pulses (#)', '/ Total IL-1', 'duration (min)'})
ylabel('AUC Nuclear RelA (fold change)')
xlim padded

% Panel 3: Scatter (exclude negative control)
nexttile(3)
excludeIdx = strcmp(FigScenarios, 'IKK_RelA_negativeControl');
AUC_NEMO_filt = [];
AUC_RelA_filt = [];
for i = 1:nGroups
    if ~excludeIdx(i)
        AUC_NEMO_filt = [AUC_NEMO_filt, allAUC(i).nemo];
        AUC_RelA_filt = [AUC_RelA_filt, allAUC(i).rela];
    end
end
scatter(AUC_NEMO_filt, AUC_RelA_filt, 25, 'filled', 'MarkerFaceColor', [0.7 0.6 0.2])
xlabel('AUC NEMO spots (#*min)')
ylabel('AUC Nuclear RelA (fold change)')
[rho, ~] = corr(AUC_NEMO_filt', AUC_RelA_filt', 'Type', 'Spearman');
title(['\rho = ' num2str(round(rho,3))])

% Save Output 
output = [cellNumber' AUC_NEMO' AUC_RelA']; 
output = num2cell(output);
output = [num2cell(displayLabels') output]
topRow = ["Pulses (#) / Total IL-1 Duration (min)","Cell Number", "AUC NEMO spots (#*min)","AUC Nuclear RelA (fold change)"]
output = [num2cell(topRow);output];
%writecell(output, filename);
end




%%
function plotMultiPanelPulseResponse(FigScenarios, FigLabels, lengthOfPulse)
    % plotMultiPanelPulseResponse - Plots multi-panel response figure
    %
    % Inputs:
    %   FigScenarios   - string array of field names in ExpData
    %   FigLabels      - labels for the top row of each column
    %   lengthOfPulse  - cell array; each entry is either:
    %                    - scalar (single pulse duration), or
    %                    - struct with fields: nPulses, pulseDuration, gap

    % Access data structure from base workspace
    figure
    allData = evalin('base', 'ExpData');

    nColumns = length(FigScenarios);
    nRows = 3;

    tiledlayout(nRows, nColumns, 'TileSpacing', 'tight', 'Padding', 'compact')

    for ith = 1:nColumns
        scenario = FigScenarios(ith);
        label = FigLabels(ith);
        pulseInput = lengthOfPulse{ith};

        data = allData.(scenario);

        %% Handle time vector based on size of data
        if isfield(data, 'RelA_FC')
            nTimepoints = size(data.RelA_FC, 2);
        elseif isfield(data, 'IKK_spotNumber')
            nTimepoints = size(data.IKK_spotNumber, 2);
        else
            nTimepoints = 46; % fallback
        end
        timeVec = 0:4:(nTimepoints - 1) * 4;

        %% Row 1: IL-1 pulses
        nexttile(ith)
        hold on
        if isnumeric(pulseInput) && isscalar(pulseInput)
            % Single pulse from 0 to duration
            plot([0 pulseInput], [10 10], 'k', 'LineWidth', 3)
            plot([pulseInput pulseInput], [10 0], 'k', 'LineWidth', 1)
        elseif isstruct(pulseInput)
            % Multi-pulse from pulse spec
            dur = pulseInput.pulseDuration;
            gap = pulseInput.gap;
            nPulses = pulseInput.nPulses;
            startTime = 0;
            for p = 1:nPulses
                t1 = startTime;
                t2 = t1 + dur;
                plot([t1 t2], [10 10], 'k', 'LineWidth', 3)
                plot([t2 t2], [10 0], 'k', 'LineWidth', 1)
                startTime = t2 + gap;
            end
        end
        xlim([0 180])
        ylim([0 12])
        set(gca, 'XTick', [])
        if ith == 1
            ylabel('[IL-1] (ng/mL)')
        else
            set(gca, 'YTickLabel', [])
        end
        title(label, 'FontWeight', 'normal')

        %% Row 2: NEMO spot numbers
        nexttile(nColumns + ith)
        if isfield(data, 'IKK_spotNumber')
            cellData = data.IKK_spotNumber;
            timeVec = 0:4:(size(cellData, 1)-1)*4;
            hold on
            plot(timeVec, cellData', 'Color', [0.3 0.8 0.3 0.2])
            plot(timeVec, mean(cellData, 2), 'Color', [0 0.5 0], 'LineWidth', 2)
            hold off
            ylim([0 500])
        end
        if ith == 1
            ylabel('NEMO spots (#)')
        else
            set(gca, 'YTickLabel', [])
        end

        %% Row 3: Nuclear RelA
        nexttile(2 * nColumns + ith)
        if isfield(data, 'RelA_FC')
            cellData = data.RelA_FC;
            timeVec = 0:4:(size(cellData, 1)-1)*4;
            hold on
            plot(timeVec, movmean(cellData',1), 'Color', [1 0.4 0.2 0.2])
            plot(timeVec, mean(cellData, 2), 'Color', [0.8 0 0], 'LineWidth', 2)
            hold off
            ylim([0 5])
        end
        if ith == 1
            ylabel('Nuclear RelA (fold change)')
            xlabel('Time (min)')
        else
            set(gca, 'YTickLabel', [])
            xlabel('Time (min)')
        end
    end
end