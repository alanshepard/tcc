%% Programmer: Frederico Bolsoni Oliveira
% Post-processes the data acquired in the 13/12/17 VT80 test.
dataRead;
close all;

%% Plots the raw data:
titlesCell = {'Uninstalled Thrust F'
              'Ambient Humidity'
              'Ambient Pressure P_0'
              'Differential Pressure P_3'
              'Ambient Temperature T_0'
              'Temperature T_9'};
legendLabels = {'10%'
                '20%'
                '30%'
                '40%'
                '50%'
                '60%'
                '80%'
                '100%'
                'Idle'};

colorRGB = {[0 0.4470 0.7410]
            [0.8500 0.3250 0.0980]
            [0.9290 0.6940 0.1250]
            [0.4940 0.1840 0.5560]
            [0.4660 0.6740 0.1880]
            [0.3010 0.7450 0.9330]
            [0.6350 0.0780 0.1840]
            [0.0000 0.5000 0.0000]
            [0.7500 0.0000 0.7500]};   

n = length(titlesCell);
figureHandles = cell(n, 1);
posConfig = [1 1 800 600];
tStudentMultiplier = 2.576;

%% Plots the properties data:
for i = 1:n
    figureHandles{i} = figure('position', posConfig);
    hold on;
    for j = 1:9
        plot(acquiredData{j}{i}(:, 1), acquiredData{j}{i}(:, 2), 'color', colorRGB{j});
    end
    grid;
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title(titlesCell{i});
    legend(legendLabels, 'location', 'best');
end

%% Plots the zero signal:
fZero = figure('position', posConfig);
hold on;
for i = 1:length(titlesCell)
    plot(acquiredData{end}{i}(:, 1), acquiredData{end}{i}(:, 2), 'color', colorRGB{i});
end
grid;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Zero Signal');
legend(titlesCell, 'location', 'best');

%% Plots the startup signal:
fStartup = figure('position', posConfig);
hold on;
for i = 1:length(titlesCell)
    plot(acquiredData{end-1}{i}(:, 1), acquiredData{end-1}{i}(:, 2), 'color', colorRGB{i});
end
grid;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Startup Signal');
legend(titlesCell, 'location', 'best');

%% Calculates the statistical parameters of the data:
meanData = cell(length(acquiredData), 1);
desvPadData = cell(length(acquiredData), 1);
for i = 1:length(acquiredData)
    desvPadAUX = zeros(n, 1);
    meanAUX = zeros(n, 1);
    for j = 1:n
        
        if j == 6
            meanAUX(j) = mean(acquiredData{i}{j}(end-200:end, 2));
            desvPadAUX(j) = std(acquiredData{i}{j}(end-200:end, 2));
        else
            meanAUX(j) = mean(acquiredData{i}{j}(:, 2));
            desvPadAUX(j) = std(acquiredData{i}{j}(:, 2));
        end
        
    end
    meanData{i} = meanAUX;
    desvPadData{i} = desvPadAUX.*tStudentMultiplier;
end
T_Mean = table(meanData{1}, meanData{2}, meanData{3}, meanData{4}, meanData{5}, ...
    meanData{6}, meanData{7}, meanData{8}, meanData{9}, 'RowNames', titlesCell);
display(T_Mean);
T_DesvPad = table(desvPadData{1}, desvPadData{2}, desvPadData{3}, desvPadData{4}, desvPadData{5}, ...
    desvPadData{6}, desvPadData{7}, desvPadData{8}, desvPadData{9}, 'RowNames', titlesCell);
display(T_DesvPad);

%% Calculates the statistical parameters of the zero refference:
meanZero = zeros(n, 1);
desvPadZero = zeros(n, 1);
for i = 1:n
    meanZero(i) = mean(acquiredData{end}{i}(:, 2));
    desvPadZero(i) = std(acquiredData{end}{i}(:, 2));
end
T_Zero = table(meanZero, desvPadZero, 'RowNames', titlesCell);
display(T_Zero);


