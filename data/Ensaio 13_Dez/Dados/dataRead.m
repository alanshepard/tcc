%% Programmer: Frederico Bolsoni Oliveira
% Reads the VT-80 test acquired data.
clear vars;
close all;
clc;

% Name structure of the data files:
dataNameStructure = {'TempCelCarga'
                     'TempHumAmb'
                     'TempPresAmb'
                     'TempPresDif'
                     'TempTempAmb'
                     'TempTermopar'};
dataNameLengths = zeros(length(dataNameStructure), 1);

rootDataFolderName = '.';

subfolderNames = {'10'
                  '20'
                  '30'
                  '40'
                  '50'
                  '60'
                  '80'
                  '100'
                  'Idle'
                  'Startup'
                  'Zero'};

dataExtension = '.txt';

emptyCell = cell(length(dataNameStructure), 1);
for i = 1:length(emptyCell)
    emptyCell{i} = [0 0];
end

% Analyzes the file structure:
rootDir = dir(rootDataFolderName);
numFileDir = length(subfolderNames);

% Searches the file structure for data:
acquiredData = cell(numFileDir, 1);
for i = 1:numFileDir
    tempCell = zeros(length(subfolderNames), 1);
    
    % Changes the current directory:
    currentDirName = strcat(rootDataFolderName, '/', subfolderNames{i});
    currentDir = dir(currentDirName);
    
    % Checks if the current directory exists:
    if isempty(currentDir)
        userInput = '';
        while ~strcmpi(userInput, 'y') && ~strcmpi(userInput, 'n')
            userInput = input(sprintf('The folder ''%s'' was not found. Do you want to skip it (y/n)?\n', subfolderNames{i}), 's');
        end
        if strcmpi(userInput, 'n')
            fprintf('The script has been stopped by the user.\n')
            return;
        else
            acquiredData{i} = emptyCell;
            continue;
        end
    end
    dataCell = emptyCell;
    
    % Opens the current file:
    for j = 1:length(dataNameStructure)
        fullSubfolderPath = strcat(rootDataFolderName, '\', subfolderNames{i});
        currentDir = dir(fullSubfolderPath);
        dataNameGoal = dataNameStructure{j};
        dataNameLength = length(dataNameGoal);
        foundData = false;
        
        
        % Searches for the correct data name:
        for k = 1:length(currentDir)
            currentDataName = currentDir(k).name;
            if length(currentDataName) < dataNameLength
                continue;
            end
            
            if strcmp(dataNameGoal, currentDataName(1:dataNameLength))
                % Checks if the file extension is correct:
                if strcmp(dataExtension, currentDataName(end-length(dataExtension)+1:end))
                    foundData = true;
                    dataFullPath = strcat(fullSubfolderPath, '\', currentDataName);
                    tempData = importdata(dataFullPath);
                    
                    % Creates the time vector:
                    dt = 1/tempData(1);
                    dataLength = length(tempData)-1;
                    t = (0:dt:(dataLength-1)*dt)';
                    
                    % Stores the desired data:
                    dataCell{j} = [t tempData(2:end)];
                    break;
                end
            end
        end
        
        % Warns the user if the data was not found:
        if ~foundData
            userInput = '';
            while ~strcmpi(userInput, 'y') && ~strcmpi(userInput, 'n')
                userInput = input(sprintf('The file ''%s'' in folder ''%s''was not found. Do you want to skip it (y/n)?\n', dataNameGoal, dataFullPath), 's');
            end
            if strcmpi(userInput, 'n')
                fprintf('The script has been stopped by the user.\n')
                return;
            else
                dataCell{j} = [0 0];
            end
        end
    end
    
    % Saves the current data in a cell:
    acquiredData{i} = dataCell;
    
end
clear tempData;
