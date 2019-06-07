% ISETTreeShewLocalHookTemplate
%
% Template for setting preferences and other configuration things, for the
% ISETTreeShrew project.

% 10/23/18  dhb   Wrote it.

%% Define project
projectName = 'ISETTreeShrew';

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end

%% Specify project location
projectBaseDir = tbLocateProject(projectName);

% Figure out where baseDir for other kinds of data files is.
%
% Can only do this when we have GetComputerInfo available.
if (exist('GetComputerInfo','file'))
    sysInfo = GetComputerInfo();
    switch (sysInfo.localHostName)
        case 'eagleray'
            % DHB's desktop
            baseDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
            
        case {'Manta', 'Manta-2'}
            % Nicolas's iMac
            baseDir = fullfile(filesep,'Volumes','DropBoxDisk/Dropbox','Dropbox (Aguirre-Brainard Lab)');
            
        otherwise
            % Some unspecified machine, try user specific customization
            switch(sysInfo.userShortName)
                % Could put user specific things in, but at the moment generic
                % is good enough.
                otherwise
                    baseDir = ['/Users/' sysInfo.userShortName '/Dropbox (Aguirre-Brainard Lab)'];
            end
    end
end

% Get path to data in project code with getpref('ISETTwoLine','dataDir');
setpref(projectName,'dataDir',fullfile(baseDir,'IBIO_Analysis',projectName));




