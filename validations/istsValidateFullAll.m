function istsValidateFullAll(varargin)
% Run unit tests for SST
%
% Syntax:
%    istsValidateFullAll
%
% Description:
%    Run the set of unit tests for ISETTreeShrew. Return
%    true or false depending on success. Designed for use with Jenkins
%    build integration.
%

%% Close all figures so that we start with a clean slate
close all;

%% Set return flag
success = true;

%% We will use preferences for this project
thisProject = 'ISETTreeShrew';

%% Run your test here.  Set success to false if any fail
%
% For example ...
% status = ReceptorIsolateDemo('validate','basichuman','whichDirectionNumber',1);
% if (~status)
%     success = false;
% end

%% Report whether we are OK
assert(success, 'One or more validations failed.');

end