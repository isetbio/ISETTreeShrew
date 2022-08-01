function rootPath = isetTreeShrewRootPath()
% Return the path to the root ISETTreeShrew directory
%
% Syntax:
%   rootPath = iisetTreeShrewRootPath();
%
% Description:
%    This points at the top level of the ISETTreeShrew tree on the Matlab path.
%
%    Examples are included within the code.
%
% Inputs:
%    None.
%
% Outputs:
%    rootPath - The root directory for ISETTreeShrew
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - This function works by using the function mfilename to
%      find itself, and then walks back up the result to the top level of
%      isetbio. Thus, you can't move this function within the tree
%      without also adjusting the number of levels in the walk to match
%      where you move it to.]
% 
% See Also:
%    isetbioDataPath, isetRootPath

% History:
%    07/27/17  dhb  Wrote an ISETTreeShew version

% Examples:
%{
    fullfile(isetTreeShrewRootPath, 'data')
%}

%% Get path to this function and then walk back up to the isetbio root.
pathToMe = mfilename('fullpath');

%% Walk back up the chain
rootPath = fileparts(fileparts(fileparts(pathToMe)));