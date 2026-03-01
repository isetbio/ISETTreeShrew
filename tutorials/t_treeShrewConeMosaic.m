function t_treeShrewConeMosaic()
% Generate treeshrew cone mosaic
%
% Description:
%   Generates a treeshrew cone mosaic

% History:
%    02/28/26  NPC  Wrote it.

  
    theConeMosaic = cMosaicTreeShrewCreate(...
        'fovDegs', [5 1], ...
        'integrationTimeSeconds', 100/1000);

    theConeMosaic.visualize();
    theConeMosaic.visualize( ...
        'domain', 'microns');
end


