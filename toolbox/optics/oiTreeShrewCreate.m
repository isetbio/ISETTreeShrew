function oi = oiTreeShrewCreate(varargin)
    
    oi.type = 'opticalimage';
    oi = oiSet(oi, 'optics', opticsTreeShrewCreate(varargin{:}));
    oi = oiSet(oi, 'name', 'tree-shrew');
    
    oi = oiSet(oi, 'bit depth', 32);
    oi = oiSet(oi, 'diffuser method', 'skip');
    oi = oiSet(oi, 'consistency', 1);
    
    if checkfields(oi.optics, 'transmittance')
        oi.optics = rmfield(oi.optics, 'transmittance');
    end
        
end
