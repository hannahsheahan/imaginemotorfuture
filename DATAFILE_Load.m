function [D,R] = DATAFILE_Load(FileName)
% *************************************************************************
% function [D,R] = DATAFILE_Load(FileName)
%
% Reads experimental data for trials from a file and saves as a MAT file
%
%   [D,R] = DATAFILE_Load(FileName)
%
%   D = Data
%   R = Return code (Boolean success)
%
% Author:   James Ingram
% Email:    jni20@cam.ac.uk
% *************************************************************************

    fmat = sprintf('%s.MAT',FileName);
    fdat = sprintf('%s.DAT',FileName);
    
    if( exist(fmat) == 2 )
        disp(sprintf('Load: %s',fmat));
        x = load(fmat);
        D = x.D;
        clear x;
        R = 1;
    elseif( exist(fdat) == 2 )
        [D,R] = DATAFILE_Read(fdat);
        if( R )
            disp(sprintf('Save: %s',fmat));
            save(fmat,'D', '-v7.3');
        end
    else
        error('File not found: %s',FileName);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
