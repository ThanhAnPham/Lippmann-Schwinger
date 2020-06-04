function out = hasandis(par,str)
%% HASANDIS Verify whether structure 'par' has the boolean field 'str' and
%           if yes, return its value. If not, send false.

out = isfield(par,str);
if out
    out = par.(str);
end

end

