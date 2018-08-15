function out = centerPos(SP,Lx,Lz)
% CENTERPOS Center the provided positions (SP.x, SP.z) such that
% SP.x \in [0, Lx] -> out.x \in [-0.5,0.5]*Lx and equivalently for .z/Lz

out = SP;

out.x = out.x - Lx/2;
out.z = out.z - Lz/2;

end