function out = createCostsArray(H,y,weights_opt,varargin)
%CREATECOSTSARRAY Create an array of L2 costs with H a cell array with
%a map in each cell and y the data with the third dimension equal to the
%array cell (ad hoc solution)

Ncosts = H.numMaps;

out = cell(Ncosts,1);

w = 1;

repel = repmat({':'},1,max(ndims(y)-1,1));

for kk = 1:Ncosts
    if (islogical(weights_opt) && weights_opt) || isnumeric(weights_opt)
        if ndims(weights_opt)==ndims(y) && all(size(weights_opt)==size(y))
            w = LinOpDiag([],weights_opt(repel{:},kk));
        elseif nargin > 3
            w = abs(y(repel{:},kk));
            w = 1/sum(w(:));
        elseif islogical(weights_opt)
            w = abs(y(repel{:},kk)).^2;
            w = 1/sum(w(:));
        else
            w = weights_opt;
        end
    end
    if nargin > 3
        out{kk} = myCostIntensity3D(y(repel{:},kk),w,varargin{1})*H.mapsCell{kk};
    else
        out{kk} = CostL2(H.mapsCell{kk}.sizeout, y(repel{:},kk),w)*H.mapsCell{kk};
    end
    
end

out = CostSummation(out,1);%could weight differently each one of them
end