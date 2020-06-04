function [out,M,varargout] = setMeasurements3D(curr_y,par,simP,noreflection,algo,M,SP,varargin)
%% SETMEASUREMENTS : Set the measurements for the algorithm

global roix_brut roiy_brut
%plane measurement or circular (or other)
isplanar = ~hasandis(simP,'fresnel') || hasandis(simP,'planar');
isintensity = strcmpi(algo(end),'i') || strcmpi(algo,'lfr');
FFTGtild = ~(strcmpi(algo,'bpm') || strcmpi(algo,'lipp') || strcmpi(algo,'lippi') ||...
    strcmpi(algo,'seagle') || strcmpi(algo,'seaglei') ||...
    strcmpi(algo,'rytj') || strcmpi(algo,'rytji') ||...
    strcmpi(algo,'bpmg') || strcmpi(algo,'bpmgi'))...
    || hasandis(par,'discreteGtild');
hasUin = nargin > 7 && ~isempty(varargin{1});
doSubtract = hasUin && par.yscat && ~isintensity;%for complex Lippmann, data measurements are the scattered fields (:= Utotal - uin)
methrsz = 'box';%'nearest';

if hasUin
    uinin = varargin{1};
end

if ~isplanar
    out = curr_y(:,par.curr_thetas);
    if par.yscat
        out = out - uinin(:,par.curr_thetas);
    end
else
    if noreflection && (hasandis(SP,'Ncam') && SP.Ncam > 1) %if there is no field of this, 1 camera plane
        out = curr_y(:,:,unique(SP.z) >= simP.Lz/2,:);
        if hasUin
            uinin = uinin(:,:,unique(SP.z) >= simP.Lz/2,:);
        end
    else
        out = curr_y;
    end
    %Reduce from brut data but useless since setFields is called after
    %smaller ROI reduction
    roix_brut = floor(simP.Lxroi/simP.dx/2);
    roiy_brut = floor(simP.Lyroi/simP.dy/2);
    if isfield(simP,'roi')
        roi_center = simP.roi;
    else
        roi_center = size(out); roi_center = ceil(roi_center(1:2)/2);
    end
    if par.SP.nCamera==1
        out = out(max(1 + roi_center(1) - roix_brut,1):roi_center(1) + roix_brut,...
            max(1 + roi_center(2) - roiy_brut,1):roi_center(2) + roiy_brut,:);
    elseif par.SP.nCamera>=2
        out = out(max(1 + roi_center(1) - roix_brut,1):roi_center(1) + roix_brut,...
            max(1 + roi_center(2) - roiy_brut,1):roi_center(2) + roiy_brut,:,:);
    end
    if hasUin
        if par.SP.nCamera==1
            uinin = uinin(max(1 + roi_center(1) - roix_brut,1):min(roi_center(1) + roix_brut,end),...
                max(1 + roi_center(2) - roiy_brut,1):min(roi_center(2) + roiy_brut,end),:);
        elseif par.SP.nCamera==2
            if noreflection
                uinin = uinin(max(1 + roi_center(1) - roix_brut,1):roi_center(1) + roix_brut,...
                    max(1 + roi_center(2) - roiy_brut,1):roi_center(2) + roiy_brut,2,:);
            else
                uinin = uinin(max(1 + roi_center(1) - roix_brut,1):roi_center(1) + roix_brut,...
                    max(1 + roi_center(2) - roiy_brut,1):roi_center(2) + roiy_brut,:,:);
            end
        end
    end
    if par.SP.nCamera > 1
        out = out(:,:,:,par.curr_thetas);
    else
        out = out(:,:,par.curr_thetas);
    end
    if FFTGtild %lipp(i),rytj(i) use Gtild
        %smaller ROI than available data measurement FOV, (avoids border artifacts)
        %Sometimes, round(par.Lxext/simP.dx/2) > size(out,1)/2 for BPM because of
        %a trick used to reduce/avoid periodic boundary effects due to uses of FFT
        if ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
            [~,out] = setFields3D([],out,par,simP,algo,'uM');
            roix = min(floor(par.Lxext/simP.dx/2),size(out,1)/2);
            roiy = min(floor(par.Lyext/simP.dy/2),size(out,2)/2);
            if par.SP.nCamera==1
                out = out(max(1 + ceil(end/2) - roix,1):ceil(end/2) + roix - mod(floor(par.Lxext/simP.dx),2),...
                    max(1 + ceil(end/2) - roiy,1):ceil(end/2) + roiy - mod(floor(par.Lyext/simP.dy),2),:);
            elseif par.SP.nCamera>=2
                out = out(max(1 + ceil(end/2) - roix,1):ceil(end/2) + roix - mod(floor(par.Lxext/simP.dx),2),...
                    max(1 + ceil(end/2) - roiy,1):ceil(end/2) + roiy - mod(floor(par.Lyext/simP.dy),2),:,:);
                
            end
            if doSubtract
                
                [~,uinin] = setFields3D([],uinin,par,simP,algo,'uM');
                if par.SP.nCamera==1
                    uinin = uinin(max(1 + ceil(end/2) - roix,1):ceil(end/2) + roix - mod(floor(par.Lxext/simP.dx),2),...
                        max(1 + ceil(end/2) - roiy,1):ceil(end/2) + roiy - mod(floor(par.Lyext/simP.dy),2),:);
                elseif par.SP.nCamera>=2
                    uinin = uinin(max(1 + ceil(end/2) - roix,1):ceil(end/2) + roix - mod(floor(par.Lxext/simP.dx),2),...
                        max(1 + ceil(end/2) - roiy,1):ceil(end/2) + roiy - mod(floor(par.Lyext/simP.dy),2),:,:);
                end
                
                if nargout > 2
                    [out,varargout{1}] = ProcAndSubtractUin(uinin,out,par,simP,algo);
                else
                    [out] = ProcAndSubtractUin(uinin,out,par,simP,algo);
                end
                fprintf('Substracting uin to data measurements : scattered field\n');
            elseif hasUin && nargout > 2
                [~,varargout{1}] = ProcAndSubtractUin(uinin,out,par,simP,algo,false);
            end
            
        end
        curr_SP = par.SP;
        
    else
        %Gtild's independent of the algorithm discretisation
        %convert the sensor position into the measurement discretisation
        curr_SP.x = SP.x;
        curr_SP.y = SP.y;
        curr_SP.z = SP.z;
        curr_SP.nCamera = par.SP.nCamera;
        curr_SP.Nmeas = par.SP.Nmeas;
        if strcmpi(algo,'lipp')
            M = LinOpIdentity(SP.Nmeas);%setM3D(simP.Nx,simP.Ny,simP.Nz,simP.dx,simP.dy,simP.dz,curr_SP);
        end
        if doSubtract
            if nargout > 2
                [out,varargout{1}] = ProcAndSubtractUin(uinin,out,par,simP,algo);
            else
                out = ProcAndSubtractUin(uinin,out,par,simP,algo);
            end
            fprintf('Substracting uin to data measurements : scattered field\n');
        elseif hasUin && nargout > 2
            [~,varargout{1}] = ProcAndSubtractUin(uinin,out,par,simP,algo,false);
        end
        
        [~,out] = setFields3D([],out,par,simP,algo,'uM');
    end
    %if only a fraction of the measurements are used
    if par.Mroi(1) < 1 && ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
        warning('not adapted for several planes measurements');
        ind_rm =  repmat(reshape(abs(curr_SP.x - simP.Lxroi/2) > par.Mroi(1)*simP.Lxroi/2,[simP.Nx,simP.Ny]),[1,1,simP.Ntheta]);
        
        out(ind_rm) = [];
        if nargout > 2
            varargout{1}(ind_rm) = [];
        end
    end
    if par.Mroi(2) < 1 && ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
        warning('not adapted for several planes measurements');
        ind_rm =  repmat(reshape(abs(curr_SP.y - simP.Lyroi/2) > par.Mroi(2)*simP.Lyroi/2,[simP.Nx,simP.Ny]),[1,1,simP.Ntheta]);
        out(ind_rm) = [];
        if nargout > 2
            varargout{1}(ind_rm) = [];
        end
    end
    if FFTGtild && ~(nargin > 8 && strcmpi(varargin{2},'noproc')) || (~FFTGtild && par.dx~=simP.dx)
        if par.SP.nCamera==1
            out = imresize3slice(out,[par.SP.Nmeas,par.Ntheta],methrsz);
            if nargout > 2
                varargout{1} = imresize3slice(varargout{1},[par.SP.Nmeas,par.Ntheta],methrsz);
            end
        elseif par.SP.nCamera>=2
            if strcmpi(algo,'bpm')
                out = imresize3volume(out,[par.Nx,par.Ny,2 - noreflection,par.Ntheta],methrsz);
            else
                out = imresize3volume(out,[par.SP.Nmeas,2 - noreflection,par.Ntheta],methrsz);
            end
            if nargout > 2
                varargout{1} = imresize3volume(varargout{1},[par.SP.Nmeas,2 - noreflection,par.Ntheta],methrsz);
            end
        end
    end
end

out = squeeze(out);
if nargout > 2
    if isplanar
        varargout{1} = squeeze(varargout{1});
    else
        varargout{1} = varargin{1};
    end
end
out = out/par.field_scale;
%if size(out,ndims(out))~=par.Ntheta
%    repel = repmat({':'},1,ndims(out) - 1);
%    out = out(repel{:},par.curr_thetas);
%end

if isintensity %intensity
    out = abs(out).^2;
end

if isa(M,'LinOpComposition')
    M_curr = M.H2;
else
    M_curr = M;
end
if isplanar && noreflection
    if FFTGtild
        zpos = arrayfun(@(x) x*ones(par.Nxext,par.Nyext),1:par.Nzext,'UniformOutput',false);
        zpos = cat(3,zpos{:});
        if par.dx < simP.dx %cannot be done with a LinOpSelectorPatch
            M = LinOpSelector(M_curr.sel{1} & (zpos - par.Nzext/2) > 0);
            M = LinOpShape(M.sizeout, [size(out,1),size(out,2)])*M;
        else
            if hasandis(simP,'realData')
                ind = find(M_curr.sel{1} & (zpos - par.distuM/par.dz) > 0);
            else
                ind = find(M_curr.sel{1} & (zpos - par.Nzext/2) >= 0);
            end
            [indx,indy,indz] = ind2sub(M.sizein,ind);
            
            SPax = find(any(any(M_curr.sel{1},1),2));
            
            if any(SPax > par.Nzext/2 - par.Nz/2 & SPax < par.Nzext/2 + par.Nz/2)
                fprintf('Measurements are inside the domain, not possible to get it directly. Rather take them at par.Nzext*par.dz\n');
                M = LinOpSelectorPatch(M.sizein,min([indx,indy,par.Nzext*ones(size(indz))]),max([indx,indy,par.Nzext*ones(size(indz))]));
            else
                M = LinOpSelectorPatch(M.sizein,min([indx,indy,indz]),max([indx,indy,indz]));
            end
        end
    elseif par.dx~=simP.dx 
        if strcmpi(algo,'lipp') %probably wrong
            roi_center = M.sizein/2;
            roix_ratio = floor(simP.Lxroi/par.dx/2);
            roiy_ratio = floor(simP.Lyroi/par.dy/2);
            M = LinOpSelectorPatch(M.sizein, max(1 + [roi_center(1) - roix_ratio,roi_center(2) - roiy_ratio],1),...
                [roi_center(1) + roix_ratio, roi_center(2) + roiy_ratio]);
        %elseif strcmpi(algo,'bpm')
        %    M = LinOpSelectorPatch(M.sizein, max(1 + [roi_center(1) - roix_ratio,roi_center(2) - roiy_ratio],1),...
        %        [roi_center(1) + roix_ratio, roi_center(2) + roiy_ratio]);
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out,Muin] = ProcAndSubtractUin(Muin,out,par,simP,algo,varargin)

doSubtract = true;
if nargin > 5
    doSubtract = false;
end

if ndims(Muin) > 2 && size(Muin,ndims(Muin)) < length(par.curr_thetas)
    curr_ind = true(size(Muin,ndims(Muin)),1);
elseif ndims(Muin)==2 && par.Ntheta==1
    curr_ind = 1;
else
    curr_ind = par.curr_thetas;
    
end
comp_sz = ones(1,2 + (par.SP.nCamera > 1) + (nnz(curr_ind)>1));

repCol = repmat({':'},1,2 + par.SP.nCamera - 1);

if ndims(Muin)~=ndims(out) || any(size(Muin)./comp_sz~=size(out))
    error('Not the good size uem vs y');
else
    Muin = Muin(repCol{:},curr_ind);
end
if doSubtract
    out = out - Muin;
end
par.distUin = par.distuM;
Muin = setFields3D(Muin,[],par,simP,algo,'uem');
end