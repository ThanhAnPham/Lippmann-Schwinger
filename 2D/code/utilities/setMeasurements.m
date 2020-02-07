function [out,M,varargout] = setMeasurements(curr_y,par,simP,noreflection,algo,M,SP,varargin)
% SETMEASUREMENTS Set the measurements for the algorithm
% Optional input : uin to subtract from the measured total field
global roix_brut
%plane measurement or circular (or other)
isplanar = ~hasandis(simP,'fresnel') || hasandis(simP,'planar');
isintensity = strcmpi(algo(end),'i') || strcmpi(algo,'lfr');
FFTGtild = ~(strcmpi(algo,'lipp') || strcmpi(algo,'lippi') ||...
    strcmpi(algo,'seagle') || strcmpi(algo,'seaglei') ||...
    strcmpi(algo,'rytj') || strcmpi(algo,'rytji') ||...
    strcmpi(algo,'bpmg') || strcmpi(algo,'bpmgi'))...
    || hasandis(par,'discreteGtild');
hasUin = nargin > 7 && ~isempty(varargin{1});
doSubtract = hasUin && par.yscat && ~isintensity;%for complex Lippmann, data measurements are the scattered fields (:= Utotal - uin)

if hasUin
    uinin = varargin{1};
end

if ~isplanar
    out = curr_y;
    if doSubtract
        if nargout > 2
            [out,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
        else
            [out] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
        end
        fprintf('Substracting uin to data measurements : scattered field\n');
    elseif hasUin && nargout > 2
        [~,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
    end
else
    if noreflection
        out = curr_y(SP.z>=par.Lzext/2,:);%no reflection measurement
    else
        out = curr_y;
    end
    %Reduce from brut data
    roix_brut = floor(simP.Lxroi/simP.dx/2);
    if par.SP.nCamera==1 || noreflection
        out = out(max(1 + ceil(end/2) - roix_brut,1):ceil(end/2) + roix_brut,:);
    elseif par.SP.nCamera==2
        out = out([max(ceil(end/4) - roix_brut,1):ceil(end/4) + roix_brut,...
            max(1 + end/2 + ceil(end/4) - roix_brut,end/2):end/2 + ceil(end/4) + roix_brut],:);
    end
    
    %roix = round(par.Lxext/simP.dx/2);
    %out = out(1 + end/2 - roix:end/2 + roix,:);
    out = out(:,par.curr_thetas);
    if FFTGtild %lipp(i),rytj(i) use Gtild
        %smaller ROI than available data measurement FOV, (avoids border artifacts)
        %Sometimes, round(par.Lxext/simP.dx/2) > size(out,1)/2 for BPM because of
        %a trick used to reduce/avoid periodic boundary effects due to uses of FFT
        if ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
            
            roix = min(floor(par.Lxext/simP.dx/2),size(out,1)/(par.SP.nCamera - noreflection)/2);
            if par.SP.nCamera==1 || noreflection
                out = out(max(1 + ceil(end/2) - roix,1):ceil(end/2) + roix - mod(floor(par.Lxext/simP.dx),2),:);
            elseif par.SP.nCamera==2
                out = out([max(ceil(end/4) - roix,1):ceil(end/4) + roix - mod(floor(par.Lxext/simP.dx),2),...
                    max(1 + end/2 + ceil(end/4) - roix,end/2):end/2 + ceil(end/4) + roix - mod(floor(par.Lxext/simP.dx),2)],:);
            end
            
            if doSubtract
                if nargout > 2
                    [out,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
                else
                    [out] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
                end
                fprintf('Substracting uin to data measurements : scattered field\n');
            elseif hasUin && nargout > 2
                [~,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
            end
            [~,out] = setFields([],out,par,simP,algo,'uM');
        end
        %M.sizein(1)*(2 - (par.noreflection & par.SP.nCamera==2)),'nearest');
        curr_SP = par.SP;
        %curr_SPorig = par.SPorig;
    else
        %Gtild's independent of the algorithm discretisation
        %convert the sensor position into the measurement discretisation
        curr_SP.x = SP.x;%par.SP.x*par.dx;
        curr_SP.z = SP.z;%par.SP.z*par.dz;
        curr_SP.nCamera = par.SP.nCamera;
        curr_SP.Nmeas = par.SP.Nmeas;
        M = setM(simP.Nx,simP.Nz,simP.dx,simP.dz,curr_SP);
        %curr_SP.x = curr_SP.x/simP.dx;
        %curr_SP.z = curr_SP.z/simP.dz;
        %curr_SPorig.x = par.SPorig.x*par.dx/simP.dx;
        %curr_SPorig.z = par.SPorig.z*par.dz/simP.dz;
        
        if doSubtract
            if nargout > 2
                [out,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
            else
                out = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
            end
            fprintf('Substracting uin to data measurements : scattered field\n');
        elseif hasUin && nargout > 2
            [~,varargout{1}] = subtractUin(uinin,out,par,simP,algo,noreflection && isplanar,SP);
        end
        
        [~,out] = setFields([],out,par,simP,algo,'uM');
    end
    %if only a fraction of the measurements are used
    if par.Mroi < 1 && ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
        ind_rm =  abs(curr_SP.x - simP.Lxroi/2) > par.Mroi*simP.Lxroi/2;
        out(ind_rm,:) = [];
        if nargout > 2
            varargout{1}(ind_rm,:) = [];
        end
    end
    if FFTGtild && ~(nargin > 8 && strcmpi(varargin{2},'noproc'))
        out = imresizecol(out,...
            par.SP.Nmeas/(1 + (noreflection && par.SP.nCamera > 1)),'box');
        if nargout > 2
            varargout{1} = imresizecol(varargout{1},...
                par.SP.Nmeas/(1 + (noreflection && par.SP.nCamera > 1)),'box');
        end
    end
end



out = out/par.field_scale;
if size(out,2)~=par.Ntheta
    out = out(:,par.curr_thetas);
end

if isintensity %intensity
    out = abs(out).^2;
end

if noreflection && isplanar && FFTGtild
    M = LinOpSelector(M.sel & repmat(1:par.Nzext,par.Nxext,1) - par.Nzext/2 > 0);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,Muin] = subtractUin(Muin,out,par,simP,algo,planarnoreflec,SP)
if size(Muin,2) < length(par.curr_thetas)
    curr_ind = true(size(Muin,2),1);
else
    curr_ind = par.curr_thetas;
end
if ndims(Muin)~=ndims(out) || any(size(Muin(:,curr_ind))./[1 + (par.SP.nCamera>1 && planarnoreflec),1]~=size(out))
    if exist('roix_brut','var')
        curr_siz = size(Muin(max(ceil(end/2) - roix_brut,1):ceil(end/2) + roix_brut,curr_ind))./[1 + (par.SP.nCamera>1 && planarnoreflec),1];
    else
        curr_siz = size(Muin(:,curr_ind))./[1 + (par.SP.nCamera>1 && planarnoreflec),1];
    end
    if any(curr_siz==size(out))
        if par.SP.nCamera>1 && planarnoreflec
            Muin = Muin(SP.z==simP.Lz,curr_ind);
            Muin = Muin(max(ceil(end/2) - roix_brut,1):ceil(end/2) + roix_brut,:);
        elseif exist('roix_brut','var')
            Muin = Muin(max(ceil(end/2) - roix_brut,1):ceil(end/2) + roix_brut,curr_ind);
        else
            Muin = Muin(:,curr_ind);
        end
    else
        try
            Muin = reshape(LinOpSelector(repmat(M.sel,[1,1,simP.Ntheta]))*Muin,size(out));
            if par.SP.nCamera > 1 && planarnoreflec
                Muin = Muin(SP.z==simP.Lz,curr_ind);%no reflection measurement
            end
        catch
            error('argument 6 wrong size');
        end
    end
else
    if par.SP.nCamera > 1 && planarnoreflec
        Muin = Muin(SP.z>=simP.Lz,curr_ind);%no reflection measurement
    else
        Muin = Muin(:,curr_ind);
    end
end

out = out - Muin;
par.distUin = par.distuM;
Muin = setFields(Muin,[],par,simP,algo,'uin');
Muin = squeeze(Muin);
end