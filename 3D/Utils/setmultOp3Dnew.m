function out = setmultOp3Dnew(uin,M,par,algo,uem,varargin)
% SETMULTOP3D : set the stack of (Linear) Operators depending on the forward
% model

Ntheta = par.Ntheta;%size(uin,ndims(uin));
out = cell(Ntheta,1);
switch algo
    case {'lipp','lippi','rytj','rytji','born','borni',...
            'seagle','seaglei','bpmgi','bpmg','helm','lippSMLM'}
        G = FreeSpaceKernel3D(par.k*par.dz,[par.Nx,par.Ny,par.Nz],par.ForwardModel);
        iscentered = false; Zcorr = 0;
        if par.discreteGtild
            Gtil = Gtild3DFFT(par.k*par.dz, par.siz, par.siz_ext,...
                2*par.siz_ext, M,par.ForwardModel,false);
            ForwPos = par.Lzext/2;
        else
            if nargin > 5 && ~isempty(varargin{1})
                if nargin > 6 && ~isempty(varargin{2})
                    ddx = varargin{2}.dx/2;
                    ddy = varargin{2}.dy/2;
                    ddz = varargin{2}.dz/2;
                else
                    ddx = 0;
                    ddy = 0;
                    ddz = 0;
                end
                curr_SP.x = varargin{1}.x;
                curr_SP.y = varargin{1}.y;
                curr_SP.z = varargin{1}.z;%dx,dz should be identical
                curr_SP.Nmeas = varargin{1}.Nmeas;
            else %imprecise because of rounderror
                curr_SP.x = par.SP.x*par.dx;% - par.dx/2;
                curr_SP.y = par.SP.y*par.dy;% - par.dy/2;
                curr_SP.z = par.SP.z*par.dz;%dx,dy,dz should be identical
                curr_SP.Nmeas = par.SP.Nmeas;
            end
            if prod(curr_SP.Nmeas)~=length(curr_SP.x)
                fprintf('Different sensor positions for each view\n');
                %Ntheta = length(curr_SP.x)/curr_SP.Nmeas;
                Gtil = cell(Ntheta,1); sensors_sets = Gtil;
                ind_shared_sensors = zeros(Ntheta,1);
                sensors_ind = 1;
                new_sensors = true;
                for kk = 1:Ntheta
                    curr_ind = find(par.curr_thetas,kk);curr_ind = curr_ind(end);
                    subSP.x = curr_SP.x(1 + (curr_ind-1)*curr_SP.Nmeas:curr_ind*curr_SP.Nmeas);
                    subSP.y = curr_SP.y(1 + (curr_ind-1)*curr_SP.Nmeas:curr_ind*curr_SP.Nmeas);
                    subSP.z = curr_SP.z(1 + (curr_ind-1)*curr_SP.Nmeas:curr_ind*curr_SP.Nmeas);
                    sensors_sets{kk} = [subSP.x,subSP.y,subSP.z];
                    for ll = 1:kk-1
                        if ~any(sensors_sets{kk}(:) - sensors_sets{ll}(:))
                            ind_shared_sensors(kk) = ind_shared_sensors(ll);
                            new_sensors = false;
                            break;
                        end
                    end
                    if new_sensors
                        fprintf('Setting a new sensor set...%i\n',sensors_ind)
                        
                        Gtil{sensors_ind} = FreeSpaceKernel3DtildeMat(par.k*par.dz,...
                            par.siz,subSP,[par.dx,par.dy,par.dz],...[varargin{2}.dx,varargin{2}.dy,varargin{2}.dz],...
                            true,[varargin{2}.dx,varargin{2}.dy,varargin{2}.dz],1,0);
                        ind_shared_sensors(kk) = sensors_ind;
                        sensors_ind = sensors_ind + 1;
                    end
                    
                    curr_pos = [varargin{2}.EmitterPos(curr_ind,1),varargin{2}.EmitterPos(curr_ind,2),...
                        varargin{2}.EmitterPos(curr_ind,3)];
                    curr_k = zeros(1,3);
                    curr_k(2:3) = varargin{2}.k*[sin(par.thetas(kk,1));sin(par.thetas(kk,2))];
                    curr_k(1) = (-1)^(varargin{2}.EmitterPos(curr_ind,3) > 0)*sqrt(varargin{2}.k^2 - sum(curr_k(2:3).^2));
                    curr_k = varargin{2}.dist*curr_k([2:3,1])/varargin{2}.k;
                    figure(3);
                    hold on;
                    scatter3(curr_pos(2),curr_pos(1),curr_pos(3),'g');
                    line(curr_pos(2) + [0,curr_k(1)],curr_pos(1) + [0,curr_k(2)],...
                        curr_pos(3) + [0,curr_k(3)]);
                    set(gca,'YDir','Reverse','XDir','Reverse');
                    %axis auto
                    pause(0.001);
                    new_sensors = true;
                end
            else %must be planar
                %Warning: assumes some sensors are at the end (in x, y & z)
                if min(curr_SP.x)>0
                    curr_SP = centerPos3D(curr_SP,max(curr_SP.x),max(curr_SP.y),0);%max(curr_SP.z));
                end
                iscentered = true;
                if nargin > 5 && hasandis(varargin{2},'realData') || strcmp(varargin{2}.ObjectType,'RBC') || strcmp(varargin{2}.ObjectType,'3Dcell')
                    curr_SP.z = par.distuM*ones(prod(curr_SP.Nmeas),1);
                elseif hasandis(par,'Zsen')
                    curr_SP.z = par.Zsen*ones(prod(curr_SP.Nmeas),1);
                end
                if par.noreflection
                    curr_SP.x = curr_SP.x(curr_SP.z >= 0);
                    curr_SP.y = curr_SP.y(curr_SP.z >= 0);
                    curr_SP.z = curr_SP.z(curr_SP.z >= 0);
                    curr_SP.Nmeas = length(curr_SP.x);
                end
                if sign(max(curr_SP.z)) > 0
                    if max(curr_SP.z) < par.Lz/2 %cannot take measurements inside the volume of reconstruction
                        Zcorr = max(curr_SP.z) - (par.Lz/2 + par.dz);
                        curr_SP.z(:) = par.Lz/2 + par.dz;
                    end
                elseif abs(max(curr_SP.z)) < par.Lz/2 %cannot take measurements inside the volume of reconstruction
                    Zcorr = max(curr_SP.z) - (par.Lz/2 + par.dz);
                    curr_SP.z(:) = (par.Lz/2 + par.dz);
                end
                if true %strcmpi(varargin{2}.ObjectType,'RBC') || (varargin{2}.realData && ~hasandis(varargin{2},'fresnel'))
                    if par.Lxext/2 > max(curr_SP.x(:)) %assume square plane 2D measurements
                        curr_SP.x = -(par.Lxext/2-par.dx):par.dx:par.Lxext/2;
                        curr_SP.y = -(par.Lyext/2-par.dy):par.dy:par.Lyext/2;
                        [curr_SP.x,curr_SP.y] = meshgrid(curr_SP.x,curr_SP.y);
                        curr_SP.x = curr_SP.x(:); curr_SP.y = curr_SP.y(:);
                        curr_SP.z = unique(curr_SP.z(:))*ones(size(curr_SP.x));
                        curr_SP.Nmeas = par.siz_ext(1:2);
                    end
                    Gtil = Vainikkotilde(par.k,par.siz,curr_SP,[par.dx,par.dy,par.dz],ones(1,3));%[varargin{2}.dx,varargin{2}.dy,varargin{2}.dz],...
                    %[par.dx/varargin{2}.dx,par.dy/varargin{2}.dy,par.dz/varargin{2}.dz]);
                else
                    Gtil = FreeSpaceKernel3Dtilde(par.k,par.siz,curr_SP,[varargin{2}.dx,varargin{2}.dy,varargin{2}.dz],...
                        [par.dx/varargin{2}.dx,par.dy/varargin{2}.dy,par.dz/varargin{2}.dz]);
                end
                %[par.dx,par.dy,par.dz]);
            end
            ForwPos = max(unique(curr_SP.z)) - ~iscentered*varargin{2}.Lz/2;
        end
        if exist('curr_SP','var')
            MeasPos = max(unique(curr_SP.z)) + Zcorr - ~iscentered*varargin{2}.Lz/2;%measurement positions wrt center
            distProp = MeasPos - ForwPos;
        else
            distProp =  max(unique(varargin{1}.z)) - varargin{2}.Lz/2 + par.distuM...
                - (max(unique(M.idxmax(3)*par.dz)) - par.Lzext/2);
        end
        %distProp = 0;%TMP
        if isfield(par,'distPropShift')
            distProp = distProp + par.distPropShift;
        end
        if distProp~=0
            if par.dx~=par.dy
                error('dx must be equal to dy')
            end
            if true
                prop = LinOpTiltedPropagator(varargin{2}.lambda0, varargin{2}.n0, distProp, par.dx, Gtil.sizeout, zeros(1,2), 'AS');
                if hasandis(par,'pupil')
                    prop = prop*LinOpTiltedPropagator(varargin{2}.lambda0, varargin{2}.n0, 0, par.dx, Gtil.sizeout, zeros(1,2), 'Pupil',par.NA);
                end
            else
                prop = cell(Ntheta,1);
                for kk = 1:Ntheta
                    prop{kk} = LinOpTiltedPropagator(varargin{2}.lambda0, varargin{2}.n0, -distProp, par.dx, [size(uem,1),size(uem,2)], deg2rad(par.thetas(kk,:)), 'AS');
                end
            end
        end
        
end




for kk = 1:Ntheta
    switch algo
        case {'rytj','rytji'}
            if iscell(Gtil)
                out{kk} = OpRyt3D(G,uin(:,:,:,kk),M,Gtil{kk},uem(:,:,kk),...
                    strcmpi(algo,'rytji') | ~par.yscat);
            else
                out{kk} = OpRyt3D(G,uin(:,:,:,kk),M,Gtil,uem(:,:,kk),...
                    strcmpi(algo,'rytji') | ~par.yscat);
            end
            if distProp~=0
                out{kk}.totalfield = true;%total field to propagate, not scat
                out{kk} = prop{kk}*out{kk};
            end
        case {'bpm','bpmi'}
            out{kk} = OpBPM3D(par.siz,par.Lx,par.Ly,par.Lz,par.k0,par.k,...
                uin(:,:,kk),par.SenDist,par.padfact,M,par.thetas(kk,:));
            if isfield(par,'mod_dist')
                out{kk}.mod_dist = par.mod_dist;
            else
                out{kk}.mod_dist = 0;
            end
        case {'lippSMLM'}
            if iscell(uin)
                iuin = uin{kk};
            else
                iuin = uin(:,:,:,kk);
            end
            if iscell(Gtil)
                out{kk} = OpLipp3DSMLM(G,Gtil{ind_shared_sensors(kk)},M,...
                    iuin,par.Niterin,par.xtolin,uem(:,kk),~par.yscat);
            else
                out{kk} = OpLipp3DSMLM(G,Gtil,M,iuin,par.Niterin,par.xtolin,...
                    uem(:,:,kk),~par.yscat);
            end
            
            out{kk}.loadUin = par.loadUinOnline;
            out{kk}.weightGrad = par.weightGrad;
            if hasandis(par,'useLocalGPU')
                out{kk}.setGPU(par.useLocalGPU);
            end
            if distProp~=0
                if iscell(prop)
                    %out{kk}.totalfield = false;%scat field to propagate
                    out{kk} = prop{kk}*out{kk};
                else
                    out{kk} = prop*out{kk};
                end
            end
            if ~iscell(Gtil) && any(M.sizeout - out{kk}.sizeout)
                out{kk} = LinOpSelectorPatch(out{kk}.sizeout, 1 + (out{kk}.sizeout - M.sizeout)/2,(out{kk}.sizeout + M.sizeout)/2)*out{kk};
            end
        case {'lipp','lippi'}
            if iscell(Gtil)
                out{kk} = OpLipp3D(G,Gtil{ind_shared_sensors(kk)},M,uin(:,:,:,kk),par.Niterin,par.xtolin,uem(:,kk),...
                    strcmpi(algo,'lippi') | ~par.yscat);
            else
                out{kk} = OpLipp3D(G,Gtil,M,uin(:,:,:,kk),par.Niterin,par.xtolin,uem(:,:,kk),...
                    strcmpi(algo,'lippi') | ~par.yscat);
            end
            if hasandis(par,'useLocalGPU')
                out{kk}.setGPU(par.useLocalGPU);
            end
            if distProp~=0
                if iscell(prop)
                    %out{kk}.totalfield = false;%scat field to propagate
                    out{kk} = prop{kk}*out{kk};
                else
                    out{kk} = prop*out{kk};
                end
            end
            if ~iscell(Gtil) && any(M.sizeout - out{kk}.sizeout)
                out{kk} = LinOpSelectorPatch(out{kk}.sizeout, 1 + (out{kk}.sizeout - M.sizeout)/2,(out{kk}.sizeout + M.sizeout)/2)*out{kk};
            end
        case 'helm'
            if iscell(Gtil)
                out{kk} = OpHelm3D(G,Gtil{ind_shared_sensors(kk)},M,uin(:,:,:,kk),par.Niterin,par.xtolin,uem(:,kk),...
                    strcmpi(algo,'lippi') | ~par.yscat);
            else
                out{kk} = OpHelm3D(G,Gtil,M,uin(:,:,:,kk),par.Niterin,par.xtolin,uem(:,:,kk),...
                    strcmpi(algo,'lippi') | ~par.yscat);
            end
            if hasandis(par,'useLocalGPU')
                out{kk}.setGPU(par.useLocalGPU);
            end
            if distProp~=0
                if iscell(prop)
                    %out{kk}.totalfield = false;%scat field to propagate
                    out{kk} = prop{kk}*out{kk};
                else
                    out{kk} = prop*out{kk};
                end
            end
        otherwise
            warning('Non existent model');
    end
    if hasandis(par,'memoizeApply')
        out{kk}.memoizeOpts.apply = par.memoizeApply;
    end
end
out = StackMap(out);
end