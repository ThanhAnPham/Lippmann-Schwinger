function [uin,varargout] = setUin3D(par,simP,algo,uinin,varargin)
%% SETUIN3D Set the 3D incident fields (simulated or numerically propagated)
if nargout > 1
    varargout{1} = [];
end

switch par.uin
    case 'sim'
        doGB = hasandis(simP,'GaussianBeam') || (nargin > 4 && strcmpi(varargin{1},'GB'));
        if strcmp(algo,'bpm') || strcmp(algo,'bpmi') || strcmp(algo,'lfr')
            Nx = par.Nxext; Ny = par.Nyext;
            if doGB
                Nxext = simP.padfact(1)*simP.Nx; Nyext = simP.padfact(2)*simP.Ny;
            else
                Nxext = Nx; Nyext = Ny;
            end
            %Nx = par.Nxext; Ny = par.Nyext; Nz = par.Nzext;
            
            uin = zeros(Nx,Ny,par.Ntheta);
        else
            Nx = par.Nx; Ny = par.Ny; Nz = par.Nz;
            Nxext = par.Nxext; Nyext = par.Nyext; %Nzext = par.Nzext;
            uin = zeros(Nx,Ny,Nz,par.Ntheta);
        end
        
        curr_k = zeros(1,3);
        if doGB
            dx = simP.dx; dy = simP.dy;
            k = par.k;
        else
            k = par.k;
            dx = par.dx; dy = par.dy; dz = par.dz;
        end
        for kk = 1:par.Ntheta
            angx = par.thetas(kk,1);
            angy = par.thetas(kk,2);
            curr_k(2:3) = k*[sin(angx);sin(angy)];%[x,y], requires dx = dy
            if isfield(simP,'EmitterPos')
                curr_k(1) = (-1)^(simP.EmitterPos(kk,3) > 0)*sqrt(k^2 - sum(curr_k(2:3).^2));
            else
                curr_k(1) = (-1)^((abs(angx) > pi/2) + (abs(angy) > pi/2))*sqrt(k^2 - sum(curr_k(2:3).^2));
            end
            if doGB
                Bx = simP.sigmaBeam*simP.dx*Nxext/simP.padfact(1); % beam scale (m)
                Bx = Bx/cos(angx); % projection of the propagation plane to XY plane, sin(pi/2 - angx)
                By = simP.sigmaBeam*simP.dy*Nyext/simP.padfact(2); % beam scale (m)
                By = By/cos(angy); % projection of the propagation plane to XY plane, sin(pi/2 - angx)
                %Bz = simP.sigmaBeam*simP.Lz;
                %Bz = -Bz/(sin(angx)*sin(angy));
            end
            if hasandis(simP,'uin_photon')
                ain = simP.uin_photon(min(kk,end));
            else
                ain = 1;
            end
            
            Vectx = ((0:Nxext-1) - Nxext/2)*dx;
            Vecty = ((0:Nyext-1) - Nyext/2)*dy;
            [Y,X] = meshgrid(Vecty,Vectx);
            if doGB
                ain = ain*exp(-((X./Bx).^2 + (Y./By).^2));
                curr_uin = ain.*exp(1i*(curr_k(2)*X + curr_k(3)*Y));
                if strcmp(algo,'bpm') || strcmp(algo,'bpmi') || strcmp(algo,'lfr')
                    curr_par.distUin = -par.Lz/2;
                else
                    curr_par = par;
                    %curr_par.distUin = -par.Lzext/2:par.dz:par.Lzext/2 - par.dz;
                    curr_par.Nz = par.Nzext;
                end
                curr_uin = setFields3D(curr_uin,[],curr_par,simP,algo,'uin');
            elseif strcmp(algo,'bpm') || strcmp(algo,'bpmi') || strcmp(algo,'lfr')
                curr_uin = ain.*exp(1i*(curr_k(2)*X + curr_k(3)*Y - curr_k(1)*par.Lz/2));
            else
                curr_uin = PlaneWave3D(curr_k([2:3,1]).*[dx,dy,dz]);
                curr_uin = ain.*curr_uin.Apply(Nx,Ny,Nz);
            end
                
            if strcmp(algo,'bpm') || strcmp(algo,'bpmi') || strcmp(algo,'lfr')
                if doGB
                    uin(:,:,kk) = setInput2D(curr_uin,simP,par,[Nx,Ny],0);
                else
                    uin(:,:,kk) = curr_uin;
                end
                if hasandis(par,'apod') && par.apod > 0
                    uin(:,:,kk) = (tukeywin(Nx,par.apod)*tukeywin(Nx,par.apod)').*uin(:,:,kk);
                end
            else
                uin(:,:,:,kk) = curr_uin;
            end
        end
        if hasandis(simP,'fresnel')
            for kk = 1:par.Ntheta
                indoi = 5 + 9*floor((par.curr_ind_thetas(kk)-1)/9);
                angx = simP.thetas(indoi,1);
                angy = simP.thetas(indoi,2);
                SP = varargin{1};
                x = SP.x(1 + (indoi -1)*SP.Nmeas:indoi*SP.Nmeas);%probably inverted
                y = SP.y(1 + (indoi -1)*SP.Nmeas:indoi*SP.Nmeas);
                z = SP.z(1 + (indoi -1)*SP.Nmeas:indoi*SP.Nmeas);
                cuin = uinin(:,indoi);
                curr_k(2:3) = k*[sin(angx);sin(angy)];%[x,y], requires dx = dy
                if isfield(simP,'EmitterPos')
                    curr_k(1) = (-1)^(simP.EmitterPos(indoi,3) > 0)*sqrt(k^2 - sum(curr_k(2:3).^2));
                else
                    curr_k(1) = (-1)^((abs(angx) > pi/2) + (abs(angy) > pi/2))*sqrt(k^2 - sum(curr_k(2:3).^2));
                end
                theta0 = angle(cuin(14)./exp(1i*(curr_k(1)*z(14) + curr_k(2)*y(14) + curr_k(3)*x(14))));
                a = max(abs(cuin(:)));
                %theta0 = 0;%-3.129078009066249e-01;
                uin(:,:,:,kk) = a*exp(1i*theta0)*uin(:,:,:,kk);
            end
        end
    case {'real','Soulez'}
        if strcmp(par.uin,'Soulez')
            %%
            uin = zeros([par.siz,par.Ntheta]);
            
            sz_uinin = [size(uinin,1),size(uinin,2)];
            %%
            if hasandis(par,'doUinApod')
                t = 0.3;
                Napod = 154;
                bounwin = tukeywin(Napod + par.siz(1),t).*tukeywin(Napod + par.siz(1),t)';
                bounwin = padarray(bounwin,(sz_uinin - size(bounwin))/2,0,'both');
            else
                bounwin = 1;
            end
            %figure,imagesc(bounwin);axis image
            %%
            if hasandis(par,'uinpad')
                pfact = par.uinpad;
            else
                pfact = 1;
            end
            if hasandis(par,'doUinApod')
                fname = sprintf('%s_uin_lipp_%i_%i_%1.2f_Ntheta_%i_sz_%i%i%i_szuinin_%i_%i_dx_%1.2e.mat',...
                    simP.ObjectType,pfact,Napod,t,par.Ntheta,par.siz,sz_uinin,simP.dx);
            else
                fname = sprintf('%s_uin_lipp_%i_Ntheta_%i_sz_%i%i%i_szuinin_%i_%i_dx_%1.2e.mat',...
                    simP.ObjectType,pfact,par.Ntheta,par.siz,sz_uinin,simP.dx);
            end
            if exist(fname,'file')
                load(fname);
            else
                if ~isfield(simP,'roi')
                    simP.roi = sz_uinin/2;
                end
                Padder = LinOpSelectorPatch(pfact*sz_uinin,1 + (pfact-1)*sz_uinin/2,(pfact+1)*sz_uinin/2)';
                PadderExt = LinOpSelectorPatch([pfact*sz_uinin,length(par.distUin)],[1 + (pfact-1)*sz_uinin/2,1],[(pfact+1)*sz_uinin/2,length(par.distUin)])';
                if par.dx~=simP.dx
                    fsz = simP.dx/par.dx;
                    simP.roi = simP.roi*fsz;
                    C = @(x) imresize3(x,round([fsz*sz_uinin,length(par.distUin)]),'nearest');
                    S = LinOpSelectorPatch(round([sz_uinin*fsz,length(par.distUin)]),...
                        round([1 + simP.roi - par.siz(1:2)/2,1]), round([simP.roi + par.siz(1:2)/2,length(par.distUin)]));
                else
                    C = @(x) x;
                    S = LinOpSelectorPatch([sz_uinin,length(par.distUin)],...
                        [1 + simP.roi - par.siz(1:2)/2,1], [simP.roi + par.siz(1:2)/2,length(par.distUin)]);
                
                end
                
                xx = (simP.dx:simP.dx:(sz_uinin(1)*simP.dx)) - (sz_uinin(1)*simP.dx)/2 - simP.dx;
                yy = (simP.dy:simP.dy:(sz_uinin(2)*simP.dy)) - (sz_uinin(2)*simP.dy)/2 - simP.dx;
                [YYB,XXB] = meshgrid(xx,yy);
                global isGPU
                tmp_isGPU = isGPU;
                useGPU(0);
                for ll = 1:par.Ntheta
                    tic
                    tiltFactorB = exp(1i*simP.k*(sin(par.thetas(ll,1))*XXB + sin(par.thetas(ll,2))*YYB));
                    curr_uinf = fft2(Padder*(bounwin.*uinin(:,:,par.curr_ind_thetas(ll))./tiltFactorB));
                    cthetas = rad2deg(par.thetas(ll,:))';%asin(simP.k*sin(par.thetas(ll,[2,1]))*simP.lambda0/(2*pi).*[simP.Lx,simP.Ly])';
                    P = LinOpTiltedPropagator(simP.lambda0,...
                        simP.n0,0,simP.dx,pfact*sz_uinin,-cthetas,'AS');
                    lf = @(z) abs(z*sqrt(P.mu)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(1)*P.dxy...
                        & abs(z*sqrt(P.mv)./(sqrt((  P.n0/ P.lambda)^2- (P.mu + P.mv)))) <= P.sizein(2)*P.dxy;
                    ephi = arrayfun(@(x) P.mod.*exp(P.phase*x).*lf(x),-par.distUin,'UniformOutput',false);
                    ephi = cat(3,ephi{:});
                    %for kk = 1:length(par.distUin)
                    uin(:,:,:,ll) = S*C((PadderExt'*ifft2(ephi.*curr_uinf)).*tiltFactorB);%*exp(1i*par.k*par.distUin(kk)*2);
                    %end
                    toc
                end
                
                save(fname,'uin','-v7.3');
                
                if tmp_isGPU
                    useGPU(1);
                end
            end
            doCrop = false;
        else
            if numel(uinin)*length(par.distUin)*16 > 4e9
                
                
                if isfield(simP,'roi') && exist(fname,'file')
                    fname = sprintf('uin_lipp_Ntheta_%i_roi_%i_%i_real.mat',par.Ntheta,simP.roi);
                    load(fname);
                else
                    uin = zeros([par.siz,par.Ntheta]);
                    for kk = 1:par.Ntheta
                        currui = setFields3D(uinin(:,:,par.curr_ind_thetas(kk)),[],par,simP,algo,'uin');
                        roix = floor(par.Nx*par.dx/simP.dx/2);
                        roiy = floor(par.Ny*par.dy/simP.dy/2);
                        if isfield(simP,'roi')
                            roi_center = simP.roi;
                        else
                            roi_center = size(currui); roi_center = ceil(roi_center(1:2)/2);
                        end
                        uin(:,:,:,kk) = currui(1 + roi_center(1) - roix:roi_center(1) + roix,...
                            1 + roi_center(2) - roiy:roi_center(2) + roiy,:);
                    end
                    fname = sprintf('uin_lipp_Ntheta_%i_roi_%i_%i_real.mat',par.Ntheta, roi_center);
                    save(fname,'uin','-v7.3');
                end
                doCrop = false;
            else
                uin = setFields3D(uinin(:,:,par.curr_thetas),[],par,simP,algo,'uin');
                doCrop = true;
            end
        end
        if (strcmp(algo,'bpm') || strcmp(algo,'bpmi')) || strcmp(algo,'lfr')
            roix = floor(par.Lxext/simP.dx/2);%floor(par.Lxext/(simP.dx*size(uin,1))*simP.Nx/2);
            roiy = floor(par.Lyext/simP.dy/2);%floor(par.Lyext/(simP.dy*size(uin,2))*simP.Ny/2);
            uin = uin(1 + ceil(end/2) - roix:ceil(end/2) + roix,...
                1 + ceil(end/2) - roiy:ceil(end/2) + roiy, :);
            uin = imresize3slice(uin,[par.Nxext,par.Nyext],'nearest');
        elseif doCrop
            if isfield(simP,'roi')
                roi_center = simP.roi;
            else
                roi_center = size(uin); roi_center = ceil(roi_center(1:2)/2);
            end
            roix = floor(par.Nx*par.dx/simP.dx/2);
            roiy = floor(par.Ny*par.dy/simP.dy/2);
            uin = uin(1 + roi_center(1) - roix:roi_center(1) + roix,...
                1 + roi_center(2) - roiy:roi_center(2) + roiy,:, :);
            expsz = [par.siz,par.Ntheta];
            expsz(expsz==1) = [];
            if any(size(uin) - expsz)
                fprintf('Should not happen');
                uin = imresize3volume(uin,par.siz,'box');
            end
        end
        
        uin = uin/par.field_scale;
end
if contains(algo,'bpmg')
    uin = squeeze(uin(:,round(1 + end/2 - 0.5*par.Lx/simP.Lx*length(uin)),:));
end
end