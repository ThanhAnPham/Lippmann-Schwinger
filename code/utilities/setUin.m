function [uin,varargout] = setUin(par,simP,algo,uinin,varargin)

if nargout > 1
    varargout{1} = [];
end

switch par.uin
    case 'sim'
        Nx = par.Nxext; Nz = par.Nzext;
        %uin_rytj = setUin(par.Nxext,par.Nzext,par.kdz,par,simP);
        uin = zeros(Nx,Nz,par.Ntheta);
        
        [ZZ,XX] = meshgrid(par.dz*(-Nz/2+1:Nz/2),...
            par.dx*(-Nx/2+1:Nx/2));
        
        doGB = hasandis(simP,'GaussianBeam') || (nargin > 4 && strcmpi(varargin{1},'GB'));
        for kk = 1:par.Ntheta
            angx = par.thetas(kk);
            curr_k = par.kdz*[cos(angx);sin(angx)];%[z,x]
            if doGB
                Bx = simP.sigmaBeam*simP.Lx; % beam scale (m)
                Bx = Bx/cos(angx); % projection of the propagation plane to XY plane
                Bz = simP.sigmaBeam*simP.Lz; Bz = -Bz/sin(angx);%"-" because from (z',x') -> (z,x)
                ain = exp(-(XX./Bx + ZZ./Bz).^2); % illumination beam amplitude
            elseif (nargin > 4 && strcmpi(varargin{1},'tukey'))
                ain = repmat(tukeywin(Nx,1/par.padfact),1,Nz);
            elseif hasandis(simP,'uin_photon')
                ain = simP.uin_photon;
            else
                ain = 1;
            end
            u = PlaneWave(curr_k);
            uin(:,:,kk) = ain.*u.Apply(Nx,Nz);
        end
        if strcmp(algo,'bpm') || strcmp(algo,'bpmi')
            uin = squeeze(uin(:,1 + round(1000*par.SenDist/par.dz)/1000,:));
        elseif strcmp(algo,'lfr')
            uin = squeeze(uin(:,unique(par.SP.z),:));
        end
        
    case 'real'
        uin = setFields(uinin,[],par,simP,algo,'uin');
        roix = floor(par.Lxext/(simP.dx*size(uin,1))*simP.Nx/2);
        %roix = floor(par.Lxext/simP.Lx*simP.Nx/2);
        if (strcmp(algo,'bpm') || strcmp(algo,'bpmi')) || strcmp(algo,'lfr')
            uin = uin(max(ceil(end/2) - roix,1):ceil(end/2) + roix,par.curr_thetas);
            uin = imresizecol(uin,par.Nxext,'nearest');
        else
            uin = uin(1 + ceil(end/2) - roix:end/2 + roix,:,par.curr_thetas);
            %roiz = floor(par.Lzext/(simP.dz*length(par.distUin))*/2);
            %uin = uin(:,max(ceil(end/2) - roiz,1):ceil(end/2) + roiz,:);
            uin = imresize(uin,par.siz_ext,'box');
        end
        
        uin = uin/par.field_scale;
    case 'cylindrical'
        %2005 Crocco Testing the contrast source extended Born inversion method against real data- the TM case.pdf
        %mu0 = 1.25663706e-6;%m kg s-2 A-2
        %P =  2.5663e+05;
        
        if nargin > 4
            SP = varargin{1};
            fact = 1;
        else
            SP = par.SP;%in discretization unit
            fact = par.dz;
        end
        if SP.Nmeas < size(uinin,1) && isfield(SP,'ind_ang')
            uinin = uinin(SP.ind_ang,:);
        end
        source_pos = simP.EmitterPos/fact;%conversion to discretization unit
        Nsensors = SP.Nmeas;
        
        if isfield(par,'Nc')
            Nc = par.Nc;
        else
            % model ARA DRG 118A from Fresnel (largest dim) 9.5? x 5.6? x 7.9?
            horn_dim = 0.2413;%0.2413 x 0.14224 x 0.20066 meter
            Nc = ceil(par.kdz*horn_dim/par.dz);
            
        end
        fprintf('%i cylindrical wave used to reconstruct Uin \n',2*Nc+1);
        [ZZ,XX] = meshgrid(((0:par.Nzext-1) - par.Nzext/2)*par.dz/fact,...
            ((0:par.Nxext-1) - par.Nxext/2)*par.dx/fact);
        uin = zeros(par.Nxext,par.Nzext,par.Ntheta);
        A = zeros(2*Nc + 1,par.Ntheta);
        
        vec_theta = find(ismember(simP.thetas,par.thetas));
        
        for kk = 1:length(vec_theta)
            curr_iter = vec_theta(kk);
            curr_x = SP.x(1 + (curr_iter-1)*Nsensors:(curr_iter-1)*Nsensors + Nsensors);
            curr_z = SP.z(1 + (curr_iter-1)*Nsensors:(curr_iter-1)*Nsensors + Nsensors);
            curr_source_pos = repmat(source_pos(curr_iter,:),[Nsensors,1]);
            diff_r = [curr_x,curr_z] - curr_source_pos;%Nsensors x 2
            curr_theta = (-1).^(diff_r*[cos(simP.thetas(curr_iter));-sin(simP.thetas(curr_iter))] < 0);
            curr_theta = curr_theta.*(pi - acos(sum(curr_source_pos.*diff_r,2)...
                ./(sqrt(sum(curr_source_pos.^2,2)).*sqrt(sum(diff_r.^2,2)))));
            r_norm = sqrt(sum(diff_r.^2,2));
            if par.doPlot
                figure(8);clf
                scatter(source_pos(curr_iter,2),source_pos(curr_iter,1),'r');hold on;
                scatter(curr_z,curr_x,10,'filled','b');
                scatter(ZZ(:),XX(:),10,'filled','MarkerFaceColor',[0.75,0.25,0]);
                title('System');legend('Source','Sensors','ROI');
                set(gca,'YDir','reverse','XDir','reverse');axis image;hold off;
            end
            Hmat = zeros(Nsensors, 2*Nc + 1);
            for n = -Nc:Nc
                curr_H = besselh(n,1,par.k*fact*r_norm).*exp(1i*n*curr_theta);
                Hmat(:,1 + n + Nc) = curr_H;%build operator
            end
            H = LinOpMatrix(Hmat);
            b = uinin(:,curr_iter);
            fprintf('%i unknowns, %i measurements\n',2*Nc+1,size(Hmat,1))
            b = H'*b;H = H'*H;CG = OptiConjGrad(H,b);
            CG.OutOp = OutputOpti(1,[],50);
            CG.maxiter = 5e3;
            CG.ItUpOut = 0;
            CG.verbose = 0;
            CG.CvOp = TestCvgStepRelative(1e-9);
            CG.run(H'*b);
            
            A(:,kk) = CG.xopt;%Coefficient
            Hx = Hmat*A(:,kk);
            fprintf('view %i/%i. Cost value : %1.2e\n',kk,par.Ntheta,norm(Hx - uinin(:,curr_iter)));
            
            if par.doPlot
                figure(9);clf
                subplot(121);
                plot(abs(Hx),'b-');hold on
                plot(abs(uinin(:,curr_iter)),'g--');
                title('Magnitude');legend('Estimated','Measured');
                subplot(122);
                plot(angle(Hx),'b-');hold on;
                plot(angle(uinin(:,curr_iter)),'g--');
                title('Phase');legend('Estimated','Measured');
            end
            %build Uin
            rx = XX - source_pos(curr_iter,1);
            rz = ZZ - source_pos(curr_iter,2);
            
            curr_theta = (-1).^(rx*cos(simP.thetas(curr_iter)) - rz*sin(simP.thetas(curr_iter)) < 0);
            curr_theta = curr_theta.*(pi - acos((source_pos(curr_iter,1)*rx + source_pos(curr_iter,2)*rz)./(norm(source_pos(curr_iter,:))*sqrt(rx.^2 + rz.^2))));
            r_norm = sqrt(rx.^2 + rz.^2);
            for n = -Nc:Nc
                uin(:,:,kk) = uin(:,:,kk) + A(1 + n + Nc,kk)*besselh(n,1,par.k*fact*r_norm).*exp(1i*n*curr_theta);
            end
            if nargout > 1
                varargout{1}(:,kk) = Hx - uinin(:,curr_iter);
            end
            %uin(:,:,kk) = conj(uin(:,:,kk));
            
            if par.doPlot
                figure(10);
                subplot(121);imagesc(abs(uin(:,:,kk)));axis image;colorbar;title('Reconstructed $\|Uin\|$');set(gca,'XDir','reverse');
                subplot(122);imagesc(angle(uin(:,:,kk)));axis image;colorbar;title('Reconstructed $<$Uin');set(gca,'XDir','reverse');
                pause(0.001);
            end
        end
end

end