classdef FreeSpaceKernel3D < LinOp
    % FreeSpaceKernel3D : Implements the 3D Helmholtz Green function (and its adjoints) in two different ways :
    %       1/ By arbitrarily cropping the 3D space kernel
    %       2/ Using the Vainikko filter in the Fourier domain [1]
    %
    % :param k: Wavenumber (scalar, 2*pi/lambda0*n0*dN, dimensionless)
    % :param N: # of voxels in the volume (scalar or (1,3)/(3,1) vector
    % :param ForwardModel: String 'CroppedGreen' or 'Vainikko'
    % :param p: Padding factor (default 4)
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: This operator corresponds to the \textbf{G} of [1].
    %           Be careful, can be memory demanding depending on the sizein
    %
    % **References**
    %
    % [1] Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.,
    %     IEEE Transactions on Computational Imaging (2020), T-a. Pham, E. Soubies, 
    %     A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
    %
    % **Example** G=FreeSpaceKernel3D(k,N,ForwardModel)
    %
    % See also :class:`Map` :class:`LinOp`
    
    properties
        k;
        mtf;
        N;
        ForwardModel;  %Chosen forward model : CroppedGreen or Vanikko
        Crp;
        NOvs;
        L;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function this = FreeSpaceKernel3D(k,N,ForwardModel,varargin)
            this.name = 'FreeSpaceKernel3D';
            if nargin < 4
                ForwardModel = 'Vainikko';
            end
            
            this.sizein = N.*ones(1,3);
            this.sizeout = this.sizein;
            L = norm(N.*ones(1,3));
            this.k = k;
            meth = 2 - strcmp(ForwardModel,'CroppedGreen');
            switch meth
                case 1 % À la Vico
                    if nargin > 3
                        p = varargin{1};
                    else
                        p = 4;
                    end
                    NOvs = p*N.*ones(1,3);
                    this.NOvs = NOvs;
                   
                    this.Crp = LinOpSelectorPatch(this.NOvs.*ones(1,3),ones(size(this.sizein)),this.sizein);
                    
                    this.ForwardModel = ForwardModel;
                    
                    Vect = (0:NOvs(1)-1) - NOvs(1)/2;
                    if true %memory economy
                        Rcoord = repmat(Vect(:).^2,[1,NOvs(2),NOvs(3)]);
                        Rcoord = Rcoord + repmat(reshape(Vect(:).^2,1,NOvs(2)),[NOvs(1),1,NOvs(3)]);
                        Rcoord = sqrt(Rcoord + repmat(reshape(Vect(:).^2,1,1,NOvs(3)),[NOvs(1),NOvs(2),1]));
                    else
                        [X,Y,Z] = meshgrid(Vect);
                        Rcoord = sqrt(X.^2 + Y.^2 + Z.^2);clear X Y Z
                    end
                    this.N = N;
                    
                    if strcmp(this.ForwardModel,'CroppedGreen')
                        thres = eps;
                        pos_0 = Rcoord <= thres; % Note : just to get the center
                        Rcoord(pos_0) = 1/(4*pi*thres);%arbitrary, but not too low, neither too high
                        otf = exp(1i*k*Rcoord)./(4*pi*Rcoord);
                        %this.otf(pos_0) = 1;%Fareol'd rather set 1 at the
                        %singularity, but doesn't work
                        this.mtf = gpuCpuConverter(fftn(ifftshift(otf)));%(0,0) coordinate at the top left
                        this.Crp = LinOpSelectorPatch(this.NOvs.*ones(1,3),ones(size(this.sizein)),this.sizein);
                    elseif strcmp(this.ForwardModel,'Vainikko')
                        
                        s = 2*pi*Rcoord/NOvs(1);
                        clear Rcoord;
                        
                        this.mtf = (-1 + exp(1i*L*k)*(cos(L*s)...
                            - 1i*k.*L.*sinc(L*s/pi)))./(k^2 - s.^2);%sinc in matlab got a pi,%- 1i*k.*sin(L*s)));%
                        
                        
                        pos0 = s.^2 ==  k^2;
                        fprintf('%i points at the singularity in Green function (this is taken care)\n',nnz(pos0(:)));
                        
                        this.mtf(pos0) = exp(1i*k*L)*((k^2*L - 1i*k)*sin(L*k) + 1i*k^2*L*cos(L*k))/k^3;
                        
                        this.L = L;
                        this.mtf = ifftshift(this.mtf);
                        
                        fprintf('Calculating the precomputated kernel for Green function by cropping\n');
                        
                        this.NOvs = 2*this.sizein;
                        
                        CrpT = LinOpSelectorPatch(NOvs.*ones(1,3),1 + NOvs/2 - this.NOvs/2,NOvs/2 + this.NOvs/2);
                        this.mtf = gpuCpuConverter(fftn(ifftshift(CrpT*fftshift(ifftn(this.mtf)))));
                        this.Crp = LinOpSelectorPatch(this.NOvs.*ones(1,3),ones(size(this.sizein)),this.sizein);
                    else
                        fprintf('Please give valid forward model choice\n\n');
                    end
                case 2 % See [1]
                    n = max(N);
                    if nargin > 3
                        p = varargin{1};
                    else
                        p = 4;
                    end
                    % -- Some precomputations
                    h=1;%unit domain (the physical size is included in k
                    %kb=2*pi*nb/lamb;
                    assert(mod(p,2)==0,'p must be a multiple of 2'); % (for the second way of constructing the kernel)
                    fourierStep=2/(2*h*n*p);
                    fr=fftshift(linspace(-1/(2*h),1/(2*h)-fourierStep,n*p));
                    gc=zeros([2*n,2*n,2*n]);
                    cc=[0:n-1,n*p-n:n*p-1];  % indexes in the cropped region (fftshift)
                    [x1,x2,x3]=meshgrid(cc,cc,cc);
                    nnzs = 0;
                    for s1=0:p/2-1
                        fr1=fr(s1+1:p/2:end);
                        for s2=0:p/2-1
                            fr2=fr(s2+1:p/2:end);
                            for s3=0:p/2-1
                                fr3=fr(s3+1:p/2:end);
                                if true
                                    normk = repmat(fr1(:).^2,[1,2*n,2*n]);
                                    normk = normk + repmat(reshape(fr2(:).^2,1,2*n),[2*n,1,2*n]);
                                    normk = normk + repmat(reshape(fr3(:).^2,1,1,2*n),[2*n,2*n,1]);
                                    normk = 2*pi*sqrt(normk);
                                else
                                    [k1,k2,k3]=meshgrid(fr1,fr2,fr3);
                                    normk=2*pi*sqrt(k1.^2+k2.^2+k3.^2);
                                end
                                
                                gchat=(-1 + exp(1i*L*k)*(cos(L*normk) - 1i * k * L *sinc(L*normk/pi)))./(k^2 - abs(normk).^2);
                                nnzs = nnzs + nnz(normk==k);
                                gchat(normk==k)=1i*(L/(2*k) -exp(1i*L*k)*sin(L*k)/(2*k^2));% singularity at kb
                                gc=gc + (fftshift(ifftn(gchat).*exp(2*1i*pi/(n*p)*(s1.*x1+s2.*x2+s3.*x3))/(p/2)^3));
                            end
                        end
                    end
                    fprintf('%i points at the singularity in Green function (this is taken care)\n',nnzs);
                    this.NOvs = 2*N.*ones(1,3);%2*N;
                    gc = gc(1 + end/2 - this.NOvs(1)/2:end/2 + this.NOvs(1)/2,...
                        1 + end/2 - this.NOvs(2)/2:end/2 + this.NOvs(2)/2,...
                        1 + end/2 - this.NOvs(3)/2:end/2 + this.NOvs(3)/2);
                    this.mtf = gpuCpuConverter(fftn(ifftshift(gc/h^2)));
                    this.Crp = LinOpSelectorPatch(this.NOvs.*ones(1,3),ones(size(this.sizein)),this.sizein);
            end
        end
    end
    methods (Access = protected)
        function out = apply_(this,x)
            out = this.Crp*ifftn(fftn(this.Crp'*(x)) .* this.mtf);
        end
        
        function out = applyAdjoint_(this,x)
            out = this.Crp*(ifftn(fftn(this.Crp'*(x)) .* conj(this.mtf)));
        end
    end
    
end