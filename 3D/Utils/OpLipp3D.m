classdef OpLipp3D <  Map
    % OpLipp3D : Lippmann-Schwinger forward operator. Refer to [1] for details
    %
    % :param G: FreeSpaceKernel3D (Green function operator, $\mathbf{G}$)
    % :param Gtil: Vainikkotilde (Green function operator for Sensor, $\mathbf{\tilde{G}}$)
    % :param M: LinOpSelectorPatch/LinOpIdentity
    % :param uin: 3D incident field
    % :param Niter: # iterations to compute the forward model/Jacobian
    % :param xtol: convergence criteria to compute the forward model/Jacobian
    % :param uem: field with f=zeros(this.sizein) (~Incident field) on the sensor
    % :param totalfield: boolean for computing the total or scattered field
    %
    %
    % **References**
    %
    % [1] Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.,
    %     IEEE Transactions on Computational Imaging (2020), T-a. Pham, E. Soubies, 
    %     A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
    %
    % **Example** H=OpLipp3D(G,Gtil,M,uin,Niter,xtol,uem,totalfield)
    %
    % Please refer to the Map superclass for documentation
    %     Copyright (C)
    %      2020 T.-A. Pham thanh-an.pham@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = public,GetAccess = public)
        G;
        CG;
        MGtild;
        uin;
        Niter = 120;
        xtol = 1e-8;
        u;
        f;
        uem;
        totalfield = true;
        a = 1;
    end
    properties (SetAccess = protected, GetAccess = public)
        useLocalGPU = false;
    end
    
    methods
        
        function this = OpLipp3D(G,Gtil,M,uin,Niter,xtol,uem,totalfield)
            this.name = 'OpLipp3D';
            this.isDifferentiable = true;
            this.sizein = G.sizein;
            uem = gpuCpuConverter(uem);
            uin = gpuCpuConverter(uin);
            G.mtf = gpuCpuConverter(G.mtf);
            
            if isa(Gtil,'FreeSpaceKernel3Dtilde') || isa(Gtil,'FreeSpaceKernel3DtildeMat') || isa(Gtil,'Vainikkotilde')
                this.sizeout = Gtil.sizeout;%Gtil directly outputs the measurements
                this.MGtild = Gtil;
                this.uem = uem;%M must be incident field at the sensors if total field is true, (can be empty otherwise)
            elseif isa(Gtil,'Gtild3DFFT')
                this.sizeout = M.sizeout;
                Gtil.G.mtf = gpuCpuConverter(Gtil.G.mtf);
                this.MGtild = M*Gtil;%Gtil outputs scattered field on a volume that includes the sensors
                this.uem = uem;
            else
                error('Gtild not correct');
            end
            
            this.G = G;
            if ~isempty(uin)
                Crp = LinOpSelectorPatch(size(uin),...
                    1 + (size(uin) - this.sizein)/2,(size(uin) + this.sizein)/2);
                this.uin = Crp*uin;
            end
            this.Niter = Niter;
            this.xtol = xtol;
            this.totalfield = totalfield;
            this.u = this.uin;
        end
        function setGPU(this,newState)
            this.useLocalGPU = newState;
        end
    end
    
    methods (Access = protected)
        function [y, utot] = apply_(this,x)
            
            iuin = this.a*this.uin;
            if this.useLocalGPU
                iuin = gpuArray(iuin);
                x = gpuArray(x);
                this.G.mtf = gpuArray(this.G.mtf);
                useGPU(1);
            end
            
            forw = LippmannUi3D(this.G, x);
            
            
            this.CG = OptiBCGS(forw,iuin);
            %this.CG.OutOp = OutputOptiSNR(1,[],1);
            this.CG.CvOp = TestCvgStepRelative(this.xtol);
            this.CG.ItUpOut = 0; this.CG.maxiter = this.Niter; this.CG.verbose = 0;
            
            this.CG.run(iuin);
            utot = this.CG.xopt;
            
            %fprintf('Niter in forward : %i. ',this.CG.niter);
            
            if any(isnan(utot(:))) %because of CG not stopping at the beginning
                fprintf('NaN in CG (if empty volume, expected)\n');
                utot = iuin;
            end
            
            this.u = utot;
            y = this.MGtild*(utot.*x);
            if this.totalfield
                y = y + this.uem;
            end
            %clear iuin x utot;%shouldn't be needed but never sure.
            if this.useLocalGPU
                y = gather(y);
                this.G.mtf = gather(this.G.mtf);
                this.u = gather(this.u);
                useGPU(0);
            end
            this.CG = [];
        end
        
        function out = applyJacobianT_(this,x,f)
            if this.useLocalGPU
                useGPU(1);
                x = gpuArray(x);
                f = gpuArray(f);
                this.G.mtf = gpuArray(this.G.mtf);
                this.u = gpuArray(this.u);
                this.MGtild.H = gpuArray(this.MGtild.H);
                
            end
            forw = LippmannUi3D(this.G,f);
            GtMtx = this.MGtild'*x;
            CGgrad = OptiBCGS(forw', f.*GtMtx);
            CGgrad.ItUpOut = 1;
            CGgrad.maxiter = this.Niter; CGgrad.CvOp = TestCvgStepRelative(this.xtol);
            CGgrad.verbose = false;
            
            CGgrad.run(CGgrad.A'*CGgrad.b);
            
            
            if ~any(isnan(CGgrad.xopt(:)))
                out = real(conj(this.u).*(this.G'*CGgrad.xopt + GtMtx));
            else
                out = real(conj(this.u).*GtMtx);
            end
            this.u = [];
            if this.useLocalGPU
                clear x f
                this.G.mtf = gather(this.G.mtf);
                this.u = gather(this.u);
                this.MGtild.H = gather(this.MGtild.H);
                out = gather(out);
                useGPU(0);
            end
        end
        
    end
end
