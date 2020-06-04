classdef Vainikkotilde < LinOp
    % Vainikkotilde : Implements the 3D Lippmann-Schwinger Green function and its adjoints.
    % This operator gets the scattered field at the sensor position SenPos
    % 
    % :param k: Wavenumber (scalar, 2*pi/lambda0*n0, not dimensionless)
    % :param N: # of voxels in the volume (scalar or (1,3)/(3,1) vector
    % :param SenPos: Sensor Positions (structure with fields x,y,z)
    % :param dN: physical size of voxel (isotropic)
    % :param ups: ratio of size of sensor pixel/voxel (>1, int), experimental
    %
    % All attributes of parent class :class:`LinOp` are inherited. 
    %
    % **Note**: This operator corresponds to the \textbf{\tilde{G}} of [1].
    %           Be careful, can be memory demanding depending on the
    %           sizein/sizeout.
    %
    % **References**
    %
    % [1] Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.,
    %     IEEE Transactions on Computational Imaging (2020), T-a. Pham, E. Soubies, 
    %     A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
    %
    % **Example** Gtild=Vainikkotilde(k,N,SenPos,dN,ups)
    %
    % See also :class:`Map` :class:`LinOp`
  
    %%    Copyright (C) 
    %     2020 Thanh-an Pham thanh-an.pham@epfl.ch
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
    
    
    properties
        k;%2*pi/lambda*n0
        dN;%step size (e.g. nm / pixel), for now MUST be isotropic
        %Sensors position wrt sample center (ROI center), Nsensor x Ndims
        %UNIT : distance unit (e.g. microns, nm) converted in discrete unit
        %       by dividing with dN (e.g. nm / (nm / pixel) = pixel
        N;
        H;
        P;
        M;
        S;
        ups;
        ind;
        Zsen;
        fname;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function this = Vainikkotilde(k,N,SenPos,dN,ups)
            this.name = 'Vainikkotilde';
            
            this.sizein = N.*ones(1,3);
            N = N.*ones(1,3);
            
            this.dN = dN;
            this.k = k;
            this.N = N;
            this.ups = ups;
            this.Zsen = unique(SenPos.z);% - 0.5*dN(3);
            XYsenMin = [min(SenPos.x/(dN(1))),min(SenPos.y/(dN(2)))];
            XYsenMax = [max(SenPos.x/(dN(1))),max(SenPos.y/(dN(2)))];
            this.sizeout = floor(single(XYsenMax - XYsenMin)/ups(1)) + 1;
            Next = N(1:2) + this.sizeout;
            sizPad = [Next,N(3)];
            
            fname = sprintf('CrpGreenFunction_Z_%1.3e_dz_%1.3e_k_%1.3e_size_%i_%i_%i_%i_%i.mat',this.Zsen,this.dN(1),this.k,this.sizein,this.sizeout);
            fprintf('%s\n',fname);
            if exist(fname,'file')
                loadM = load(fname);
                this.H = flip(Sfft(ifftshift(ifftshift(loadM.tmp,1),2),3),3);
            else
                if 1%(this.sizeout(1)*4)^3*16 <= 30e9%3*this.sizeout(1)*16*4 <= 8e9
                    fprintf('Computing Gtild\n');
                    tmp = FreeSpaceKernel3D(k*dN(1),...
                        max(max(max(this.sizeout),this.sizein),2*this.Zsen/dN(3)),'Vainikko');
                    tmp = fftshift(iSfft(tmp.mtf,[]));
                    t = this.sizein/2;
                    t0 = this.Zsen/dN(end);
                    tmp = tmp(1 + end/2 - this.sizeout(1)/2 - t(1):end/2 + this.sizeout(1)/2 + t(1),...
                        1 + end/2 - this.sizeout(2)/2 - t(2):end/2 + this.sizeout(2)/2 + t(2),...
                        1 + end/2 + t0 - t(3):end/2 + t0 + t(3));
                    save(fname,'tmp','-v7.3');
                    this.H = flip(Sfft(ifftshift(ifftshift(tmp,1),2),3),3);
                else
                    fprintf('Check that the matrix was computed on the elephant\n');
                    %rethrow(ME);
                end
            end
            this.M = LinOpSelectorPatch(Next, round(Next/2 + XYsenMin/ups(1)),...
                    round(Next/2 + XYsenMax/ups(1)));
            this.P = LinOpSelectorPatch(sizPad,1 + [sizPad(1:2)/2 - this.sizein(1:2)/2 , 0],...
                    [sizPad(1:2)/2  - this.sizein(1:2)/2, 0] + this.sizein)';
            this.S = LinOpSum([Next,N(3)],3);
            this.fname = fname;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Access = protected) % all the other underscore methods
        function out = apply_(this,x)
            
            out = this.M*iSfft(this.S*(this.H.*Sfft(this.P*x,3)),[]);
            
        end
        function out = applyAdjoint_(this,x)
            out = this.P'*iSfft(conj(this.H).*(this.S'*Sfft(this.M'*x,[])),3);
        end
        
    end
    
end