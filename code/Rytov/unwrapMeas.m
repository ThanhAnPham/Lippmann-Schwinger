function [uMtot,loguMtot,uintot]  = unwrapMeas(uintot,uMtot,thetatot,doPlot,varargin)

for kk = 1:size(uMtot,2)
    uM = uMtot(:,kk);
    theta = thetatot(kk);
    curr_uin = uintot(:,kk);
    
    if round(1e3*theta)/1e3 >= 0
        %diff_phi = unwrap(angle(uM./curr_uin));
        phi_in = unwrap(angle(curr_uin),[],1);
        phi_tot = unwrap(angle(uM),[],1);
        diff_phi = phi_tot - phi_in;
    else %if round(1e3*theta)/1e3~=0
        %diff_phi = flipud(unwrap(flipud(angle(uM./curr_uin))));
        phi_in = flipud(unwrap(flipud(angle(curr_uin)),[],1));
        phi_tot =  flipud(unwrap(flipud(angle(uM)),[],1));
        diff_phi = phi_tot - phi_in;
        %diff_phi = angle(uM./curr_uin);
    end
    
    if doPlot
        figure(12);
        if nargin > 4
            subplot(121);
        end
        
        yyaxis left
        plot(angle(uM));hold on;
        title(num2str(rad2deg(theta)));
        plot(angle(curr_uin));
        plot(diff_phi,'LineWidth',2);hold off;
        yyaxis right
        plot(phi_tot); hold on;plot(phi_in);hold off;

        legend('wrapped u_M','wrapped uin','difference','unwrapped u_M','unwrapped uin');
        if nargin <=4
            pause();
        end
        if nargin > 4
            subplot(222);
            imagesc(angle(varargin{1}{1}(:,:,kk)));set(gca,'YDir','normal');
            subplot(224);
            imagesc(angle(varargin{1}{2}(:,:,kk)));set(gca,'YDir','normal');
            pause();
        end
    end
    
    uMtot(:,kk) = log(abs(uM./curr_uin)) + 1i*diff_phi;
    
end
%sometimes wrong starting value when unwrapped
shifted_uM = circshift(imag(uMtot(:,:)),-floor(size(uMtot,2)/2),2);
ind_corr = median(abs(uintot).*diff(shifted_uM(:,[1,1:end]),[],2));
unwrap_diff = circshift(cumsum(sign(ind_corr).*(abs((ind_corr)) > pi/2),2),floor(size(uMtot,2)/2),2);
loguMtot = real(uMtot) + 1i*(imag(uMtot) - repmat(2*pi*unwrap_diff,size(uMtot,1),1));
%+ repmat(round(unwrap_diff/(2*pi))*2*pi,size(uMtot,1),1));
if max(imag(loguMtot(:))) < 0
    fprintf('Unwrapped difference shifted by - 2 pi\n');
    loguMtot = real(loguMtot) + 1i*(imag(loguMtot) + 2*pi);
end
uMtot = uintot.*loguMtot;
end