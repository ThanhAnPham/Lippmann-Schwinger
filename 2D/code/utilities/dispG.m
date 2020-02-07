function dispG(f_res,n_gt,uin,y,G,varargin)

if nargin > 5
    Gtild = varargin{1};
end

if ~isempty(G)
    t = 3;
else
    t = 2;
end
figure(4);clf
subplot(t,2,1);imagesc(f_res);
axis image;colorbar;title('Ground Truth f/dn : ROI');
%if max(f_res(:)) <= 0
%    tmp = colormap;
%    colormap(flipud(tmp));
%end
subplot(t,2,2);imagesc(n_gt);
axis image;colorbar;title('Ground Truth n: sensors included');
subplot(t,2,3);
if ismatrix(uin)
    imagesc(abs(uin));
else
    imagesc(abs(squeeze(uin(:,:,round(end/2)))));
end
colorbar;

subplot(t,2,4);imagesc(abs(y));colorbar;

if ~isempty(G)
    subplot(t,2,5);imagesc(log(1 + fftshift(abs(G.mtf))));
    axis image;title('Green function Fourier Transform G');
    if exist('Gtild','var')
        subplot(t,2,6);imagesc(log(1 + fftshift(abs(Gtild.mtf))));
        axis image;title('Green function Fourier Transform Gtild');
    end
end
drawnow
end