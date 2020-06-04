function compareFields(u_gt,y,varargin)
%% COMPAREFIELDS Compare fields u_gt vs y (magnitude, phase, real, imag)

if nargin > 2
    curr_f = varargin{1};
else
    curr_f = figure;
end
subplot(531);imagesc(abs(u_gt));title('abs predicted');
subplot(532);imagesc(abs(y));title('abs measured');
subplot(533);imagesc(abs(u_gt) - abs(y));title('abs diff');
subplot(534);imagesc(angle(u_gt));title('phase predicted');
subplot(535);imagesc(angle(y));title('phase measured');
subplot(536);imagesc(angle(u_gt./y));title('phase diff');
subplot(537);imagesc(real(u_gt));title('real predicted');
subplot(538);imagesc(real(y));title('real measured');
subplot(539);imagesc(real(u_gt) - real(y));title('real diff');
subplot(5,3,10);imagesc(imag(u_gt));title('imag predicted');
subplot(5,3,11);imagesc(imag(y));title('imag measured');
subplot(5,3,12);imagesc(imag(u_gt) - imag(y));title('imag diff');
subplot(5,3,13);imagesc(real(u_gt - y));title('real (predicted - measured)');
subplot(5,3,14);imagesc(imag(u_gt - y));title('imag (predicted - measured)');
subplot(5,3,15);imagesc(abs(u_gt - y));title('abs (predicted - measured)');

h = findall(curr_f,'type', 'axes');
drawnow
linkaxes(h,'xy');

for kk = 1:length(h)
    colorbar(h(kk));
    if size(u_gt,1)/size(u_gt,2) < 20
        axis(h(kk),'image');
    end
end
end

