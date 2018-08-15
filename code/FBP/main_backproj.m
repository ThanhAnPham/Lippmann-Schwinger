%% Back projection

if exist('log_imag_ur2D','var')
    shift = hasandis(simP,'shift');
    curr_log_imag_ur2D = (-1)^par.negative*circshift(log_imag_ur2D(:,par.curr_thetas),[shift,0]);
    n_hat.fbp = circshift((-1)^par.negative*FBProj(curr_log_imag_ur2D,rad2deg(par.thetas))'/(simP.k0*simP.dx) + simP.n0,[-shift,0]);
    par.fbp.Lx = simP.Nx*simP.dx; par.fbp.Lz = par.fbp.Lx;
    par.fbp.dx = simP.dx; par.fbp.dz = par.fbp.dx;
    
    figure(13);
    imagesc(simP.Lx_vec,simP.Lz_vec,n_hat.fbp);axis image;colorbar;
    if par.negative
        colormap(flipud(gray));
    else
        colormap gray;
    end
    dispGTover(simP,n);
else
    main_FBP;
end
