dir=sprintf('%s/spe_606080_2/',getenv('SCRATCH'));

rhoB_file = sprintf('%s/rhoB_160', dir)
vp_file   = sprintf('%s/vp_160', dir)
vs_file   = sprintf('%s/vs_160', dir)

scale_from_files=1.0e+3;

rhoB_fid = fopen(rhoB_file, 'r');
if (rhoB_fid == -1)
    disp(['Cannot open file ', rhoB_file]);
    return;
end
vp_fid = fopen(vp_file, 'r');
if (vp_fid == -1)
    disp(['Cannot open file ', vp_file]);
    return;
end
vs_fid = fopen(vs_file, 'r');
if (vs_fid == -1)
    disp(['Cannot open file ', vs_file]);
    return;
end

Nx=60;
Ny=60;
Nz=80;

rhoB = fread(rhoB_fid,[Nx*Ny*Nz 1],'single');
rhoB = scale_from_files*rhoB;
rhoB = reshape(rhoB, Nx, Ny, Nz);
fclose(rhoB_fid);

vp = fread(vp_fid,[Nx*Ny*Nz 1],'single');
vp = scale_from_files*vp;
vp = reshape(vp, Nx, Ny, Nz);
fclose(vp_fid);

vs = fread(vs_fid,[Nx*Ny*Nz 1],'single');
vs = scale_from_files*vs;
vs = reshape(vs, Nx, Ny, Nz);
fclose(vs_fid);

disp(['rhoB: min ', num2str(min(rhoB(:))), ', max ', num2str(max(rhoB(:)))]);
disp(['vp: min ', num2str(min(vp(:))), ', max ', num2str(max(vp(:)))]);
disp(['vs: min ', num2str(min(vs(:))), ', max ', num2str(max(vs(:)))]);

lambda = rhoB .* (vp.^2 - 2.0*vs.^2);
mu     = rhoB .* (vs.^2);

disp(['lambda: min ', num2str(min(lambda(:))), ', max ', num2str(max(lambda(:)))]);
disp(['mu: min ', num2str(min(mu(:))), ', max ', num2str(max(mu(:)))]);

N = 10;

C = zeros(size(vs,1), size(vs,2), size(vs,3)/N);
L = zeros(size(vs,1), size(vs,2), size(vs,3)/N);
R = zeros(size(vs,1), size(vs,2), size(vs,3)/N);
VP= zeros(size(vs,1), size(vs,2), size(vs,3)/N);
VS= zeros(size(vs,1), size(vs,2), size(vs,3)/N);
for ii=1:(size(vs,3)/N);
    for jj=1:N;
        C(:,:,ii) = C(:,:,ii) + ones(size(vs,1), size(vs,2)) ./ (lambda(:, :, (ii-1)*N+jj) + 2.0*mu(:, :, (ii-1)*N+jj));
        L(:,:,ii) = L(:,:,ii) + ones(size(vs,1), size(vs,2)) ./ mu(:, :, (ii-1)*N+jj);
        R(:,:,ii) = R(:,:,ii) + rhoB(:, :, (ii-1)*N+jj);
        VP(:,:,ii)= VP(:,:,ii)+ vp(:, :, (ii-1)*N+jj);
        VS(:,:,ii)= VS(:,:,ii)+ vs(:, :, (ii-1)*N+jj);
    end
    C(:,:,ii) = C(:,:,ii) / N;
    L(:,:,ii) = L(:,:,ii) / N;
    R(:,:,ii) = R(:,:,ii) / N;
    VP(:,:,ii)= VP(:,:,ii)/ N;
    VS(:,:,ii)= VS(:,:,ii)/ N;
end

C = ones(size(C,1), size(C,2), size(C,3)) ./ C;
L = ones(size(L,1), size(L,2), size(L,3)) ./ L;
    
vp_aver = sqrt(C ./ R);
vs_aver = sqrt(L ./ R);

disp(['rhoB_aver: min ', num2str(min(R(:))), ', max ', num2str(max(R(:)))]);
disp(['vp_aver: min ', num2str(min(vp_aver(:))), ', max ', num2str(max(vp_aver(:)))]);
disp(['vs_aver: min ', num2str(min(vs_aver(:))), ', max ', num2str(max(vs_aver(:)))]);
disp(['VP: min ', num2str(min(VP(:))), ', max ', num2str(max(VP(:)))]);
disp(['VS: min ', num2str(min(VS(:))), ', max ', num2str(max(VS(:)))]);


rhoB_orig_file = sprintf('%s_orig', rhoB_file)
vp_orig_file   = sprintf('%s_orig', vp_file)
vs_orig_file   = sprintf('%s_orig', vs_file)

rhoB_fid_orig = fopen(rhoB_orig_file, 'w');
if (rhoB_fid_orig == -1)
    disp(['Cannot open file ', rhoB_orig_file]);
    return;
end
vp_fid_orig = fopen(vp_orig_file, 'w');
if (vp_fid_orig == -1)
    disp(['Cannot open file ', vp_orig_file]);
    return;
end
vs_fid_orig = fopen(vs_orig_file, 'w');
if (vs_fid_orig == -1)
    disp(['Cannot open file ', vs_orig_file]);
    return;
end



rhoB_av_file = sprintf('%s_backus', rhoB_file)
vp_av_file   = sprintf('%s_backus', vp_file)
vs_av_file   = sprintf('%s_backus', vs_file)


rhoB_fid = fopen(rhoB_av_file, 'w');
if (rhoB_fid == -1)
    disp(['Cannot open file ', rhoB_av_file]);
    return;
end
vp_fid = fopen(vp_av_file, 'w');
if (vp_fid == -1)
    disp(['Cannot open file ', vp_av_file]);
    return;
end
vs_fid = fopen(vs_av_file, 'w');
if (vs_fid == -1)
    disp(['Cannot open file ', vs_av_file]);
    return;
end

Nz1=51;
Nz2=21;

rhoB_above = 2.0e+3*ones(Nx, Ny, Nz1);
rhoB_below = 2.0e+3*ones(Nx, Ny, Nz2);

vp_above = 3.0e+3*ones(Nx, Ny, Nz1);
vp_below = 3.0e+3*ones(Nx, Ny, Nz2);

vs_above = 1.6e+3*ones(Nx, Ny, Nz1);
vs_below = 1.6e+3*ones(Nx, Ny, Nz2);


for iy=1:size(vs,2)
    fwrite(rhoB_fid_orig, rhoB(:,iy,:), 'single');
    fwrite(vp_fid_orig, vp(:,iy,:), 'single');
    fwrite(vs_fid_orig, vs(:,iy,:), 'single');
        
    fwrite(rhoB_fid, rhoB_above(:,iy,:), 'single');
    fwrite(rhoB_fid, R(:,iy,:), 'single');
    fwrite(rhoB_fid, rhoB_below(:,iy,:), 'single');
    
    fwrite(vp_fid, vp_above(:,iy,:), 'single');
    fwrite(vp_fid, vp_aver(:,iy,:), 'single');
    fwrite(vp_fid, vp_below(:,iy,:), 'single');
    
    fwrite(vs_fid, vs_above(:,iy,:), 'single');
    fwrite(vs_fid, vs_aver(:,iy,:), 'single');
    fwrite(vs_fid, vs_below(:,iy,:), 'single');
end

fclose(rhoB_fid);
fclose(vp_fid);
fclose(vs_fid);

fclose(rhoB_fid_orig);
fclose(vp_fid_orig);
fclose(vs_fid_orig);
