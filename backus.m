dir=sprintf('%s/spe_z85/',getenv('SCRATCH'));

rhoB_file = sprintf('%s/rhoB_0', dir)
vp_file   = sprintf('%s/vp_0', dir)
vs_file   = sprintf('%s/vs_0', dir)

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

rhoB = fread(rhoB_fid,[60*220*80 1],'single');
rhoB = scale_from_files*rhoB;
rhoB = reshape(rhoB, 60, 220, 80);
fclose(rhoB_fid);

vp = fread(vp_fid,[60*220*80 1],'single');
vp = scale_from_files*vp;
vp = reshape(vp, 60, 220, 80);
fclose(vp_fid);

vs = fread(vs_fid,[60*220*80 1],'single');
vs = scale_from_files*vs;
vs = reshape(vs, 60, 220, 80);
fclose(vs_fid);

lambda = rhoB .* (vp.^2 - vs.^2);
mu     = rhoB .* (vs.^2);

%    for iy=1,size(vs,2)
%        for ix=1,size(vs,1)
%            mu_aver = zeros(size(vs,3)/10, 1);
%            for ii=1,size(vs,3)/10
%                for jj=1,10
%                    mu_aver(ii) = mu_aver(ii) + mu(ix, iy, (ii-1)*10+jj);
%                end
%            end
%        end
%    end

%    lam_aver = zeros(size(vs,1), size(vs,2), size(vs,3)/10);
%    mu_aver  = zeros(size(vs,1), size(vs,2), size(vs,3)/10);
%    for ii=1,size(vs,3)/10
%        for jj=1,10
%            lam_aver(:,:,ii) = lam_aver(:,:,ii) + lambda(:, :, (ii-1)*10+jj);
%            mu_aver(:,:,ii)  = mu_aver(:,:,ii)  + mu(:, :, (ii-1)*10+jj);
%        end
%    end

C = zeros(size(vs,1), size(vs,2), size(vs,3)/10);
L = zeros(size(vs,1), size(vs,2), size(vs,3)/10);
R = zeros(size(vs,1), size(vs,2), size(vs,3)/10);
for ii=1,size(vs,3)/10
    for jj=1,10
        C(:,:,ii) = C(:,:,ii) + ones(size(vs,1), size(vs,2)) ./ (lambda(:, :, (ii-1)*10+jj) + 2.0*mu(:, :, (ii-1)*10+jj));
        L(:,:,ii) = L(:,:,ii) + ones(size(vs,1), size(vs,2)) ./ mu(:, :, (ii-1)*10+jj);
        R(:,:,ii) = R(:,:,ii) + rhoB(:, :, (ii-1)*10+jj);
    end
    C(:,:,ii) = C(:,:,ii) / 10.0;
    L(:,:,ii) = L(:,:,ii) / 10.0;
    R(:,:,ii) = R(:,:,ii) / 10.0;
 end

C = ones(size(C,1), size(C,2), size(C,3)) ./ C;
L = ones(size(L,1), size(L,2), size(L,3)) ./ L;
    
vp_aver = sqrt(C ./ R);
vs_aver = sqrt(L ./ R);


rhoB_av_file = sprintf('%s_aver', rhoB_file)
vp_av_file   = sprintf('%s_aver', vp_file)
vs_av_file   = sprintf('%s_aver', vs_file)


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


for iy=1:size(vs,2)
    fwrite(rhoB_fid, R(:,iy,:), 'single');
    fwrite(vp_fid, vp_aver(:,iy,:), 'single');
    fwrite(vs_fid, vs_aver(:,iy,:), 'single');
end

fclose(rhoB_fid);
fclose(vp_fid);
fclose(vs_fid);

