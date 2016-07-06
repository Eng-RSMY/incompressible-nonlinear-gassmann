dir=sprintf('%s/spe_606080_2/',getenv('SCRATCH'));

fnames={'rhoB_160','vp_160','vs_160'};
above_prop=[2e+3 3e+3 1.6e+3];
below_prop=[2e+3 3e+3 1.6e+3];
scale_from_files=1.0e+3;

n_cells_above = 120; %210; %120; %210
n_cells_insid = 80;
n_cells_below = 70; %110; %70; %110

nX = 60;
nY = 60;

for i=1:size(fnames(:))
    
    fname=sprintf('%s/%s', dir, fnames{i});
    disp(['Input file ', fname]);
    fid=fopen(fname, 'r');
    if (fid == -1)
        disp(['Cannot open file ', fname]);
        return;
    end
    arr=fread(fid,[nX*nY*n_cells_insid 1],'single');
    arr=scale_from_files*arr;
    arr=reshape(arr,nX,nY,n_cells_insid);
    fclose(fid);

    above=above_prop(i)*ones(nX,nY,n_cells_above);
    below=below_prop(i)*ones(nX,nY,n_cells_below);

    fname_full=sprintf('%s_smooth_3D_f10',fname);
    disp(['Output file ', fname_full]);
    fid=fopen(fname_full, 'w');
    if (fid == -1)
        disp(['Cannot open file ', fname_full]);
        return;
    end

    tic;
    for iy=1:nY
        fwrite(fid,above(:,iy,:),'single');
        fwrite(fid,  arr(:,iy,:),'single');
        fwrite(fid,below(:,iy,:),'single');
    end
    
    toc;
    fclose(fid);

end

