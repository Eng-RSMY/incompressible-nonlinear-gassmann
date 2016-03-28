dir=sprintf('%s/spe_z85/',getenv('SCRATCH'));

fnames={'rhoB_0','vp_0','vs_0'};
above_prop=[2.2e+3 4.2e+3 2.3e+3];
below_prop=[2.5e+3 5.9e+3 3.8e+3];
scale_from_files=1.0e+3;

for i=1:size(fnames(:))
    
    fname=sprintf('%s/%s', dir, fnames{i});
    disp(['Input file ', fname]);
    fid=fopen(fname, 'r');
    if (fid == -1)
        disp(['Cannot open file ', fname]);
        return;
    end
    arr=fread(fid,[60*220*80 1],'single');
    arr=scale_from_files*arr;
    arr=reshape(arr,60,220,80);
    fclose(fid);

    above=above_prop(i)*ones(60,220,100);
    below=below_prop(i)*ones(60,220,20);

    fname_full=sprintf('%s_full',fname);
    disp(['Output file ', fname_full]);
    fid=fopen(fname_full, 'w');
    if (fid == -1)
        disp(['Cannot open file ', fname_full]);
        return;
    end

    for iy=1:size(arr,2)
        fwrite(fid,above(:,iy,:),'single');
        fwrite(fid,  arr(:,iy,:),'single');
        fwrite(fid,below(:,iy,:),'single');
    end
    fclose(fid);

end

