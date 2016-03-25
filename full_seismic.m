dir=getenv('SCRATCH');

fnames=['rhoB_2' 'vp_2' 'vs_2'];
above_prop=[2.2e+3 4.2e+3 2.3e+3];
below_prop=[2.2e+3 4.2e+3 2.3e+3];
scale_from_files=1.0e+3;

for i=1:size(fnames(:))

    fname=sprintf('%s/%s', dir, fnames(i));
    disp(['Input file ', fname]);
    fid=fopen(fname, 'r');
    if (fid == -1)
        disp(['Cannot open file ', fname]);
        return;
    end
    arr=fread(fid,[60*220*10 1],'single');
    arr=scale_from_files*arr;
    arr=reshape(arr,60,220,10);
    fclose(fid);

    above=above_prop*ones(60,220,50);
    below=below_prop*ones(60,220,100);

    fname_full=sprintf('%s_full',fname);
    disp(['Output file ', fname_full]);
    fid=fopen(fname_full, 'w');
    if (fid == -1)
        disp(['Cannot open file ', fname_full]);
        return;
    end

    for iy=1:size(arr,2)
%        for iz=1:size(above,3)
%            for ix=1:size(above,1)
                fwrite(fid,above(:,iy,:),'single');
%            end
%        end
%        for iz=1:size(arr,3)
%            for ix=1:size(arr,1)
                fwrite(fid, arr(:,iy,:),'single');
%            end
%        end
%        for iz=1:size(below,3)
%            for ix=1:size(below,1)
                fwrite(fid,below(:,iy,:),'single');
%            end
%        end
   end
   fclose(fid);

end

