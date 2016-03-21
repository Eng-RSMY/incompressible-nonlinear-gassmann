fname='data/rhoB_2';
fid=fopen(fname, 'r');
if (fid == -1)
    disp('Cannot open file');
    return;
end
rhoB=fread(fid,[60*220*10 1],'single');
rhoB=(1.0e+3)*rhoB;
rhoB=reshape(rhoB,60,220,10);
fclose(fid);

above=(2.2e+3)*ones(60,166,220);
under=(2.2e+3)*ones(60,84,220);

fname='data/rhoB_2_full';
fid=fopen(fname, 'w');
if (fid == -1)
    disp('Cannot open file');
    return;
end
fwrite(fid,above,'single');
%fwrite(fid,rhoB,'single');
for iz=1:size(rhoB,3)
    for iy=1:size(rhoB,2)
        for ix=1:size(rhoB,1)
            fwrite(fid,rhoB(ix,iy,iz),'single');
        end
    end
end
fwrite(fid,under,'single');
fclose(fid);

