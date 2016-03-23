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

above=(2.2e+3)*ones(60,220,86);
below=(2.2e+3)*ones(60,220,84);

fname='data/rhoB_2_full';
fid=fopen(fname, 'w');
if (fid == -1)
    disp('Cannot open file');
    return;
end

% we take an XY layer and put it as an XZ layer


for iy=1:size(rhoB,2)
    for iz=1:size(above,3)
       for ix=1:size(above,1)
           fwrite(fid,above(ix,iy,iz),'single');
       end
    end
    for iz=1:size(rhoB,3)
        for ix=1:size(rhoB,1)
            fwrite(fid, rhoB(ix,iy,iz),'single');
        end
    end
    for iz=1:size(below,3)
       for ix=1:size(below,1)
           fwrite(fid,below(ix,iy,iz),'single');
       end
    end
end
fclose(fid);

