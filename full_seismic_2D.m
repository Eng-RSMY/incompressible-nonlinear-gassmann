dir=sprintf('%s/spe_606080_2/',getenv('SCRATCH'));

fnames={'rhoB_160','vp_160','vs_160'};
above_prop=[2e+3 3e+3 1.6e+3];
below_prop=[2e+3 3e+3 1.6e+3];
scale_from_files=1.0e+3;

for i=1:size(fnames(:))
    
    fname=sprintf('%s/%s', dir, fnames{i});
    disp(['Input file ', fname]);
    fid=fopen(fname, 'r');
    if (fid == -1)
        disp(['Cannot open file ', fname]);
        return;
    end
    arr=fread(fid,[60*60*80 1],'single');
    arr=scale_from_files*arr;
    arr=reshape(arr,60,60,80);
    fclose(fid);

    above=above_prop(i)*ones(60,60,120); %210);
    below=below_prop(i)*ones(60,60,70); %110);

    fname_full=sprintf('%s_fine_2D_specfem4_new',fname);
    disp(['Output file ', fname_full]);
    fid=fopen(fname_full, 'w');
    if (fid == -1)
        disp(['Cannot open file ', fname_full]);
        return;
    end
    
    iy = 30; % middle
    out_matrix = zeros(size(arr, 1), size(above, 3) + size(arr, 3) + size(below, 3));
    out_matrix(:, 1:size(above, 3)) = reshape(above(:, iy, :), size(above, 1), size(above, 3));
    out_matrix(:, size(above, 3)+1:size(above, 3)+size(arr, 3)) = reshape(arr(:, iy, :), size(arr, 1), size(arr, 3));
    out_matrix(:, size(above, 3)+size(arr, 3)+1:size(above, 3)+size(arr, 3)+size(below, 3)) = reshape(below(:, iy, :), size(below, 1), size(below, 3));
    
    disp(['Output array of size: ', num2str(size(out_matrix))]);
    
%    for k = 1:size(out_matrix,2)
%        for j = 1:size(out_matrix,1)
%            % In each cell we have NGLLX x NGLLZ points, where NGLL? is
%            % the number of points. For the 4-th SEM order, there are
%            % 5 NGLLX and 5 NGLLZ
%            for p = 1:25
%                fwrite(fid, out_matrix(j,k), 'single');
%            end
%        end
%    end
    
%    fwrite(fid, out_matrix, 'single');

%     fwrite(fid,above(:,iy,:),'single');
%     fwrite(fid,  arr(:,iy,:),'single');
%     fwrite(fid,below(:,iy,:),'single');

    fclose(fid);
    
    
%     rot_out_matrix = rot90(rot90(out_matrix));
%     
%     fname_full=sprintf('%s_fine_2D_rotated',fname);
%     disp(['Output file ', fname_full]);
%     fid=fopen(fname_full, 'w');
%     if (fid == -1)
%         disp(['Cannot open file ', fname_full]);
%         return;
%     end
%     
%     fwrite(fid, rot_out_matrix, 'single');
% 
%     fclose(fid);

end

