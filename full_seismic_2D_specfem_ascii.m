dir=sprintf('%s/spe_606080_2/',getenv('SCRATCH'));

fnames={'rhoB_160','vp_160','vs_160'};
above_prop=[2e+3 3e+3 1.6e+3];
below_prop=[2e+3 3e+3 1.6e+3];
scale_from_files=1.0e+3;

n_cells_above = 120; %210
n_cells_insid = 80;
n_cells_below = 70; %110

nX = 60;
nY = 60;

out_matrix = zeros(nX, n_cells_above + n_cells_insid + n_cells_below, size(fnames(:), 1));

for i=1:size(fnames(:), 1)
    
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

    iy = 30; % middle
    out_matrix(:, 1:n_cells_above, i) = reshape(above(:, iy, :), size(above, 1), size(above, 3));
    out_matrix(:, n_cells_above+1:n_cells_above+n_cells_insid, i) = reshape(arr(:, iy, :), size(arr, 1), size(arr, 3));
    out_matrix(:, n_cells_above+n_cells_insid+1:n_cells_above+n_cells_insid+n_cells_below, i) = reshape(below(:, iy, :), size(below, 1), size(below, 3));

end

fname_out = sprintf('%s/spe2d_large_fine_specfem_rho_vp_vs.txt', dir);
disp(['Output file ', fname_out]);
fid = fopen(fname_out, 'w');
if (fid == -1)
    disp(['Cannot open file ', fname_out]);
    return;
end

for k = 1:size(out_matrix,2)
    for j = 1:size(out_matrix,1)
        % In each cell we have NGLLX x NGLLZ points, where NGLL? is
        % the number of points. For the 4-th SEM order, there are
        % 5 NGLLX and 5 NGLLZ
        for p = 1:25
            fprintf(fid, '0.0 0.0 %.12e %.12e %.12e\n', ...
                out_matrix(j, k, 1), out_matrix(j, k, 2), out_matrix(j, k, 3));
        end
    end
end

fclose(fid);

