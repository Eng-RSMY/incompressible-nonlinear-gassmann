function save_bin(tp,name,S)
dir=getenv('SCRATCH');
fname=sprintf('%s/spe_606080_2/%s_%i',dir,name,tp);
fid=fopen(fname,'w');
if (fid==-1) disp(['Cannot open file ', fname]); return; end
fwrite(fid,S,'single');
fclose(fid);

