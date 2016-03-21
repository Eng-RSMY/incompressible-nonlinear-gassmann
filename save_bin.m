function save_bin(tp,name,S)
fname=sprintf('output/%s_%i',name,tp);
fid=fopen(fname,'w');
if (fid==-1) disp(['Cannot open file ', fname]); return; end
fwrite(fid,S,'single');
fclose(fid);

