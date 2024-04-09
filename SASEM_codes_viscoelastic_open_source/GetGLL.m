function [x,w,h]=GetGLL(ngll,kind)

if nargin>1, 
  prefix=kind(1:3);
else
  prefix = 'gll';
end

name = sprintf('gll/%s_%0.2u.tab',prefix,ngll);
if ~exist(name,'file')
  error(sprintf('Data file %s does not exist',name))
end

fid=fopen(name);
data=fscanf(fid,'%f',[ngll,ngll+2]);
fclose(fid);

x=data(:,1);
w=data(:,2);
h=data(:,3:end)';
