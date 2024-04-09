function [x,w,h,h2] = GetGRL(ngrl) 

%% load GRL points
name = sprintf('grl/grl_%0.3d.mat',ngrl);
if ~exist(name,'file')
  error(sprintf('Data file %s does not exist',name))
end
% disp(name)
load(name);
x=real(x);
w=real(w);
h=real(h);

if nargout>3
    h2=real(h2);
end
end

