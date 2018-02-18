%load objective function values and compute gradient

p = load('piVals');

g = (p(1:length(p)-1)-p(length(p)))/1E-7;

fid = fopen('gradVals','w');

for ii=1:length(g)
    fprintf(fid,'%1.5e\n',g(ii));
end

fclose(fid);
