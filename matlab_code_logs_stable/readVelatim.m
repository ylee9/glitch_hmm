fid=fopen([dname d(1).name]);
x=[];
s=fgets(fid);
while 1
	s=fgets(fid)
	if s==-1
		break
	end
	data = sscanf(s,'%s%g%g%g%s');
	n=find(data>5e4);
	x=[x; data(n) data(n+1)];
end
fclose(fid);
F0=11.1946499395;
dF=-1.5666E-11;
G=[57733.88938623312   57734.53626163834];
t0=51559.319;