
P = [];


Pk = [xk yk dr1dy];  
P = [P; Pk]; 

% generate elememnt connectivity
nloc = [];

for k = 1:n
    loc = [(k-1)*2+1 (k-1)*2+2 (k-1)*2+3];
    nloc = [nloc; loc];
end
