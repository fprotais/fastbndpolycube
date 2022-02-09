fid=fopen('build/matrix_ref.mm');
e=textscan(fid,'%f %f %f','headerlines',1);
I = e{1};
J = e{2};
V = e{3};
S = sparse(I,J,V);
figure
spy(S)

