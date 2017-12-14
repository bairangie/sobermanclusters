function f = readbinfile(filename);

fid = fopen(filename,'r');
c = fread(fid,69); %skip 69 bytes header
num = fread(fid,1,'int32'); %read number of molecule in frame 0
A = fread(fid,18*num,'float32');
fclose(fid);

mol.x = A(1:18:end);
mol.y = A(2:18:end);
mol.xc = A(3:18:end);
mol.yc = A(4:18:end);
mol.h = A(5:18:end);
mol.area = A(6:18:end);
mol.width = A(7:18:end);
mol.phi = A(8:18:end);
mol.Ax = A(9:18:end);
mol.bg = A(10:18:end);
mol.I = A(11:18:end);
mol.cat = floor(A(12:18:end)/1.401e-45);
mol.valid = floor(A(13:18:end)/1.401e-45);
mol.frame = floor(A(14:18:end)/1.401e-45);
mol.length = floor(A(15:18:end)/1.401e-45);
mol.link = floor(A(16:18:end)/1.401e-45);
mol.z = A(17:18:num*18);
mol.zc = A(18:18:num*18);

clear A;

f = mol;