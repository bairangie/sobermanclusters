% read molecule list from bin or txt file

[readfile,path1] = uigetfile('*.bin','Read  Bin File');
filename=strcat(path1,readfile);
r = readbinfile(filename);
pixel_size = 160; % in nm
x_size = 256+10; % 10 extra pixel for margin of drift correction
y_size = 256+10;


%assign each localization to pixel
%count the number of localization in each pixel
bin = cell(x_size,y_size,2);
count = zeros(x_size,y_size,2);
for x=1:length(r.xc)
    x_position = ceil(r.xc(x));
    if x_position < 1
        x_position = 1;
    end
    y_position = ceil(r.yc(x));
    if y_position < 1
        y_position = 1;
    end
    if r.cat(x) == 1
        bin{x_position,y_position,1} = [bin{x_position,y_position,1} x];
        count(x_position,y_position,1) = count(x_position,y_position,1)+1;
    elseif r.cat(x) == 2
        bin{x_position,y_position,2} = [bin{x_position,y_position,2} x];
        count(x_position,y_position,2) = count(x_position,y_position,2)+1;
    end
end

figure(1)
imshow(count(:,:,1),colormap(parula(50)));
colorbar
outfile1=strcat(path1,'Fig1.png');
saveas(gcf,outfile1)

figure(2)
imshow(count(:,:,2),colormap(hot(50)));
colorbar
outfile2=strcat(path1,'Fig2.png');
saveas(gcf,outfile2)