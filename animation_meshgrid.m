filename = 'outp_newSpin.txt';
delimiterIn = ' ';
X_trans = importdata(filename, delimiterIn);
X =  X_trans';
area = 1;

for i=1:size(X,1),
   for j=1:size(X,2),
      if(X(i,j)~= 0)
          %T(i,j) = X(i,j);
          T(i,j) = log(X(i,j)/area);
          if(T(i,j) < 0)
              T(i,j) = 0;
          end    
          %T(i,j) = log10(X(i,j));
      else
          T(i,j) = 0;
      end
      area = area + 2; 
   end
    
end
%---------------------------------------
r = 0.08;  % radius of circle

r_res = 0.0005;

rx = -r:r_res:r;
ry = rx;

[x, y] = meshgrid(rx, ry);

rxy2 = x.^2+y.^2;
z=ones(size(rxy2))*NaN;

%---------------------------------------

Nshells = size(T,2);
r = [0:1/Nshells:1]*r;
r2 = r.^2;

figure('Color', 'w');
%colormap hot
writerObj = VideoWriter('logOfDensity_GPU_g0.avi');
writerObj.FrameRate = 2.5;
open(writerObj);
for ind_t = 1:30
    for ii = 1:Nshells
        ir_find = find(rxy2<=r2(ii+1) & rxy2>r2(ii));
        z(ir_find) = T(ind_t,ii);
    end
    z(161,161) = T(ind_t,1);
    hax = axes('Position',[0 0 1 1]);
    h = surf(x,y,z);  % sphere centered at origin

    shading interp
    set(h, 'EdgeColor', 'None');
    %colorbar;
    view(0,90);
    axis equal;
    set(hax, 'Visible', 'Off', 'CLim', [min(T(:)) max(T(:))]);
    pause(0.4);
    f = getframe;
    writeVideo(writerObj, f);
end
close(writerObj);