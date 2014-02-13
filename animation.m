
Rmax = 0;
r = [];
for i=0:200,
    Rmax = Rmax + 10;
    r(end + 1) = Rmax;    
end

[x,y,z] = sphere(100);
x=x*Rmax;
y=y*Rmax;

rxy2 = x.^2+y.^2;


r2 = r.^2;

figure('Color', 'w');
%colormap('hot');
filename = 'outp_newSpin.txt';
delimiterIn = ' ';
X_trans = importdata(filename, delimiterIn);
X =  X_trans';

for i=1:size(X,1),
   for j=1:size(X,2),
      if(X(i,j)~= 0)
          T(i,j) = log10(X(i,j)/(i*2*pi*0.0125));
      else
          T(i,j) = 0;
      end
       
   end
    
end

%title('Visual Representation of Photon Exit Simulated with the GPU');

for ind_t = 1:size(T,1)
    for ii = 1:length(r2)-1
        ir_find = find(rxy2<=r2(ii+1) & rxy2>r2(ii));
        z(ir_find) = T(ind_t,ii);
    end

    hax = axes('Position',[0 0 1 1]);
   
    h = surf(x,y,z); % sphere centered at origin
    colorbar;
    
    shading interp
    set(h, 'EdgeColor', 'None');
  
    view(0,90);
    axis equal;
    %set(hax, 'Visible', 'Off', 'CLim', [min(T(:)) max(T(:))]);
    set(hax, 'Visible', 'Off', 'CLim', [1 5]);
    pause(0.15);
end