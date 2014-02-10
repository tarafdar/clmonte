%T = [0 0 0 0 0;
%T = [0   0   0   300   0 0 0 0 0 0];
%T = [300   105   110   118   128 100 150 150 100 0];
%    109  110   117   124   134 100 150 130 1000 200;
%    114  118   120   130   138 100 1000 150 120 200];
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

filename = 'ece496/outp_newSpin.txt';
delimiterIn = ' ';
T_trans = importdata(filename, delimiterIn);
X =  T_trans';
% T =  T_trans';
for i=1:size(X,1),
   for j=1:size(X,2),
      if(X(i,j)~= 0)
          T(i,j) = X(i,j);
      else
          T(i,j) = 0;
      end
       
   end
    
end

colormap('hot');
%C = real2rgb(log(T), 'jet');
%colormap jet;
%T = log(T);
for ind_t = 1:size(T,1)
    for ii = 1:length(r2)-1
        ir_find = find(rxy2<=r2(ii+1) & rxy2>r2(ii));
        z(ir_find) = T(ind_t,ii);
    end

    hax = axes('Position',[0 0 1 1]);
    %t = real2grb(log(z), 'jet');
    %h = surf(x,y,z, t); % sphere centered at origin
    h = surf(x,y,z); % sphere centered at origin
    colorbar;
    %colormap jet;
    shading interp
    set(h, 'EdgeColor', 'None');
  
    view(0,90);
    axis equal;
    set(hax, 'Visible', 'Off', 'CLim', [min(T(:)) max(T(:))]);
    pause(0.4);
end