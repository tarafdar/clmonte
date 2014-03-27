
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

%x = (1:1:size(T,2));
%y = (1:1:size(T,1));
%h= surf(x, y, T, gradient(T));
Tnew = T(1:40, 1:40);
[x,y] = meshgrid([0:1:39]);
%hax = axes('Position',[0 0 1 1]);
h= surf(x, y, Tnew, Tnew);
%shading interp
set(h, 'EdgeColor', 'None');
%axis equal;
%set(hax, 'Visible', 'Off', 'CLim', [min(T(:)) max(T(:))]);
colorbar;
% [x,y] = meshgrid([-2:.2:2]);
% Z = x.*exp(-x.^2-y.^2);
% surf(x,y,Z,gradient(Z))
% colorbar


