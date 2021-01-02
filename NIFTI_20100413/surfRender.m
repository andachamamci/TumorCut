% visualize the surface
%function surfRender(a,isovalue,surfColor,transparencyVal)
%Ex: surfRender(A,0.5,'green');
function surfRender(a,isovalue,surfColor,transparencyVal)

if (nargin < 4)  transparencyVal = 1; end
if (nargin < 3)  surfColor = 'green'; end

p = patch(isosurface(a,isovalue));

[xsize,ysize,zsize] = size(a);

isonormals(a, p)

set(p, 'FaceColor', surfColor, 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
camlight; lighting phong
axis([1 xsize 1 ysize 1 zsize]);
alpha(transparencyVal);