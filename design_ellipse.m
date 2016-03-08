function [ell_mat]=design_ellipse(a, b, L)

if nargin==0
    a=2.5;
    b=3;
    L=30;
end

xlim=floor(a*sqrt(L));
x=[-xlim:xlim];
y=zeros(1,length(x)+2);
y(2:end-1)=round((b/a)* (a*a*L-x.*x).^0.5);

x2=[1 x+xlim+2 2*xlim+3];
% plot(x2,y)


ell_mat=repmat((1:max(y)+1)', [1 length(y)]);
y_mat=repmat(y, [max(y)+1, 1]);
ell_mat=ell_mat<=y_mat;

ell_mat2=flipud(ell_mat);
ell_mat=[ell_mat2; ell_mat];
% imtool(ell_mat)


