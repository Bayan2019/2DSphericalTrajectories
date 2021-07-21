% exponential map on unit sphere
function newp=Exp_Sphere(p,v)

d=norm(v,'fro');
if d>0.000000000001
    newp=cos(d)*p+sin(d)*v/d;
else
    newp=p;
end