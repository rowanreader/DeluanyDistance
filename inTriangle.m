function inside = inTriangle(p,a,b,c)
if sameSide(p,a,b,c) && sameSide(p,b,a,c) && sameSide(p,c,a,b)
    inside = 1;
else
    inside = 0;
end