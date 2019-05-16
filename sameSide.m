function same = sameSide(p1, p2, a,b)
cp1 = cross(b-a, p1-a);
cp2 = cross(b-a,p2-a);
if dot(cp1, cp2) > 0
    same = 1;
else
    same = 0;
end