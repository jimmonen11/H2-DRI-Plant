function scalex = scalefun(a, b, x)

maxx = max(x);
minx = min(x);

scalex = (b-a)*((x-minx)/(maxx-minx)) + a;
