function ps = as_polygon(a,b,nPoly)
    % Returns the polygon chain that corresponds to the given state.
    
    as = (0:360/(nPoly-1):360) * pi / 180;
    ps = [a;b].*[cos(as); sin(as)];
end