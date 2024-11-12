function tf = maskSA(X,Y,C)
    [inep,onep] = inpolygon(X,Y,C{1}(:,1),C{1}(:,2));
    [inen,onen] = inpolygon(X,Y,C{2}(:,1),C{2}(:,2));
    tf = (inep & ~inen) | onep | onen;
end