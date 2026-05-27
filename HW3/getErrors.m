function [glb_errL2 ,glb_errH1] = getErrors(P,T,C,fun,Area)
errL2 = 0;
errH1 = 0;
for e = 1:size(T,1)
    nodes = T(e,:);
    Pe = P(nodes,:);
    Ce = C(nodes);

    x1 = Pe(1,1); y1 = Pe(1,2);
    x2 = Pe(2,1); y2 = Pe(2,2);
    x3 = Pe(3,1); y3 = Pe(3,2);
    area = Area(e);

    %quadrature points
    a1 = [x1,y1]; a2 = [x2,y2]; a3 = [x3,y3];
    a12 = (a1+a2)/2; a13 = (a1+a3)/2; a23 = (a2 + a3)/2;
    a123 = (a1 + a2 + a3)/3;
    
    quadpts = [a1;a2;a3;a12;a13;a23;a123];
    quadwts = [3;3;3;8;8;8;27];

    grad_phi = 1/(2*area)*[y2-y3 x3-x2;y3-y1 x1-x3;y1-y2 x2-x1];
    grad_uh  = Ce' * grad_phi;

    Ele_errL2 = 0;
    Ele_errH1 = 0;
    for q = 1:7
        xq = quadpts(q,1);
        yq = quadpts(q,2);
        
        u  = fun.UTRUE(xq,yq);
        ux = fun.UTRUEX(xq,yq);
        uy = fun.UTRUEY(xq,yq);
        uh = UhElement(Pe,Ce,[xq,yq],area);
        
        Ele_errL2 = Ele_errL2 + quadwts(q) * (uh - u)^2;
        Ele_errH1 = Ele_errH1 + quadwts(q) ...
             * ((grad_uh(1)-ux)^2 + (grad_uh(2)-uy)^2);
    end 
errL2 = errL2 + area/60 * Ele_errL2;
errH1 = errH1 + area/60 * Ele_errH1;
end
glb_errL2 = sqrt(errL2);
glb_errH1 = sqrt(errH1 + errL2);
return
end