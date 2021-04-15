
% Franke function
% NOT FULLY IMPLEMENTED

classdef F2d < Function2d     

    methods(Static)
        function obj = F2d(), obj@Function2d(); end
        
        function f = F(x,y)
            f = 0.75.*exp(-0.25.*(9.*x-2).^2 - 0.25.*(9.*y-2).^2) + 0.75.*exp(-((9.*x+1).^2)./49 - ((9.*y+1).^2)./10) + ...
                    0.5.*exp(-0.25.*(9.*x-7).^2-0.25.*(9.*y-3).^2) -  0.2.*exp(-(9.*x-4).^2-(9.*y-7).^2);
        end
        
        function f = x1(x,y)
            f = (-9/4).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2).*( ...
  (-7)+9.*x)+(18/5).*exp(1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y) ...
  .^2).*((-4)+9.*x)+(-27/8).*exp(1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+( ...
  -81/4).*y.^2).*((-2)+9.*x)+(-27/98).*exp(1).^((-1/49).*(1+9.*x) ...
  .^2+(-1/10).*(1+9.*y).^2).*(1+9.*x);
        end
        
        function f = x2(x,y)
            f = (81/8).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2).*( ...
  47+(-126).*x+81.*x.^2)+(243/16).*exp(1).^((-2)+9.*x+(-81/4).*x.^2+ ...
  9.*y+(-81/4).*y.^2).*(2+(-36).*x+81.*x.^2)+(-162/5).*exp(1).^((-1) ...
  .*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*(31+(-144).*x+162.*x.^2)+( ...
  243/4802).*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*(( ...
  -47)+36.*x+162.*x.^2);
        end
        
        function f = x3(x,y)
            f = (-729/16).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2) ...
  .*((-7)+9.*x).*(43+(-126).*x+81.*x.^2)+(-2187/32).*exp(1).^((-2)+ ...
  9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*((-2)+9.*x).*((-2)+(-36).* ...
  x+81.*x.^2)+(-2187/117649).*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10) ...
  .*(1+9.*y).^2).*(1+9.*x).*((-145)+36.*x+162.*x.^2)+(2916/5).*exp( ...
  1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*((-116)+837.*x+( ...
  -1944).*x.^2+1458.*x.^3);  
        end
        
        function f = x4(x,y)
            f = (6561/32).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2) ...
  .*(1825+(-10836).*x+22842.*x.^2+(-20412).*x.^3+6561.*x.^4)+( ...
  19683/64).*exp(1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*( ...
  (-20)+144.*x+972.*x.^2+(-5832).*x.^3+6561.*x.^4)+(-26244/5).*exp( ...
  1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*(835+(-8352).*x+ ...
  30132.*x.^2+(-46656).*x.^3+26244.*x.^4)+(19683/5764801).*exp(1).^( ...
  (-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*(6619+(-10440).*x+( ...
  -45684).*x.^2+11664.*x.^3+26244.*x.^4);
        end
        
        
        function f = y1(x,y)
            f = (9/40).*(exp (1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2).* ...
  (30+(-90).*y)+16.*exp (1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y) ...
  .^2).*((-7)+9.*y)+(-15).*exp (1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+( ...
  -81/4).*y.^2).*((-2)+9.*y)+(-6).*exp (1).^((-1/49).*(1+9.*x).^2+( ...
  -1/10).*(1+9.*y).^2).*(1+9.*y));
        end
        
        function f = y2(x,y),
            f = (81/400).*(50.*exp (1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y) ...
  .^2).*(7+(-54).*y+81.*y.^2)+75.*exp (1).^((-2)+9.*x+(-81/4).*x.^2+ ...
  9.*y+(-81/4).*y.^2).*(2+(-36).*y+81.*y.^2)+12.*exp (1).^((-1/49).*( ...
  1+9.*x).^2+(-1/10).*(1+9.*y).^2).*((-4)+18.*y+81.*y.^2)+(-160).* ...
  exp (1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*(97+(-252).* ...
  y+162.*y.^2));  
        end
        
        function f = y3(x,y)
            f = (729/4000).*((-2250).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+( ...
  -3).*y).^2).*((-1)+3.*y).*(1+(-18).*y+27.*y.^2)+(-375).*exp(1).^(( ...
  -2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*((-2)+9.*y).*((-2)+( ...
  -36).*y+81.*y.^2)+(-24).*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10).*( ...
  1+9.*y).^2).*(1+9.*y).*((-14)+18.*y+81.*y.^2)+3200.*exp(1).^((-1) ...
  .*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*((-665)+2619.*y+(-3402).* ...
  y.^2+1458.*y.^3));
        end
        
        function f = y4(x,y)
            f = (6561/40000).*(3750.*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+( ...
  -3).*y).^2).*((-5)+(-108).*y+1134.*y.^2+(-2916).*y.^3+2187.*y.^4)+ ...
  1875.*exp(1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*((-20) ...
  +144.*y+972.*y.^2+(-5832).*y.^3+6561.*y.^4)+48.*exp(1).^((-1/49).* ...
  (1+9.*x).^2+(-1/10).*(1+9.*y).^2).*(46+(-504).*y+(-1944).*y.^2+ ...
  2916.*y.^3+6561.*y.^4)+(-32000).*exp(1).^((-1).*(4+(-9).*x).^2+( ...
  -1).*(7+(-9).*y).^2).*(9019+(-47880).*y+94284.*y.^2+(-81648).* ...
  y.^3+26244.*y.^4)); 
        end
        
        function f = G(x,y), f = (9.0/1960)*(784*exp(-(4 - 9*x).^2 - (7-9*y).^2).*(-11 + 9*x + 9*y) - ...
                                490*exp(-(1.0/4)*(7 - 9*x).^2 - (9/4.0)*(1 - 3*y).^2).*(-10 + 9*x + 9*y) - ...
                                735*exp(-2 + 9*x - (81*x.^2)/4.0 + 9*y - (81*y.^2)/4.0).*(-4 + 9*x + 9*y) - ...
                                6*exp(-(1.0/49)*(1 + 9*x).^2 - (1.0/10)*(1 + 9*y).^2).*(59 + 90*x + 441*y));
        
        end
        
        function f = L(x,y)
            f =  (1/20).*((1215/4).*exp(1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).* ...
  y.^2).*((-36).*x+81.*x.^2+(2+(-9).*y).^2)+(3645/2).*exp(1).^(( ...
  -1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2).*(6+(-14).*x+9.* ...
  x.^2+(-6).*y+9.*y.^2)+(-1296).*exp(1).^((-1).*(4+(-9).*x).^2+(-1) ...
  .*(7+(-9).*y).^2).*(64+(-72).*x+81.*x.^2+(-126).*y+81.*y.^2)+( ...
  243/12005).*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*( ...
  (-11954)+1800.*x+8100.*x.^2+43218.*y+194481.*y.^2));
        end
        
        function f = B(x,y)
            f = (1/20).*((32805/8).*exp (1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+( ...
  -3).*y).^2).*(1825+(-10836).*x+22842.*x.^2+(-20412).*x.^3+6561.* ...
  x.^4)+(98415/16).*exp (1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).* ...
  y.^2).*((-20)+144.*x+972.*x.^2+(-5832).*x.^3+6561.*x.^4)+(-104976) ...
  .*exp (1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*(835+( ...
  -8352).*x+30132.*x.^2+(-46656).*x.^3+26244.*x.^4)+(393660/5764801) ...
  .*exp (1).^((-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*(6619+( ...
  -10440).*x+(-45684).*x.^2+11664.*x.^3+26244.*x.^4)+(98415/8).*exp ( ...
  1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2).*((-5)+(-108) ...
  .*y+1134.*y.^2+(-2916).*y.^3+2187.*y.^4)+(98415/16).*exp (1).^((-2) ...
  +9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*((-20)+144.*y+972.*y.^2+( ...
  -5832).*y.^3+6561.*y.^4)+(19683/125).*exp (1).^((-1/49).*(1+9.*x) ...
  .^2+(-1/10).*(1+9.*y).^2).*(46+(-504).*y+(-1944).*y.^2+2916.*y.^3+ ...
  6561.*y.^4)+(-104976).*exp (1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9) ...
  .*y).^2).*(9019+(-47880).*y+94284.*y.^2+(-81648).*y.^3+26244.* ...
  y.^4)+40.*((6561/32).*exp (1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+( ...
  -3).*y).^2).*(47+(-126).*x+81.*x.^2).*(7+(-54).*y+81.*y.^2)+( ...
  19683/64).*exp (1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*( ...
  2+(-36).*x+81.*x.^2).*(2+(-36).*y+81.*y.^2)+(19683/120050).*exp(1) ...
  .^((-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*((-47)+36.*x+162.* ...
  x.^2).*((-4)+18.*y+81.*y.^2)+(-26244/5).*exp (1).^((-1).*(4+(-9).* ...
  x).^2+(-1).*(7+(-9).*y).^2).*(31+(-144).*x+162.*x.^2).*(97+(-252) ...
  .*y+162.*y.^2)));
        end
        
        function f = p12(x,y)
            f = (729/39200).*((-2450).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+ ...
  (-3).*y).^2).*((-7)+9.*x).*(7+(-54).*y+81.*y.^2)+(-3675).*exp(1) ...
  .^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*((-2)+9.*x).*(2+( ...
  -36).*y+81.*y.^2)+(-48).*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10).*( ...
  1+9.*y).^2).*(1+9.*x).*((-4)+18.*y+81.*y.^2)+31360.*exp(1).^((-1) ...
  .*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*((-4)+9.*x).*(97+(-252).* ...
  y+162.*y.^2));  
        end
        
        function f = p21(x,y)
            f = (729/384160).*((-72030).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*( ...
  1+(-3).*y).^2).*(47+(-126).*x+81.*x.^2).*((-1)+3.*y)+307328.*exp( ...
  1).^((-1).*(4+(-9).*x).^2+(-1).*(7+(-9).*y).^2).*(31+(-144).*x+ ...
  162.*x.^2).*((-7)+9.*y)+(-36015).*exp(1).^((-2)+9.*x+(-81/4).* ...
  x.^2+9.*y+(-81/4).*y.^2).*(2+(-36).*x+81.*x.^2).*((-2)+9.*y)+(-48) ...
  .*exp(1).^((-1/49).*(1+9.*x).^2+(-1/10).*(1+9.*y).^2).*((-47)+36.* ...
  x+162.*x.^2).*(1+9.*y));
        end
        
        function f = p22(x,y)
            f = (6561/32).*exp(1).^((-1/4).*(7+(-9).*x).^2+(-9/4).*(1+(-3).*y).^2) ...
  .*(47+(-126).*x+81.*x.^2).*(7+(-54).*y+81.*y.^2)+(19683/64).*exp( ...
  1).^((-2)+9.*x+(-81/4).*x.^2+9.*y+(-81/4).*y.^2).*(2+(-36).*x+81.* ...
  x.^2).*(2+(-36).*y+81.*y.^2)+(19683/120050).*exp(1).^((-1/49).*(1+ ...
  9.*x).^2+(-1/10).*(1+9.*y).^2).*((-47)+36.*x+162.*x.^2).*((-4)+ ...
  18.*y+81.*y.^2)+(-26244/5).*exp(1).^((-1).*(4+(-9).*x).^2+(-1).*( ...
  7+(-9).*y).^2).*(31+(-144).*x+162.*x.^2).*(97+(-252).*y+162.*y.^2); 
        end
        
    end % Static methods
    
end % class