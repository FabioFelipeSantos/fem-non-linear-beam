%%
clear
clc

syms L x v(x) a1 a2 a3 a4 v1 o1 v2 o2 xi
N1(x) = a1*x^3 + a2*x^2 + a3*x + a4;
dN1(x) = diff(N1,x,1);

cond = [N1(-1) == 1;dN1(-1) == 0;N1(1) == 0;dN1(1) == 0];

S = solve(cond, [a1,a2,a3,a4]);

N1 = S.a1*x^3 + S.a2*x^2 + S.a3*x + S.a4;

simplify(N1)

N2(x) = a1*x^3 + a2*x^2 + a3*x + a4;
dN2(x) = diff(N2,x,1);

cond = [N2(-1) == 0;dN2(-1) == 1;N2(1) == 0;dN2(1) == 0];

S = solve(cond, [a1,a2,a3,a4]);

N2 = S.a1*x^3 + S.a2*x^2 + S.a3*x + S.a4;

simplify(N2)

N3(x) = a1*x^3 + a2*x^2 + a3*x + a4;
dN3(x) = diff(N3,x,1);

cond = [N3(-1) == 0;dN3(-1) == 0;N3(1) == 1;dN3(1) == 0];

S = solve(cond, [a1,a2,a3,a4]);

N3 = S.a1*x^3 + S.a2*x^2 + S.a3*x + S.a4;

simplify(N3)

N4(x) = a1*x^3 + a2*x^2 + a3*x + a4;
dN4(x) = diff(N4,x,1);

cond = [N4(-1) == 0;dN4(-1) == 0;N4(1) == 0;dN4(1) == 1];

S = solve(cond, [a1,a2,a3,a4]);

N4 = S.a1*x^3 + S.a2*x^2 + S.a3*x + S.a4;

simplify(N4)

N5(x) = a1*x + a2;
dN5(x) = diff(N5,x,1);

cond = [N5(-1) == 1;N5(1) == 0];

S = solve(cond, [a1,a2]);

N5 = S.a1*x + S.a2;

simplify(N5)

N6(x) = a1*x + a2;
dN6(x) = diff(N6,x,1);

cond = [N6(-1) == 0;N6(1) == 1];

S = solve(cond, [a1,a2]);

N6 = S.a1*x + S.a2;

simplify(N6)

N = [N1, N2, N3, N4, N5, N6];
dN = diff(N,x,1)
d2N = diff(N, x, 2)

int(transpose(d2N(1:4)) * d2N(1:4),x,-1,1)
%%
clear
clc

syms N1(x) N2(x) a1 a2 L
N1(x) = a1*x + a2;

cond = [N1(-1) == 1;N1(1) == 0];

S = solve(cond, [a1,a2]);

N1(x) = S.a1*x + S.a2;

pretty(simplify(N1))

N2(x) = a1*x + a2;

cond = [N2(-1) == 0;N2(1) == 1];

S = solve(cond, [a1,a2]);

N2(x) = S.a1*x + S.a2;

pretty(simplify(N2))

%%
clear
clc

syms N1(x) N2(x) N3(x) N4(x) N5(x) N6(x) L

N1(x) = (1 - x) / 2;
N2(x) = (x + 1) / 2;
N3(x) = ((x - 1)^2 * (x + 2)) / 4;
N4(x) = ((x - 1)^2 * (x + 1)) / 4;
N5(x) = -((x + 1)^2 * (x - 2)) / 4;
N6(x) = ((x - 1) * (x + 1)^2) / 4;

disp('Derivadas Primeiras')
diff(N1)
diff(N2)
simplify(diff(N3))
simplify(diff(N4))
simplify(diff(N5))
simplify(diff(N6))

disp('Derivadas Segundas')
simplify(diff(N3,x,2))
simplify(diff(N4,x,2))
simplify(diff(N5,x,2))
simplify(diff(N6,x,2))

%%
clear
clc

syms x L
N1 = (x - 1)^2*(x + 2) / 4;
N2 = L*(x - 1)^2*(x+1) / 8;
N3 = -(x+1)^2*(x-2)/4;
N4 = L*(x-1)*(x+1)^2/8;

N = [N1, N2, N3, N4];
transpose(simplify(diff(N,x,1)))
transpose(simplify(diff(N,x,2)))
d2N = diff(N,x,2);

(8)*int(transpose(d2N)*d2N, x, -1, 1)

%%
clear
clc

syms x L a b c d xi

N1(x) = a*x^3 + b*x^2 + c*x + d;
dN1(x) = diff(N1,x,1);
cond = [N1(0)==1;dN1(0)==0;N1(L)==0;dN1(L)==0];
S = solve(cond,[a,b,c,d]);
N1(x) = S.a*x^3 + S.b*x^2 + S.c*x + S.d;
simplify(N1((L/2)*(1+xi)))
simplify(diff(N1((L/2)*(1+xi)), xi, 1))
simplify(diff(N1((L/2)*(1+xi)), xi, 2))

N2(x) = a*x^3 + b*x^2 + c*x + d;
dN2(x) = diff(N2,x,1);
cond = [N2(0)==0;dN2(0)==1;N2(L)==0;dN2(L)==0];
S = solve(cond,[a,b,c,d]);
N2(x) = S.a*x^3 + S.b*x^2 + S.c*x + S.d;
simplify(N2((L/2)*(1+xi)))
simplify(diff(N2((L/2)*(1+xi)), xi, 1))
simplify(diff(N2((L/2)*(1+xi)), xi, 2))

N3(x) = a*x^3 + b*x^2 + c*x + d;
dN3(x) = diff(N3,x,1);
cond = [N3(0)==0;dN3(0)==0;N3(L)==1;dN1(L)==0];
S = solve(cond,[a,b,c,d]);
N3(x) = S.a*x^3 + S.b*x^2 + S.c*x + S.d;
simplify(N3((L/2)*(1+xi)))
simplify(diff(N3((L/2)*(1+xi)), xi, 1))
simplify(diff(N3((L/2)*(1+xi)), xi, 2))

N4(x) = a*x^3 + b*x^2 + c*x + d;
dN4(x) = diff(N4,x,1);
cond = [N4(0)==0;dN4(0)==0;N4(L)==0;dN4(L)==1];
S = solve(cond,[a,b,c,d]);
N4(x) = S.a*x^3 + S.b*x^2 + S.c*x + S.d;
simplify(N4((L/2)*(1+xi)))
simplify(diff(N4((L/2)*(1+xi)), xi, 1))
simplify(diff(N4((L/2)*(1+xi)), xi, 2))

%%
clear
clc

syms x y L xi

A = [x xi 1;0 -1 1;L 1 1];
S = solve(det(A)==0, x)

N1(x) = 1 - x / L;
N2(x) = x/L;
N3(x) = 1 - (3 / L^2) * x^2 + (2 / L^3) * x^3;
N4(x) = x - (2 / L) * x^2 + (1 / L^2) * x^3;
N5(x) = (3 / L^2) * x^2 - (2 / L^3) * x^3;
N6(x) = -(1 / L) * x^2 + (1 / L^2) * x^3;

N(x) = [N1(x);N2(x);N3(x);N4(x);N5(x);N6(x)];
dN(x) = diff(N,x,1);


Nxi(xi) = simplify(N(S))
dNxi(xi) = simplify(diff(Nxi,xi,1))
d2Nxi(xi) = simplify(diff(dNxi,xi,1))

%%
clear
clc

syms x y L xi
A = [xi y 1;-1 0 1;1 L 1];
S(x) = solve(det(A)==0, y);
pretty(S)




















