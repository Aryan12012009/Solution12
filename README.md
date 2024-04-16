# Solution12
% Conique mixti-linéaire externe

\[
\text{clc, clear all, close all}
\]

\[
\text{syms } a\ b\ c\ \text{real} % Lengths of the sides of triangle ABC
\]

% Conway's notations
\[
S_a=\frac{{b^2+c^2-a^2}}{2};\ S_b=\frac{{c^2+a^2-b^2}}{2};\ S_c=\frac{{a^2+b^2-c^2}}{2}
\]

\[
A=[1;\ 0;\ 0];\ B=[0;\ 1;\ 0];\ C=[0;\ 0;\ 1]; % Vertices of triangle ABC
\]
\[
BC=[1,\ 0,\ 0];\ CA=[0,\ 1,\ 0];\ AB=[0,\ 0,\ 1]; % Sides of triangle ABC
\]

%-----------------------------------------------------------------------

% Center of the circumcircle and square of its radius
\[
O=[a^2S_a;\ b^2S_b;\ c^2S_c]
\]
\[
R^2=\frac{{a^2b^2c^2}}{{(a+b+c)(a+b-c)(a-b+c)(b-a+c)}}
\]

% Center Ja of the A-excircle and square of its radius
\[
J_a = [-a;\ b;\ c]
\]
\[
r_a^2=\frac{{(a+b+c)(a+b-c)(a-b+c)}}{4(-a+b+c)}
\]

% Line DA passing through Ja orthogonal to the bisector (A Ja)
\[
A Ja=Wedge(A,J_a) % Line (A Ja): AJa=[0, c, -b]
\]
\[
DA=PgcdBary(DroiteOrthogonaleBary(J_a,AJa,a,b,c)) % DA=[2bc, c(a+b-c), b(a-b+c)]
\]

\[
N_{ab}=Wedge(AB,DA) % Nab=[c(-a-b+c); 2bc; 0]
\]
\[
N_{ac}=Wedge(CA,DA) % Nac=[b(-a+b-c); 0; 2bc]
\]

%-----------------------------------------------------------------------

\[
N_{ab}=[c(-a-b+c);\ 2bc;\ 0]
\]
\[
N_{ac}=[b(-a+b-c);\ 0;\ 2bc]
\]
% Similarly, by circular permutation:
\[
N_{bc}=[0;\ a(-b-c+a);\ 2ca]
\]
\[
N_{ba}=[2ca;\ c(-b+c-a);\ 0]
\]
\[
N_{ca}=[2ab;\ 0;\ b(-c-a+b)]
\]
\[
N_{cb}=[0;\ 2ab;\ a(-c+a-b)]
\]

% Conic passing through the first 5 points
\[
CoE=Conique5PointsBary(N_{ab},N_{ac},N_{bc},N_{ba},N_{ca})
\]
\[
Fact=(a^2-2ab+6ac+b^2-2bc+c^2)b(a+b-c)(b-a+c)
\]
\[
CoE=FactorT(Fact*CoE)
\]

% We find the coefficients:
\[
Cx^2=2bc(a+b-c)(a-b+c)(b-a+c)
\]
\[
Cy^2=2ac(a+b-c)(a-b+c)(b-a+c)
\]
\[
Cz^2=2ab(a+b-c)(a-b+c)(b-a+c)
\]
\[
Cxy=c(a-b+c)(b-a+c)(a^2+b^2+c^2 - 2bc - 2ca + 6ab)
\]
\[
Cyz=a(a+b-c)(a-b+c)(a^2+b^2+c^2 + 6bc - 2ca - 2ab)
\]
\[
Czx=b(a+b-c)(b-a+c)(a^2+b^2+c^2 - 2bc + 6ca - 2ab)
\]

% Equation of the conic
\[
G(x,y,z)=Cx^2x^2+Cy^2y^2+Cz^2z^2+Cxyxy+Cyzyz+Czxzx
\]
\[
N_{cb}=Factor(G(N_{cb}(1),N_{cb}(2),N_{cb}(3))) % Verification, Ncb is on the conic
\]

\[
Om=SimplifieBary(CentreConiqueBary(CoE(1),CoE(2),CoE(3),CoE(4)/2,CoE(5)/2,CoE(6)/2))
\]
\[
f(a,b,c)=-a^7+ 3(b+c)a^6 - (b^2+18bc+c^2)a^5 -(b+c)(5b^2-18bc+5c^2)a^4 + (b-c)^2(5b^2+46bc+5c^2)a^3 + (b+c)(b-c)^2(b^2-50bc+c^2)a^2 - (b-c)^2(b+3c)(3b+c)(b^2-6bc+c^2)a + (b+c)(b-c)^4(b^2+6bc+c^2)
\]
\[
Omega=[a(a+b-c)(a-b+c)f(a,b,c);\ b(a+b-c)(b-a+c)f(b,c,a);\ c(a-b+c)(b-a+c)f(c,a,b)]
\]
% Verification:
\[
NulOm=FactorT(Om-Omega) % We find NulOm[0; 0; 0]
\]
