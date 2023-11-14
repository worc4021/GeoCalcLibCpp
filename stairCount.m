clear

i = 2;
n = 5;
Aeq = ones(1,i);
beq = n;
Aineq = -eye(i) + diag(ones(i-1,1),1);
bineq = -ones(i,1);


A = [Aeq;-Aeq;Aineq];
b = [beq;-beq;bineq];

[V,R] = vertexEnumeration(Aineq,bineq)

Aeq*V'-beq;

[A,b] = facetEnumeration(V,R)