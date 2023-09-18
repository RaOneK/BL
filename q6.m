e= 0.1

syms y(t)
[V] = odeToVectorField(diff(y, 2) ==-(1+e*3+cos(2*t))*y)

M = matlabFunction(V,'vars', {'t','Y'})
t=20
sol = ode45(M,[0 t],[1 0]);

fplot(@(x)deval(sol,x,1), [0, t])

%fplot(@(x) 