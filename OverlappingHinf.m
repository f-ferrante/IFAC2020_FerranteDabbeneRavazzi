%% Problem Definition
clear all;
A=[1 4 0; 1 2 2; 0 -2 3];
B=[1 0; 0 0; 0 1];
C=[0 1 0];
E=[0 1 0]';
np=max(size(A));
nu=min(size(B));
ny=min(size(C));
nw=min(size(E));
D=zeros(ny,nw);
S=[1 1 0; 0 1 1];
S=S(:);
le=max(size(S));
for i=1:le
    
    s=S(i);
    if s==1
         Ht=zeros(le,1);
         Ht(i)=1;
         H{i}=Ht;
    end
    
end
H=H(~cellfun('isempty',H));
for i=1:max(size(H))
T=H{i}; 
Lt{i}=reshape(T,[2,np]);   
end

k=max(size(Lt));
L=Lt{1};
for i=2:k
    L=[L,Lt{i}];
end
S=[1 1 0; 0 1 1];
%% Symbolic maniuplations to generate the structure of the matrix X
Lambda=sym('l',[k,k]);
Q=sym('q',[np,np]);
Left=L*kron(eye(k),Q);
Right=L*kron(Lambda,eye(np));
Eq=Left-Right;
Eq=Eq(:);
Asys=jacobian(Eq, Lambda(:));
Bsys=jacobian(Eq, Q(:))*Q(:);
Condition=Asys*pinv(Asys)*Bsys-Bsys;
ASysX=jacobian(Condition, Q(:));
Solution=null(ASysX);

Solution=Solution';
nc=min(size(Solution));
Coeff=sym('c',[nc,1]);
X=Coeff(1)*reshape(Solution(1,:),[np,np]);
for i=2:min(size(Solution)) 
 X=X+Coeff(i)*reshape(Solution(i,:),[np,np]);
end

XX=reshape(Solution(1,:),[np,np]);
for i=2:min(size(Solution)) 
 XX=XX+reshape(Solution(i,:),[np,np]);
end
%% Solution to the optimazion problem with line search on alpha
alpha_v=linspace(0.01,2,20);
for i=1:length(alpha_v)
Coeff=sdpvar(np,np,'full');
X=Coeff.*double(XX);
RR=sdpvar(nu,np,'full');
Y=RR.*S;
P=sdpvar(np,np,'symmetric');
gamma=sdpvar(1,1,'symmetric');
alpha=alpha_v(i);
MM=[A*X+B*Y, alpha*(A*X+B*Y)+P, zeros(np,ny), E;
    -X, -alpha*X, zeros(np,ny), zeros(np,nw);
    C*X, alpha*C*X, -gamma/2*eye(ny),D;
    zeros(nw, np),zeros(nw, np), zeros(nw, ny),-gamma/2*eye(nw); 
    ];
problem=[MM+MM'>=-200*eye(size(MM)), MM+MM'<=-1e-8*eye(size(MM)), P>=1e-5*eye(np)];

options=sdpsettings('solver','mosek','verbose',0);
solution=solvesdp(problem,gamma,options);
 if solution.problem==0
 gamma_v(i)=double(gamma);
 K=double(Y)*inv(double(X));
 alpha_o(i)=alpha;
 n_K(i)=norm(K);
 end
end

