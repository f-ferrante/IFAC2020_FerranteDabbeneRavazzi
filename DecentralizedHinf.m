%This code solves the Network Decentralized Control Problem in Section 5.1
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3
%% Problem Definition
clear all;
beta1=0;
beta2=0;
beta3=12;
beta4=0;
beta5=22;

alpha1=15;
alpha2=20;
alpha3=16;
alpha4=16.7;
alpha5=14;


A1=[-alpha1 beta1 0;
    alpha1 -beta1 0;
    0       1     0;
];

A2=[-alpha2 beta2 0;
    alpha2 -beta2 0;
    0       1     0;
];

A3=[-alpha3 beta3 0;
    alpha3 -beta3 0;
    0       1     0;
];

A4=[-alpha4 beta4 0;
    alpha4 -beta4 0;
    0       1     0;
];

A5=[-alpha5 beta5 0;
    alpha5 -beta5 0;
    0       1     0;
];


Bu=[1 0 0]';

Bd=[0 1 0]';


B=[Bu, -Bd, zeros(3,1), zeros(3,1),zeros(3,1),zeros(3,1);
   zeros(3,1), Bu, -Bd, zeros(3,1),zeros(3,1),-Bd;
   zeros(3,1),zeros(3,1), Bd, -Bu, zeros(3,1),zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),Bd, Bu, zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),-Bu, Bd;
    
];

A=blkdiag(A1, A2, A3, A4, A5);
np=max(size(A));
nu=min(size(B));
C=eye(np);
E=eye(np,np);
D=zeros(np,np);
%% Definition of the set\matchcal{S}(B') and contruction of a basis for \matchcal{S}(B')
Bup=[1 1 1]';
Bdp=[1 1 1]';
S=[Bup, -Bdp, zeros(3,1), zeros(3,1),zeros(3,1),zeros(3,1);
   zeros(3,1), Bup, -Bdp, zeros(3,1),zeros(3,1),-Bdp;
   zeros(3,1),zeros(3,1), Bdp, -Bup, zeros(3,1),zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),Bdp, Bup, zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),-Bup, Bdp;
    
];
S=abs(S');
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
Lt{i}=reshape(T,[6,15]);   
end

k=max(size(Lt));
L=Lt{1};
for i=2:k
    L=[L,Lt{i}];
end
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
S=[Bup, -Bdp, zeros(3,1), zeros(3,1),zeros(3,1),zeros(3,1);
   zeros(3,1), Bup, -Bdp, zeros(3,1),zeros(3,1),-Bdp;
   zeros(3,1),zeros(3,1), Bdp, -Bup, zeros(3,1),zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),Bdp, Bup, zeros(3,1);
   zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),-Bup, Bdp];
S=abs(S');
%% Solution to the optimazion problem with line search on alpha
nw=min(size(E));
ny=min(size(C));
alpha_v=linspace(0.01,1,50);
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
