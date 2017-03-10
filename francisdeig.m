function [H,Eig,iter] = francisdeig(A)
%francisdeig is a function that will compute the francis double shift
%eigenvalue algorithm.
%   Detailed explanation goes here

[H]=Hessenberg(A);  %Upper Hessenberg form of A.
m=min(size(A));
n=m; %lol hi
Eig=zeros(m,1);
tolerance=1e-9;
% for k=2:m
%     if (abs(H(k,k-1))<tolerance)
%         B=H(1:k-1,1:k-1);
%         C=H(k:m,k:m);  %Breaking H into seperate pieces if any of the
%         %values are sufficiently small.
%         [Eig,iter]=francisdeig(B);
%         E=Eig;
%         [Eig,iter]=francisdeig(C);
%         F=Eig;
%     end
% end

iter=0;
sig=@(u) sign(u) + (u==0);  %Changing the sign function to give sign(0)=1.

while (m>2)
    G=H(1:m,1:m);
    t=G(m-1,m-1)+G(m,m);  %Define trace of trailing 2x2 submatrix.
    d=G(m-1,m-1)*G(m,m)-G(m-1,m)*G(m,m-1);  %Define determinant of 2x2
    
    I=eye(m);
    e=I(1:m,1);  %e1 standard basis vector.
    x=G*(G*e)-t*(G*e)+d*e;  %p(A)*e1.
    y=x(1:3,1); %These are the nonzero entries of x.
    e1=e(1:3,1);  %  %Will use this to determine Q1~.
    v=sig(y(1))*norm(y)*e1+y;  %Vector defining the relfection.
    v=v/norm(v);  %Normalizing.
    Qnot=eye(3)-2*(v*v');  %Creating 3x3 reflector.
    
    G(1:3,1:m)=Qnot*G(1:3,1:m);
    G(1:m,1:3)=G(1:m,1:3)*Qnot;  %This completes the similarity 
                                 %transformation that creates the bulge.
   [G]=Hessenberg(G);                            
   H(1:m,1:m)=G;
%     for j=2:m-2 
%         %This loop returns H to upper Hessenberg form, taking advantage of 
%         %the fact that each reflector first acts on 3 rows, then on 3
%         %columns.
%         y=G(j:j+2,j);
%         v=sig(y(1))*norm(y)*e1+y;
%         v=v/norm(v);
%         Qref=eye(3)-2*(v*v');
%         G(j:j+2,1:m)=Qref*G(j:j+2,1:m);
%         G(j:m,j:j+2)=G(j:m,j:j+2)*Qref;  
%     end
%     
%     H(1:m,1:m)=G;
    
    Eig(m)=G(m,m);
    if (abs(G(m,m-1))<tolerance)
        G=G(1:m-1,1:m-1);
        m=m-1;  %Deflation
    end
    
    if m==2
       [Heig,itert]=ShiftQr(G);
       Eig(2)=Heig(2);
       Eig(1)=Heig(1);
    end
    
    iter=iter+1;
    if iter>100000 break; end
    
end

if m==2
       [Heig,itert]=ShiftQr(H(1:2,1:2));
       Eig(2)=Heig(2);
       Eig(1)=Heig(1);
end

end

