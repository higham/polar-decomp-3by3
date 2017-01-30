function [Q, H] = polar_quaternion(A)
%POLAR_QUATERNION  Polar decompoition of 3-by-3 matrix.
%   For a real, 3-by-3 matrix A, [Q, H] = polar_quaternion(A)
%   computes the polar decompoition A = Q*H, where Q is orthogonal
%   and H is symmetric positive semidefinite.
%   This function uses a method based on the connection between orthogonal
%   matrices and quaternions.
    
%   Reference:
%   N. J. Higham and V. Noferini.
%   An algorithm to compute the polar decomposition of a $3 \times 3$
%   matrix. Numer. Algorithms, 73(2):349-369, 2016.

if ~all(size(A) == [3 3]), error('Matrix must be 3-by-3'), end
if ~isreal(A), error('Matrix must be real'), end

n = realsqrt(A(1,1)*A(1,1)+A(1,2)*A(1,2)+A(1,3)*A(1,3)...
             + A(2,1)*A(2,1)+A(3,1)*A(3,1)+A(2,2)*A(2,2)+A(3,2)*A(3,2)...
             + A(2,3)*A(2,3)+A(3,3)*A(3,3));
A = A/n(1); % Scale to norm 1.

subspa = 0;

m1 = (A(2,2)*A(3,3)-A(2,3)*A(3,2));b = m1(1)*m1(1);
m1 = (A(2,1)*A(3,3)-A(2,3)*A(3,1));b = b(1)+m1(1)*m1(1);
m1 = (A(2,1)*A(3,2)-A(2,2)*A(3,1));b = b(1)+m1(1)*m1(1);
m1 = (A(1,1)*A(3,2)-A(1,2)*A(3,1));b = b(1)+m1(1)*m1(1);
m1 = (A(1,1)*A(3,3)-A(1,3)*A(3,1));b = b(1)+m1(1)*m1(1);
m1 = (A(1,2)*A(3,3)-A(1,3)*A(3,2));b = b(1)+m1(1)*m1(1);
m1 = (A(1,2)*A(2,3)-A(1,3)*A(2,2));b = b(1)+m1(1)*m1(1);
m1 = (A(1,1)*A(2,3)-A(1,3)*A(2,1));b = b(1)+m1(1)*m1(1);
m1 = (A(1,1)*A(2,2)-A(1,2)*A(2,1));b = b(1)+m1(1)*m1(1);

b = b(1)*(-4)+1;
quick = 1;
if (b -1 + 1e-4) > 0
    quick = 0;
    % LU (full).
    r = 1;c = 1;AA = A;dd = 1;
    if abs(A(2,1))>abs(A(1,1))
        r = 2;
    end
    if abs(A(3,1))>abs(A(r(1),c(1)))
        r = 3;
    end
    if abs(A(1,2))>abs(A(r(1),c(1)))
        r = 1;c = 2;
    end
    if abs(A(2,2))>abs(A(r(1),c(1)))
        r = 2;c = 2;
    end
    if abs(A(3,2))>abs(A(r(1),c(1)))
        r = 3;c = 2;
    end
    if abs(A(1,3))>abs(A(r(1),c(1)))
        r = 1;c = 3;
    end
    if abs(A(2,3))>abs(A(r(1),c(1)))
        r = 2;c = 3;
    end
    if abs(A(3,3))>abs(A(r(1),c(1)))
        r = 3;c = 3;
    end

    if r(1)>1
        AA([1 r(1)],:) = AA([r(1) 1],:);dd = -1;
    end
    if c(1)>1
        AA(:,[1 c(1)]) = AA(:,[c(1) 1]);dd = -dd;
    end
    U = [0 0 0];U(1) = AA(1,1);
    m1 = AA(1,2)/AA(1,1);m2 = AA(1,3)/AA(1,1);
    AA = [AA(2,2)-AA(2,1)*m1(1) AA(2,3)-AA(2,1)*m2(1);...
          AA(3,2)-AA(3,1)*m1(1) AA(3,3)-AA(3,1)*m2(1)];

    r = 1;c = 1;
    if abs(AA(2,1))>abs(AA(1,1))
        r = 2;
    end
    if abs(AA(1,2))>abs(AA(r(1),c(1)))
        r = 1;c = 2;
    end
    if abs(AA(2,2))>abs(AA(r(1),c(1)))
        r = 2;c = 2;
    end
    if r(1) == 2
        dd = -dd;
    end

    if c(1)>1
        dd = -dd;
    end
    U(2) = AA(r(1),c(1));
    if U(2) == 0
        U(3) = 0;
    else
        U(3) = AA(3-r(1),3-c(1))-AA(r(1),3-c(1))*AA(3-r(1),c(1))/U(2);
    end
    
    d = dd(1);
    dd = dd(1)*U(1)*U(2)*U(3);
    
    if U(1)<0
        d = -d(1);
    end
    if U(2)<0
        d = -d(1);
    end
    if U(3)<0
        d = -d(1);
    end
    AU = abs(U(2));
    if AU>6.607e-8
        nit = 16.8+2*log10(AU);
        nit = ceil(15/nit);
    else
        subspa = 1;
    end
else
    %LU (partial).
    if abs(A(2,1))>abs(A(3,1))
        if abs(A(1,1))>abs(A(2,1))
            AA = A;
            dd = 1;
        else
            AA = A([2 1 3],:);
            dd = -1;
        end
    else
        if abs(A(1,1))>abs(A(3,1))
            AA = A;
            dd = 1;
        else
            AA = A([3 2 1],:);
            dd = -1;
        end
    end
    d = dd;
    U = [0 0 0];U(1) = AA(1,1);if U(1)<0, d(1) = -d(1);end
    m1 = AA(1,2)/AA(1,1);m2 = AA(1,3)/AA(1,1);
    AA = [AA(2,2)-AA(2,1)*m1(1) AA(2,3)-AA(2,1)*m2(1);...
          AA(3,2)-AA(3,1)*m1(1) AA(3,3)-AA(3,1)*m2(1)];
    
    if abs(AA(1,1))<abs(AA(2,1));
        U(2) = AA(2,1);
        U(3) = AA(1,2) - AA(1,1)*AA(2,2)/AA(2,1);
        dd = -dd;d = -d;
        if U(2)<0, d = -d(1);end
        if U(3)<0, d = -d(1);end
    elseif AA(1,1) == 0
        U(2) = 0;U(3) = 0;
    else
        U(2) = AA(1,1);
        U(3) = AA(2,2) - AA(2,1)*AA(1,2)/AA(1,1);
        if U(2)<0, d = -d(1);end
        if U(3)<0, d = -d(1);end
    end
    
    dd = dd(1)*U(1)*U(2)*U(3);
end


if d(1) == 0,d = 1;end

dd = 8*d(1)*dd(1);
t = A(1,1) + A(2,2) + A(3,3); %calling trace is slow

B = [t(1) A(2,3)-A(3,2) A(3,1)-A(1,3) A(1,2)-A(2,1)
    0 2*A(1,1)-t(1) A(1,2)+A(2,1) A(1,3)+A(3,1)
    0 0 2*A(2,2)-t(1) A(2,3)+A(3,2)
    0 0 0 2*A(3,3)-t(1)];

B = d(1)*B;
B(2,1) = B(1,2);B(3,1) = B(1,3);B(4,1) = B(1,4);B(3,2) = B(2,3);B(4,2) = B(2,4);B(4,3) = B(3,4);


%keyboard

% Find largest eigenvalue by analytic formula

if (b >= -0.3332) % && b <= 0.9999)
    Delta0  =  1 + 3*b(1);
    Delta1  =  -1 + (27/16)*dd(1)*dd(1) + 9* b(1);
    phi = Delta1(1)/Delta0(1);phi = phi/sqrt(Delta0(1));
    SS  =  (4/3)*(1+ cos(acos(phi(1))/3)*sqrt(Delta0(1)));
    S = sqrt(SS(1))/2;
    x = S(1)+0.5*realsqrt(max(0,-SS(1)+4+dd(1)/S(1)));
else
    x = sqrt(3); % When analytic approach is ill conditioned do Newton.
    xold = 3;
    while (xold-x)>1e-12;
        xold = x(1);
        px = x(1)*(x(1)*(x(1)*x(1) - 2)-dd(1))+b(1);
        dpx = x(1)*(4*x(1)*x(1)-4)-dd(1);
        x = x(1)-px(1)/dpx(1);
    end
end


if quick
    %LDL
    BB = -B;
    BB(1,1) = x(1)+BB(1,1);BB(2,2) = x(1)+BB(2,2);BB(3,3) = x(1)+BB(3,3);BB(4,4) = x(1)+BB(4,4);
    p = 1:4;
    L = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
    D = [0;0;0;0];
    
    %First step
    
    r = 4;if BB(4,4)<BB(3,3), r = 3; end; if BB(r(1),r(1))<BB(2,2), r = 2; end;
    if BB(r(1),r(1))>BB(1,1), p([1 r(1)]) = [r(1) 1];BB = BB(p,p); end
    D(1) = BB(1,1);
    L(2,1) = BB(2,1)/D(1);L(3,1) = BB(3,1)/D(1);L(4,1) = BB(4,1)/D(1);
    BB(2,2) = BB(2,2)-L(2,1)*BB(1,2);
    BB(3,2) = BB(3,2)-L(2,1)*BB(1,3);BB(2,3) = BB(3,2);
    BB(4,2) = BB(4,2)-L(2,1)*BB(1,4);BB(2,4) = BB(4,2);
    BB(3,3) = BB(3,3)-L(3,1)*BB(1,3);
    BB(4,3) = BB(4,3)-L(3,1)*BB(1,4);BB(3,4) = BB(4,3);
    BB(4,4) = BB(4,4)-L(4,1)*BB(1,4);
    
    
    %Second step
    
    r = 4;if BB(4,4)<BB(3,3), r = 3; end;
    if BB(r(1),r(1))>BB(2,2)
        p([2 r(1)]) = p([r(1) 2]);
        BB([2 r(1)],:) = BB([r(1) 2],:);
        BB(:,[2 r(1)]) = BB(:,[r(1) 2]);
        L([2 r(1)],:) = L([r(1) 2],:);
        L(:,[2 r(1)]) = L(:,[r(1) 2]);
    end
    
    D(2) = BB(2,2);
    L(3,2) = BB(3,2)/D(2);L(4,2) = BB(4,2)/D(2);
    BB(3,3) = BB(3,3)-L(3,2)*BB(2,3);
    BB(4,3) = BB(4,3)-L(3,2)*BB(2,4);BB(3,4) = BB(4,3);
    BB(4,4) = BB(4,4)-L(4,2)*BB(2,4);
    
    
    %Third step
    
    if BB(3,3)<BB(4,4)
        D(3) = BB(4,4);
        BB([3 4],:) = BB([4 3],:);
        BB(:,[3 4]) = BB(:,[4 3]);
        L([3 4],:) = L([4 3],:);
        L(:,[3 4]) = L(:,[4 3]);
        p([3 4]) = p([4 3]);
    else
        D(3) = BB(3,3);
    end
    L(4,3) = BB(4,3)/D(3);
    v = [L(2,1)*L(4,2)+L(3,1)*L(4,3)-L(2,1)*L(4,3)*L(3,2)-L(4,1);L(4,3)*L(3,2)-L(4,2);-L(4,3);1];
    nv = realsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)+v(4)*v(4));
    v = v/nv(1);
    v(p) = v;
else
    
    BB = -B;
    BB(1,1) = x(1)+BB(1,1);BB(2,2) = x(1)+BB(2,2);BB(3,3) = x(1)+BB(3,3);BB(4,4) = x(1)+BB(4,4);
    p = 1:4;
    L = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
    D = [0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
    
    %First step
    
    r = 4;if BB(4,4)<BB(3,3), r = 3; end; if BB(r(1),r(1))<BB(2,2), r = 2; end;
    if BB(r(1),r(1))>BB(1,1), p([1 r(1)]) = [r(1) 1];BB = BB(p,p); end
    D(1) = BB(1,1);
    L(2,1) = BB(2,1)/D(1);L(3,1) = BB(3,1)/D(1);L(4,1) = BB(4,1)/D(1);
    BB(2,2) = BB(2,2)-L(2,1)*BB(1,2);
    BB(3,2) = BB(3,2)-L(2,1)*BB(1,3);BB(2,3) = BB(3,2);
    BB(4,2) = BB(4,2)-L(2,1)*BB(1,4);BB(2,4) = BB(4,2);
    BB(3,3) = BB(3,3)-L(3,1)*BB(1,3);
    BB(4,3) = BB(4,3)-L(3,1)*BB(1,4);BB(3,4) = BB(4,3);
    BB(4,4) = BB(4,4)-L(4,1)*BB(1,4);
    
    %Second step
    
    r = 3;if BB(3,3)<BB(2,2), r = 2; end;
    if BB(r(1),r(1))>BB(1,1)
        p([2 r(1)]) = p([r(1) 2]);
        BB([2 r(1)],:) = BB([r(1) 2],:);
        BB(:,[2 r(1)]) = BB(:,[r(1) 2]);
        L([2 r(1)],:) = L([r(1) 2],:);
        L(:,[2 r(1)]) = L(:,[r(1) 2]);
        
    end
    
    D(2,2) = BB(2,2);
    L(3,2) = BB(3,2)/D(2,2);L(4,2) = BB(4,2)/D(2,2);
    D(3,3) = BB(3,3)-L(3,2)*BB(2,3);
    D(4,3) = BB(4,3)-L(3,2)*BB(2,4);D(3,4) = D(4,3);
    D(4,4) = BB(4,4)-L(4,2)*BB(2,4);
    
    DD = D(3,3)*D(4,4)-D(3,4)*D(3,4);
    if DD == 0, %treat specially
        if max(abs(D(3:4,3:4))) == 0, v = [L(2,1)*L(4,2)-L(4,1);-L(4,2);0;1];v = v/norm(v);
        else
            v = L'\[0;0;null(D(3:4,3:4))];v = v/norm(v);
        end        
    else
    ID = [D(4,4) -D(3,4); -D(3,4) D(3,3)];
    
    if subspa
        v = [L(2,1)*L(3,2)-L(3,1) L(2,1)*L(4,2)-L(4,1);-L(3,2) -L(4,2);1 0;0 1];
        IL = [1 0 0 0;-L(2,1) 1 0 0;v'];        
        [v ~] = qr(v,0);%->cost in flops if implemented by hand: 37 M+24 A+4 O
        v = IL*v;%it looks faster to multiply than to solve lin syst (even though should be same flops)
        v(1,:) = v(1,:)/D(1,1);
        v(2,:) = v(2,:)/D(2,2);
        v(3:4,:) = ID*v(3:4,:)/DD(1);
        v = v'*IL;v = v';
        v = IL*v;
        v(1,:) = v(1,:)/D(1,1);
        v(2,:) = v(2,:)/D(2,2);
        v(3:4,:) = ID*v(3:4,:)/DD(1);
        v = v'*IL;v = v';
                [v ~] = qr(v,0);
        H = v'*L;H = -H*D*H';%Cheaper
        if abs(H(1,2))<1e-15
            if H(1,1)>H(1,2)
                v = v(:,1);
            else
                v = v(:,2);
            end
        else
            r = (H(1,1)-H(2,2))/(2*H(1,2));
            v = v*[r+sign(H(1,2))*sqrt(1+r(1)*r(1));1];
            v = v/norm(v);
        end
        
    else
        v = [L(2,1)*L(4,2)+L(3,1)*L(4,3)-L(2,1)*L(4,3)*L(3,2)-L(4,1);
             L(4,3)*L(3,2)-L(4,2); -L(4,3) ;1];
        IL = [1 0 0 0; -L(2,1) 1 0 0; L(2,1)*L(3,2)-L(3,1) -L(3,2) 1 0; v'];
        nv = realsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)+v(4)*v(4));
        v = v/nv(1);
        
        for it = 1:nit
            v = IL*v;
            v(1) = v(1)/D(1,1);v(2) = v(2)/D(2,2);v(3:4) = ID*v(3:4)/DD(1);
            v = v'*IL;v = v';
            nv = realsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)+v(4)*v(4));
            v = v/nv(1);
        end
    end
    end
    v(p) = v;
end

% Polar factor (up to sign).

v22 = 2*v(2)*v(2);v33 = 2*v(3)*v(3);v44 = 2*v(4)*v(4);v23 = 2*v(2)*v(3);
v14 = 2*v(1)*v(4);v24 = 2*v(2)*v(4);v13 = 2*v(1)*v(3);v12 = 2*v(1)*v(2);
v34 = 2*v(3)*v(4);
Q = [1-v33(1)-v44(1) v23(1)+v14(1) v24(1)-v13(1);
    v23(1)-v14(1) 1-v22(1)-v44(1) v12(1)+v34(1); 
    v13(1)+v24(1) v34(1)-v12(1) 1-v22(1)-v33(1)]; % Slower in practice.

if d(1) == -1,Q = -Q;end
H = Q'*A;
H = n(1)*H;
%H = (H+H')/2; %->optional; adds some flop cost
end
