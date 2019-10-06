
function varargout = Hua_2019_Cosine(P,para,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main function to implement image cipher
% P:    the input image;
% para: operation type, 'en' or 'de';
% K:    the key, when para = 'en', it can be given or can not be given; 
%       when para = 'de', it must be given;
% C:    the encryption result
% K_d:  the deceyption key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% to get the key
if ~exist('K','var') && strcmp(para,'en')
    K = round(rand(1,256));
    OutNum = 2;
elseif ~exist('K','var')  && strcmp(para,'de')
    error('Can not dectrypted without a key');
else
    OutNum = 1;
end

%%extract the key to obtain the initial values for 4 round
IV = KeyGen(K);

%% 
if max(P(:)>1)
    F = 256;
else
    F = 2;
end
F=256;
 %% To do the encryption/decryption
[M,N] = size(P);
L = floor(min(M,N)^0.5);
C = double(P);
switch para
    case 'en'
        for m = 1:4
            S = ChaoticSeq(IV(:,m),4*L*L);
            C = HEScrambling(C,S,para);   
            C = rot90(C,m);
            S = ChaoticSeq(IV(:,m),M*N);
            C = ROSubstitution(C,S,para);
        end
    case 'de'
        for m = 4:-1:1
            S = ChaoticSeq(IV(:,m),M*N);
            C = ROSubstitution(C,S,para);
            C = rot90(C,4-m);
            S = ChaoticSeq(IV(:,m),4*L*L);
            C = HEScrambling(C,S,para);
        end
end

if F == 256
    C = uint8(C);
else
    C = logical(C);
end

%% 
if OutNum == 1
    varargout{1} = C;
else
    varargout{1} = C;
    varargout{2} = K;
end
end % end of the function ImageCipher


%%
function C = HEScrambling(P,S,para)
[M,N] = size(P);
L = floor(min(M,N)^0.5);
A = S(1:L*L); B = S(L*L+1:2*L*L);
Y = S(2*L*L+1:3*L*L); Z = S(3*L*L+1:4*L*L);
[A1,IA] = sort(A); [B1,IB] = sort(B);
[Y1,IY] = sort(Y); [Z1,IZ] = sort(Z);

O = zeros(L*L,L*L);
Q = O;
for i = 1:L*L
    for j = 1:L*L
        m = mod(i+IB(j)-1,L*L) + 1;
        O(i,j) = IA(m);
        n = mod(i+IZ(j)-1,L*L) + 1;
        Q(i,j) = IY(n);
    end
end

C = P;

switch para
    case 'en'
%         T = zeros(L,L,L*L);
%         for x = 1:L*L
%             for y = 1:L*L
%                 c = O(x,y);
%                 a = floor((double(Q(c,y))-1)/L)+1;
%                 b = mod((Q(c,y)-1),L)+1;
%                 
%                 T(a,b,c) = C(x,y);
%             end
%         end
%         T = reshape(T,[L,L*L*L]);
%         for m = 1:L
%              C((m-1)*L+1:m*L,1:L*L) = T(:,(m-1)*L*L+1:m*L*L);
%         end
        for j = 1:L*L
            for i = 1:L*L
                a = O(i,j); b = Q(a,j);
                m1 = floor((a-1)/L);n1 = mod((a-1),L);
                m2 = floor((b-1)/L)+1;n2 = mod((b-1),L)+1;
                
                x = m1*L+m2; y = n1*L+n2;
                C(x,y) = P(i,j);
            end
        end

    case 'de'
        for i = 1:L*L
            for j = 1:L*L
                a = O(i,j); b = Q(a,j);
                m1 = floor((a-1)/L);n1 = mod((a-1),L);
                m2 = floor((b-1)/L)+1;n2 = mod((b-1),L)+1;
                
                x = m1*L+m2; y = n1*L+n2;
                C(i,j) = P(x,y);
            end
        end
%         T = zeros(L,L*L*L);
%         for m = 1:L
%             T(:,(m-1)*L*L+1:m*L*L) = P((m-1)*L+1:m*L,1:L*L) ;
%         end
%         T = reshape(T,[L,L,L*L]);
% 
%          for x = 1:L*L
%             for y = 1:L*L
%                 c = O(x,y);
%                 a = floor((double(Q(c,y))-1)/L)+1;
%                 b = mod((Q(c,y)-1),L)+1;        
%                 C(x,y) = T(a,b,c);
%             end
%          end
end
end %% end of the function 

%% to do the random order substitution 
function C = ROSubstitution(P,S,para)
    [M,N] = size(P);
    S = reshape(S,[M,N]);
    [A1,I] = sort(S,1);
    
    P = double(P);
    [r,c] = size(P);

%     if (max(P(:))) > 1
%         F = 256;
%     else
%         F = 2;
%     end

    F=256;

    S = floor(S.*2^32);
    C = zeros(r,c);
    switch para
        case 'en'
            for x = 1:M
                for y = 1:N
                    if y == 1
                       if x == 1
                           C(I(x,y),y) = mod(P(I(x,y),y) + P(I(M,N),N) + S(I(x,y),y), F);
                       else
                           C(I(x,y),y) = mod(P(I(x,y),y) + C(I(x-1,N),N) + S(I(x,y),y), F);
                       end
                    else
                       C(I(x,y),y) = mod(P(I(x,y),y) + C(I(x,y-1),y-1) + S(I(x,y),y), F);
                    end
                end
            end
            
        case 'de'
            for x = M:-1:1
                for y = N:-1:1
                    if y == 1
                       if x == 1
                           C(I(x,y),y) = mod(P(I(x,y),y) - C(I(M,N),N) - S(I(x,y),y), F);
                       else
                           C(I(x,y),y) = mod(P(I(x,y),y) - P(I(x-1,N),N) - S(I(x,y),y), F);
                       end
                    else
                       C(I(x,y),y) = mod(P(I(x,y),y) - P(I(x,y-1),y-1) - S(I(x,y),y), F);
                    end
                end
            end
            
    end
    
end %end of the function 

%% To obtain the four piars of initial values for encryption 
function IV = KeyGen(K)
    tran = @(K,low,high) sum(K(low:high).*2.^(-(1:(high-low+1))));
    x0 = tran(K,1,32);
    r0 = tran(K,33,64);
    p = tran(K,65,96);
    g = bi2de(K(97:128));
    T = blkproc(K(129:256),[1,32],@(x) bi2de(x));
    
    IV = zeros(2,4);
    for m = 1:4
        x = mod(g*x0+T(m)*p, 1);
        r = mod(g*r0+T(m)*p, 1);
        IV(1,m) = x;
        IV(2,m) = r;
    end   
end %% end of the function 

%% To obtain the chaotic sequences
function S = ChaoticSeq(IV,Len)
    x = IV(1); r = IV(2);
    S = zeros(1,Len);
    for m = 1:1000
        x = cos(pi*(4*r*x*(1-x)+(1-r)*sin(pi*x)-0.5));
    end
    for m = 1:Len
        x = cos(pi*(4*r*x*(1-x)+(1-r)*sin(pi*x)-0.5));
        S(m) = x;
    end
end


