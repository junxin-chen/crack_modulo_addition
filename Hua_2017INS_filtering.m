%%================================================================================
%This functionto do image encryption using the reference in
%         [1]. Hua, Zhongyun, et al. "Design of image cipher using block-based scrambling and
%              image filtering.", Information Sciences 396 (2017): 97-113.
%All copyrights are reserved by Zhongyun Hua. E-mial:huazyum@gmail.com
%All following source code is free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
%If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper. 
%If you have any questions, please email me and I will try to response you ASAP.
%It worthwhile to note that all following source code is written under MATLAB R2010a
%    and that files may call built-in functions from specific toolbox(es).
%%================================================================================
%%

function varargout = Hua_2017INS_filtering(P,para,K)
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


K1 = KeyGen1(K);

%% 
% if max(P(:)>1)
%     F = 256;
% else
%     F = 2;
% end
F=256;
 %% To do the encryption/decryption
 C = double(P);
switch para
    case 'en'
        for m = 1:4
            C = BlkScrambling(C,K1(m),para);   
            C = rot90(C,m);
            C = ImageFiltering(C,K1(m),para);
        end
    case 'de'
        for m = 4:-1:1
            C = ImageFiltering(C,K1(m),para);
            C = rot90(C,4-m);
            C = BlkScrambling(C,K1(m),para);
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

function  kk = KeyGen1(K)
sK = zeros(4,32);
S = zeros(4,32);
kk = zeros(1,4);
for m = 1:4
    sK(m,:) = K((m-1)*32+1:m*32);
    S(m,:) = K((m-1)*32+129:m*32+128);
end
for m = 1:4
    t1 = xor(sK(m,:),S(m,:));
    kk(m) = bi2de(t1);
end
end 

function  kk = KeyGen2(K)
sK = zeros(2,16);
S = zeros(2,16);
kk = zeros(1,2);
for m = 1:2
    sK(m,:) = K((m-1)*16+193:m*16+192);
    S(m,:) = K((m-1)*16+225:m*16+224);
end
for m = 1:2
    t1 = xor(sK(m,:),S(m,:));
    kk(m) = bi2de(t1);
end
end 
%%
function R = PesudoRNG(seed,Len)
rand('seed',seed);
R = randi(2^32,[1,Len]);
end

%%
function C = BlkScrambling(P,kk,para)
[M,N] = size(P);
L = floor(min(M,N)^0.5);
R = PesudoRNG(kk,2*L*L);
A = R(1:L*L); B = R(L*L+1:2*L*L);
[A1,I] = sort(A); [B1,J] = sort(B);

O = zeros(L*L,L*L);
for j = 1:L*L
    for i = 1:L*L
        m = mod(i-J(j)-1,L*L) + 1;
        O(i,j) = I(m);
    end
end

C = P;

switch para
    case 'en'
        T = zeros(L,L*L*L);
        for m = 1:L
            T(:,(m-1)*L*L+1:m*L*L) = P((m-1)*L+1:m*L,1:L*L) ;
        end
        T = reshape(T,[L,L,L*L]);

        for i = 1:L*L
            for m = 1:L
                for n = 1:L
                    x = (m-1)*L+n;
                    C(x,O(i,x)) = T(m,n,i);
                end
            end
        end
    case 'de'
        T = zeros(L,L,L*L);
        for i = 1:L*L
            for m = 1:L
                for n = 1:L
                    x = (m-1)*L+n;
                    T(m,n,i) = C(x,O(i,x));
                end
            end
        end

        T = reshape(T,[L,L*L*L]);
        for m = 1:L
             C((m-1)*L+1:m*L,1:L*L) = T(:,(m-1)*L*L+1:m*L*L);
        end
end
end


%%
function C = ImageFiltering(P,kk,para)
[M,N] = size(P);
% if max(P(:)>1)
%     F = 256;
% else
%     F = 2;
% end

F = 256;

RR = PesudoRNG(kk,9+M*N);
R1 = mod(RR(1:M*N),F);R1 = reshape(R1,[M,N]);


R2 = mod(RR(M*N+1:M*N+9),F); A = reshape(R2,[3,3]);
A(3,3) = 1;

C = P;
switch para 
    case 'en'
        C = mod(C+R1,F);
         for m = 1:M
            for n = 1:N  
                T = zeros(3,3);
                if n == 1
                    if m == 1
                        T(:,1:2) = C((M-2):M,(N-1):N);
                        T(1:2,3) = C((M-1):M,n);
                        T(3,3) = C(m,n);
                    elseif m == 2
                        T(1:2,1:2) = C((M-1):M,(N-1):N);
                        T(3,1:2) = C(1,(N-1):N);
                        T(1,3) = C(M,1);
                        T(2:3,3) = C((m-1):m,1); 
                    elseif m == 3
                        T(1,1:2) = C(M,(N-1):N);
                        T(2:3,1:2) = C(1:2,(N-1):N);
                        T(:,3) = C((m-2):m,1);
                    else
                        T(:,1:2) = C((m-3):(m-1),(N-1):N);
                        T(:,3) = C((m-2):m,1);
                    end
                elseif n == 2
                    if m == 1
                        T(:,1) = C((M-2):M,N);
                        T(1:2,2:3) = C((M-1):M,1:2);
                        T(3,2:3) = C(m,(n-1):n);
                    elseif m == 2
                        T(1:2,1) = C((M-1):M,N);
                        T(3,1) = C(1,N);
                        T(1,2:3) = C(M,1:2);
                        T(2:3,2:3) = C((m-1):m,(n-1):n);
                    elseif m == 3
                        T(1,1) = C(M,N);
                        T(2:3,1) = C((m-2):(m-1),N);
                        T(:,2:3) = C((m-2):m,(n-1):n);
                    else
                        T(:,1) = C((m-3):(m-1),N);
                        T(:,2:3) = C((m-2):m,(n-1):n);
                    end
                elseif m == 1
                    T(1:2,:) = C((M-1):M,(n-2):n);
                    T(3,:) = C(m,(n-2):n);
                elseif m == 2
                    T(1,:) = C(M,(n-2):n);
                    T(2:3,:) = C((m-1):m,(n-2):n);
                else
                    T = C((m-2):m,(n-2):n);
                end
                TT = A.*T;
% 
                C(m,n) = mod(sum(TT(:)),256);      
           end
        end
    case 'de'
        for m = M:-1:1
            for n = N:-1:1
                 T = zeros(3,3);
                if n == 1
                    if m == 1
                        T(:,1:2) = C((M-2):M,(N-1):N);
                        T(1:2,3) = C((M-1):M,n);
                        T(3,3) = C(m,n);
                    elseif m == 2
                        T(1:2,1:2) = C((M-1):M,(N-1):N);
                        T(3,1:2) = C(1,(N-1):N);
                        T(1,3) = C(M,1);
                        T(2:3,3) = C((m-1):m,1); 
                    elseif m == 3
                        T(1,1:2) = C(M,(N-1):N);
                        T(2:3,1:2) = C(1:2,(N-1):N);
                        T(:,3) = C((m-2):m,1);
                    else
                        T(:,1:2) = C((m-3):(m-1),(N-1):N);
                        T(:,3) = C((m-2):m,1);
                    end
                elseif n == 2
                    if m == 1
                        T(:,1) = C((M-2):M,N);
                        T(1:2,2:3) = C((M-1):M,1:2);
                        T(3,2:3) = C(m,(n-1):n);
                    elseif m == 2
                        T(1:2,1) = C((M-1):M,N);
                        T(3,1) = C(1,N);
                        T(1,2:3) = C(M,1:2);
                        T(2:3,2:3) = C((m-1):m,(n-1):n);
                    elseif m == 3
                        T(1,1) = C(M,N);
                        T(2:3,1) = C((m-2):(m-1),N);
                        T(:,2:3) = C((m-2):m,(n-1):n);
                    else
                        T(:,1) = C((m-3):(m-1),N);
                        T(:,2:3) = C((m-2):m,(n-1):n);
                    end
                elseif m == 1
                    T(1:2,:) = C((M-1):M,(n-2):n);
                    T(3,:) = C(m,(n-2):n);
                elseif m == 2
                    T(1,:) = C(M,(n-2):n);
                    T(2:3,:) = C((m-1):m,(n-2):n);
                else
                    T = C((m-2):m,(n-2):n);
                end
                No = 0;
                for i = 1:3
                    for j = 1:3
                        if ~(i== 3 && j ==3)
                            No = No + T(i,j)*A(i,j);
                        end
                    end
                end
                C(m,n) = mod(C(m,n)-No,256);
            end
        end
        C = mod(C-R1,F);
end
end




