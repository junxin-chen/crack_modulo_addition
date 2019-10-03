% This is the code implementation of Hua's cipher proposed in Signal Processing 144 (2018) 134–144

function varargout = Hua_Hua_Medical_MA_Cipher(P,para,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main function to simulte the medical image encryption using modulo arithmetic (MIE-MA)
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

%%extract the key
tran = @(K,low,high) sum(K(low:high).*2.^(-(1:(high-low+1))));
x0 = tran(K,1,52);
r = tran(K,53,104);
T1 = bi2de(K(105:128));
T2 = bi2de(K(129:152));
R1 = tran(K,153,204);
R2 = tran(K,205,256);

In1(1) = mod((x0 + R1)*T1,1);
In1(2) = mod((r + R1)*T1,4);

In2(1) = mod((x0 + R2)*T2,1);
In2(2) = mod((r + R2)*T2,4);

%% 
% if max(P(:)) < 2 
%     F = 1;
% elseif (max(P(:)) > 2^1) && (max(P(:)) < 2^8)
%     F = 8;
% elseif (max(P(:)) > 2^8) && (max(P(:)) < 2^16)
%     F = 16;
% elseif (max(P(:)) > 2^16) && (max(P(:)) < 2^24)
%     F = 24;
% elseif (max(P(:)) > 2^24) && (max(P(:)) < 2^32)
%     F = 32;
% else
%     error('Can not encryption for the pixel bit is bigger than 32');
% end
F=8;
 %% To do the encryption/decryption
 C = double(P);
switch para
    case 'en'
        
        C = ImageExpand(C,F);
        C = BlkScrambling(C,In1,para);
        C = ImageDiffusion(C,In1,para,F);
        
        C = BlkScrambling(C,In2,para);
        C = ImageDiffusion(C,In2,para,F);

    case 'de'
        C = ImageDiffusion(C,In2,para,F);
        C = BlkScrambling(C,In2,para);
        
        
        C = ImageDiffusion(C,In1,para,F);
        C = BlkScrambling(C,In1,para);
        
        C = C(2:(end-1),2:(end-1));
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
function R = PesudoRNG(In, Len)
R = zeros(1,Len);
R(1) = In(1);
r = In(2);
for m = 2:Len
    R(m) = mod(r*R(m-1)*(1-R(m-1))+(4-r)*sin(pi*R(m-1))/4, 1);
end
end
%%
function C = ImageExpand(P,F)
[r,c] = size(P);
rand('seed',sum(100*clock)*rand(1));
Ro = randi((2^F-1),[2,c]);
Co = randi((2^F-1),[r+2,2]);
C = zeros(r+2,c+2);
C(2:(1+r),2:(1+c)) = P;
C(1,2:(1+c)) = Ro(1,:); C(r+2,2:(1+c)) = Ro(2,:);
C(:,1) = Co(:,1); C(:,c+2) = Co(:,2);
end

%%
function C = BlkScrambling(P,In,para)
[rr,cc] = size(P);
L = PesudoRNG(In,rr+cc);
Ro = L(1:rr); Co = L(rr+1:rr+cc);
[Ro1,I] = sort(Ro); [Co1,J] = sort(Co);

S = zeros(rr,cc);
for m = 1:rr
    for n = 1:cc
        it = mod(n+I(m)-1,cc) + 1;
        S(m,n) = J(it);
    end
end

C = P;

switch para
    case 'en'
        for j = 1:cc
            for i = 1:rr
                r = i;
                c = S(i,j);
                m = mod(-S(1,j) + i - 1,rr) + 1;
                n = S(mod(-S(1,j)+i-1,rr)+1, j);
                C(m,n) = P(r,c);
            end
        end
    case 'de'
        for j = 1:cc
            for i = 1:rr
            r = i;
            c = S(i,j);
            m = mod(S(1,j) + i - 1,rr) + 1;
            n = S(mod(S(1,j)+i-1,rr)+1, j);
            C(m,n) = P(r,c);
            end
        end

end
end


%%
function C = ImageDiffusion(P,In,para,F)

%常见句柄
add_mod=@(a,b)uint8(mod(double(a)+double(b),256));%模和
inv_add_mod=@(c_sum,b)uint8(mod(double(c_sum)-double(b),256)); %已知模和结果，和一个变量，求另一个
sub_mod=@(a,b)uint8(mod(double(a)-double(b),256));%模减
inv_sub_mod_a=@(c_sub,b)uint8(mod(double(c_sub)+double(b),256));%已知模减结果和被减数b，求减数a
inv_sub_mod_b=@(c_sub,a)uint8(mod(double(c_sub)-double(a),256));%已知模减结果和减数a，求被减数b
%常见句柄

[r,c] = size(P);
FF = 2^F;
RR = PesudoRNG(In,r*c);
RR = mod(RR.*2^32, FF);
RR = floor(RR);
S = reshape(RR,[r,c]);
C = zeros(r,c);


switch para
    case 'en'
        for n = 1:c
            for m = 1:r
                if m == 1 && n == 1
%                     C(m,n) = mod(P(m,n)+P(r,c)+S(m,n), FF);
                    C(m,n) = add_mod(P(m,n),double(P(r,c))+double(S(m,n)));
                elseif m == 1
%                     C(m,n) = mod(P(m,n)+C(r,n-1)+S(m,n), FF);
                    C(m,n) = add_mod(P(m,n),double(C(r,n-1))+double(S(m,n)));
                else
%                     C(m,n) = mod(P(m,n)+C(m-1,n)+S(m,n), FF);
                    C(m,n) = add_mod(P(m,n),double(C(m-1,n))+double(S(m,n)));
                end
            end
        end
    case 'de'
        for n = c:-1:1
            for m = r:-1:1
                if m == 1 && n == 1
%                     C(m,n) = mod(P(m,n)-C(r,c)-S(m,n),FF);
                    C(m,n) = inv_add_mod(P(m,n),double(C(r,c))+double(S(m,n)));
                elseif m == 1
%                     C(m,n) = mod(P(m,n)-P(r,n-1)-S(m,n),FF);
                    C(m,n) = inv_add_mod(P(m,n),double(P(r,n-1))+double(S(m,n)));
                else
%                     C(m,n) = mod(P(m,n)-P(m-1,n)-S(m,n),FF);
                    C(m,n) = inv_add_mod(P(m,n),double(P(m-1,n))+double(S(m,n)));
                end
            end
        end
end
end




