%%================================================================================
% This function is to implement the data encryption algorith (TL-DEA) in
%        Y. Zhou, Z. Hua, C. M. Pun, and C. L. P. Chen, “Cascade chaotic system
%         with applications,” IEEE Trans. Cybern., vol. 45, no. 9, pp. 2001–2012,
%         Sep. 2015.
%All copyrights are reserved by Zhongyun Hua. E-mial:huazyum@gmail.com
%All following source code is free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
%If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper. 
%If you have any questions, please email me and I will try to response you ASAP.
%It worthwhiles to note that all following source code is written under MATLAB R2012a
%    and that files may call built-in functions from specific toolbox(es).
%%================================================================================
%%
function varargout = Zhou_TC_Cipher(P,para,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P: The input image, with 2-bit, 8-bit, 16-bit, 24-bit and 32-bit
% para: The operation model, either 'encryption' or 'decryption'
% K: the security key with length of 256 bits, if not provided, a random 
%    key will be generated and return
% varargout: when K is not given, return the result and the randomly
%            generated key; when K is given, return the result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% to get the key
if ~exist('K','var') && strcmp(para,'encryption')
    K = round(rand(1,256));
    OutNum = 2;
elseif ~exist('K','var')  && strcmp(para,'decryption')
    error('Can not dectrypted without a key');
else
    OutNum = 1;
end

%%extract the key
tran = @(K,low,high) sum(K(low:high).*2.^(-(1:(high-low+1))));


x0 = tran(K,1,52);
u0 = tran(K,53,104);
a0 = tran(K,105,156);
T = tran(K,157,208);

R = blkproc(K(209:256),[1,24],@(x) bi2de(x));


%% 
if max(P(:)) < 2
    F = 2;
elseif max(P(:)) < 2^8
    F = 8;
elseif max(P(:)) < 2^16
    F = 16;
elseif max(P(:)) < 2^24
    F = 24;
elseif max(P(:)) < 2^32
    F = 32;
else
    error('Can not support pixel format!');
end
%% To do the encryption algorithm
C = double(P);
[row, column] = size(C);
iter = 2;
switch para
    case 'encryption'
        for m = 1:iter
            x = mod(x0 + R(m)*T,1);
            u = 1.8 + mod(u0 + R(m)*T,0.2);
            a = 3.8 + mod(a0 + R(m)*T,0.2);
            X = ChaoticSeq(x,u,a,row*column);
            C = Substitution(para,C,X,F);
            C = Permutation(para,C,X);      
            C=uint8(C);
        end
    case 'decryption'
        for m = iter:-1:1
           x = mod(x0 + R(m)*T,1);
            u = 1.8 + mod(u0 + R(m)*T,0.2);
            a = 3.8 + mod(a0 + R(m)*T,0.2);
            X = ChaoticSeq(x,u,a,row*column);
            C = Permutation(para,C,X);
            C = Substitution(para,C,X,F);
            C=uint8(C);
            
        end
end

%% 
if OutNum == 1
    varargout{1} = C;
else
    varargout{1} = C;
    varargout{2} = K;
end
end % End of the main function
 

function X = ChaoticSeq(x,u,a,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (x,u,a): the initial state of the Tent-Logistic map 
% len: the length 
% X: the generated chaotic sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(1,len);
for m =1:len
    if x < 0.5
        x = u*x;
    else
        x = u*(1-x);
    end
    x = a*x*(1-x);
    X(m) = x;
end
end % The end of ChaoticSeq

function C = Substitution(para,P,X,F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para: The operation model, either 'encryption' or 'decryption'
% P: the input image
% X: the generated chaotic matrix with the same with of P
% F: the pixel format
% C: the operation output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%常见句柄
add_mod=@(a,b)uint8(mod(double(a)+double(b),256));%模和
inv_add_mod=@(c_sum,b)uint8(mod(double(c_sum)-double(b),256)); %已知模和结果，和一个变量，求另一个
sub_mod=@(a,b)uint8(mod(double(a)-double(b),256));%模减
inv_sub_mod_a=@(c_sub,b)uint8(mod(double(c_sub)+double(b),256));%已知模减结果和被减数b，求减数a
inv_sub_mod_b=@(c_sub,a)uint8(mod(double(c_sub)-double(a),256));%已知模减结果和减数a，求被减数b
%常见句柄

S = floor(X.*2^20);
S = mod(S,256);
[row,column] = size(P);
len = row * column;
P = reshape(P,[1,len]);
C = zeros(1,len);
switch para
    case 'encryption'
        cf1=P(len-1);
        cf0=P(len);
        for m = 1:len       
            if m == 1
%                 C(m) = mod(S(1) + P(1) + P(len) + P(len-1), 2^F);
                C(m) = add_mod(P(1),double(S(1))+double(P(len))+double(P(len-1)));
            elseif m == 2
%                 C(m) = mod(S(2) + P(2) + C(1) + P(len), 2^F);
                C(m) = add_mod(P(2),double(S(2))+double(C(m-1))+double(P(len)));
            else
%                  C(m) = mod(S(m) + P(m) + C(m-1) + C(m-2), 2^F);
                 C(m) = add_mod(P(m),double(S(m))+double(C(m-1))+double(C(m-2)));
            end
        end

                
    case 'decryption'
        for m = len:-1:1            
            if m == 1
%                 C(m) = mod(P(1) - S(1) - C(len) - C(len-1), 2^F);
                C(m) = inv_add_mod(P(1), double(S(1))+double(C(len-1))+double(C(len)));
            elseif m == 2
%                 C(m) = mod(P(2) - S(2) - P(1) - C(len), 2^F);
                C(m) = inv_add_mod(P(2), double(S(2))+double(P(1))+double(C(len)));
            else
%                 C(m) = mod(P(m) - S(m) - P(m-1) - P(m-2), 2^F);
                C(m) = inv_add_mod(P(m), double(S(m))+double(P(m-1))+double(P(m-2)));
            end
        end

end
C = reshape(C,[row,column]);
end % End of the Substitution function

function C = Permutation(para,P,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para: The operation model, either 'encryption' or 'decryption'
% P: the input image
% X: the generated chaotic matrix with the same with of P
% C: the operation output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row,column] = size(P);
X = reshape(X,[row,column]);

[X_sorted, I] = sort(X,2);
C = zeros(row,column);
switch para
    case 'encryption'
        for m = 1:row
            for n = 1:column
                x = mod(m - I(m,n)-1, row) + 1;
                y = I(x,:)==I(m,n);
                C(x,y) = P(m,n);
            end
        end
    case 'decryption'
        for m=1:row
            for n=1:column
                x = mod(m - (row-I(m,n))-1, row) + 1;
                y = I(x,:)==I(m,n);
                C(x,y) = P(m,n);                
            end
        end
end
end %End of the Permutation

