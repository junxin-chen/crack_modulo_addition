% This is part of the source code for a chosen-ciphertext attack which is given in
% 'Universal chosen-ciphertext attack for a family of image encryption
% schemes' (IEEE Transactions on Multimedia, vol **, no **, pp **-**, 2019).
% Preliminary results can also be found in: https://arxiv.org/abs/1903.11987


% This file is the code implementation of the basic encryption model
% (cryption part)  mentioned in this paper


% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.


function p5 = basic_enc_modadd(p)

% handles used in the encryption
add_mod=@(a,b)uint8(mod(double(a)+double(b),256));% modulo addition
inv_add_mod=@(c_sum,b)uint8(mod(double(c_sum)-double(b),256)); % inverse of modulo addition
sub_mod=@(a,b)uint8(mod(double(a)-double(b),256));% modulo substraction
inv_sub_mod_a=@(c_sub,b)uint8(mod(double(c_sub)+double(b),256));% inverse of modulo substraction, solving the subtractor
inv_sub_mod_b=@(c_sub,a)uint8(mod(double(c_sub)-double(a),256));% inverse of modulo substraction, solving the minuend
mul_mod=@(lam,mm)uint8(mod(lam*double(mm),256));

% the key of the encryption, it is assumed unchanged in the attack
[M,N]=size(p);
miu=3.99976;
x0=0.92837471;
cat_a=33;
cat_b=44;
cat_counts=4;

count=10;

% prepositive iterations to promote the chaotic behavior of the employed chaos system  
for j=1:300
    x0=miu*x0*(1-x0);
end

% produce the subsitution masks
for i=1:M
    for j=1:N
        x0=miu*x0*(1-x0);
        key(i,j)=mod(floor(x0*10^14),256);        
    end
end
key=uint8(key);
key2=reshape(key.',1,M*N);

% begin encryption
for k=1:count
%     c0=123;
    p2=arnold_trans(p,cat_a,cat_b,cat_counts); % permutation using 3 rounds cat map
    p3=reshape(p2.',1,M*N);
    for j=1:M*N
       p4(j)=add_mod(p3(j),key2(j));  % substitution using modulo addition
    end
    p5=reshape(p4,M,N).';
    p=p5;
end

