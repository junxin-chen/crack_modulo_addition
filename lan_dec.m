% This is part of the source code for a chosen-ciphertext attack which is given in
% 'Universal chosen-ciphertext attack for a family of image encryption
% schemes' (IEEE Transactions on Multimedia % Technology, vol **, no **, pp **-**, 2019).
% Preliminary results can also be found in: https://arxiv.org/abs/1903.11987


% This file is Junxin Chen's implementation (decryption) of Lan's cipher proposed in 
% Lan' cipher Signal Processing 147 (2018) 133¨C145
% Please refer to the original paper for more detail


% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

function decipher=lan_dec(C)
%
tentmap=@(x)(x<0.5)*2*0.96*x+(x>=0.5)*2*0.96*(1-x);
logimap=@(x)4*0.96*x*(1-x);
sinemap=@(x)0.96*sin(pi*x);
chaos_map=@(x,fn,hn)mod((fn<0.5)*logimap(tentmap(x))+(fn>=0.5)*tentmap(tentmap(x))+(hn<0.5)*logimap(x)+(hn>=0.5)*sinemap(x),1);
sub_mod=@(a,b)uint8(mod(double(a)-double(b),256));%modular substraction

% Different keys are used in each of 4 encryption rounds
xx=[0.3330;0.556698;0.236547;0.366599];
ff=[0.25412;0.129874;0.3985234;0.785641];
hh=[0.932145;0.323652;0.785246;0.463289];
count=4;
[H,W]=size(C);
% decryption
decipher=C;
for kk=count:-1:1
    x0=xx(kk);
    f0=ff(kk);
    h0=hh(kk);    
    
for j=1:1:H+W+H*W
    chaos_var(j)=chaos_map(x0,f0,h0);
    x0=chaos_var(j);
    f0=sinemap(f0);
    h0=tentmap(h0);
end

% Produce T_r and T_c
% For T_r, the first step is to convert H chaotic variables to [0,H]
% For T_c, the first step is to convert W chaotic variables to [0,W]
% As this step is not described in detail in the original paper, so we use
% chaotic mapping to do it
% After obtaining T_r and T_c, then obtain W_r and W_c according to Eq.
% (13)
tr_1=chaos_var(1:H);
tr_2=sort(tr_1);
W_r=zeros(H,H);
for j=1:H
    T_r(j)=find(tr_2==tr_1(j));
    W_r(j,T_r(j))=1;
end
tc_1=chaos_var(H+1:H+W);
tc_2=sort(tc_1);
W_c=zeros(W,W);
for j=1:H
    T_c(j)=find(tc_2==tc_1(j));
    W_c(j,T_c(j))=1;
end

% Produce diffusion mask
mask_x=reshape(chaos_var(H+W+1:H+W+H*W),H,W);
mask_x=uint8(mask_x*256);

%Encryption
decipher=sub_mod(mask_x,decipher);
decipher=uint8(inv(W_r)*double(decipher)*inv(W_c));
end





