% This is part of the source code for a chosen-ciphertext attack which is given in
% 'Universal chosen-ciphertext attack for a family of image encryption
% schemes' (IEEE Transactions on Multimedia % Technology, vol **, no **, pp **-**, 2019).
% Preliminary results can also be found in: https://arxiv.org/abs/1903.11987


% This file is the core cracker algorithm
% mentioned in this paper


% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.


%% Start
clc
clear all


%%  handles used in the encryption

add_mod=@(a,b)uint8(mod(double(a)+double(b),256));% modulo addition
inv_add_mod=@(c_sum,b)uint8(mod(double(c_sum)-double(b),256)); % inverse of modulo addition
sub_mod=@(a,b)uint8(mod(double(a)-double(b),256));% modulo substraction
inv_sub_mod_a=@(c_sub,b)uint8(mod(double(c_sub)+double(b),256));% inverse of modulo substraction, solving the subtractor
inv_sub_mod_b=@(c_sub,a)uint8(mod(double(c_sub)-double(a),256));% inverse of modulo substraction, solving the minuend
mul_mod=@(lam,mm)uint8(mod(lam*double(mm),256));

%% function handles of the ciphers to be cryptanalyzed
% if you want to breack other ciphers, please change this part 

% % for the basic cipher
% encrypt=@(m)basic_enc_modadd(m);
% decrypt=@(m)basic_dec_modadd(m);

% % for Lan's cipher Signal Processing 147 (2018) 133每145
% encrypt=@(m)lan_enc(m);
% decrypt=@(m)lan_dec(m);

% % for Hua's cipher MIE-MA Signal Processing 144 (2018) 134每144
% load Hua_Medical_K
% encrypt=@(m)uint8(Hua_Hua_Medical_MA_Cipher(m,'en',K));
% decrypt=@(m)uint8(Hua_Hua_Medical_MA_Cipher(m,'de',K));

% for Zhou's cipher IEEE Trans. Cybern., vol. 45, no. 9, pp. 2001每2012, 2015. 
% load Zhou_TC_K
% encrypt=@(m)Zhou_TC_Cipher(m,'encryption',K);
% decrypt=@(m)Zhou_TC_Cipher(m,'decryption',K);

% for Hua's cipher Information Sciences 297 (2015): 80-94.
load Hua_K;
encrypt=@(m)Hua_ImageCipher(m,'encryption',K);
decrypt=@(m)Hua_ImageCipher(m,'decryption',K);



%% the ciphertext to be recovered
mm=imread('lenna32.bmp');  % the original file to be encrypted. 
% mm=imread('lenna32.bmp');  % the original file to be encrypted. 
% mm=imread('lenna256.bmp');  % the original file to be encrypted. 
% With respected to your computer's computational capacity, you can choose
% files with various sizes. 

cc=encrypt(mm); % this is the ciphertext to be recovered

%% the proposed chosen-ciphertext attack

[M,N]=size(cc); % get the ciphertext's size 
[M_plain,N_plain]=size(mm); %get the plaintext's size 
% The sizes of the ciphertext and plaintext maybe different, for example
% Hua's cipher in Z. Hua, S. Yi, and Y. Zhou, ※Medical image encryption using high-speed
% scrambling and pixel adaptive diffusion,§ Signal Processing, vol. 144,
% pp. 134 每 144, 2018.

c0=uint8(zeros(M,N));
atom_m0=uint8(decrypt(c0));
atom_delta_m=uint8(zeros(M,N,M_plain,N_plain));

% Steps 1-3: generate the atoms
for i=1:M
    i
    for j=1:N  
        current_c=uint8(zeros(M,N));
        current_c(i,j)=1;
        current_m=decrypt(current_c);
        atom_delta_m(i,j,:,:)=sub_mod(current_m,atom_m0);
    end
end

% Step 4: use these atoms to calculate the differential of the plaintext
delta_m=uint8(zeros(M_plain,N_plain));
lam=0;
for i=1:M
    i
    for j=1:N  
        lam=double(cc(i,j));
        current_atom(:,:)=atom_delta_m(i,j,:,:);
        delta_m=add_mod(delta_m,mul_mod(lam,current_atom));
    end
end
% Step 5: recover the plaintext
crack_m=add_mod(delta_m,atom_m0);

% check whether the plaintext has been accurately recovered
dd=double(mm)-double(crack_m);



