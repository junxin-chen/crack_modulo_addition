% This is part of the source code for a chosen-ciphertext attack which is given in
% 'Universal chosen-ciphertext attack for a family of image encryption
% schemes' (IEEE Transactions on Multimedia, vol **, no **, pp **-**, 2019).
% Preliminary results can also be found in: https://arxiv.org/abs/1903.11987


% This file is the code implementation of cat map permutation
% (encryption part)  mentioned in this paper


% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

function img_transed = arnold_trans(img,a,b,count)

% image size
[M,N] = size(img);
img_transed = zeros(M,N);

% permutation using cat map
for loop = 1:count
    for x=1:M
        for y=1:N
            x1=mod((x-1)+a*(y-1),M)+1;
            y1=mod(b*(x-1)+(a*b+1)*(y-1),N)+1;
            img_transed(x1,y1)=img(x,y);
        end
    end
    img = img_transed;
end

% ¸ñÊ½×ª»»
img_transed = uint8(img_transed);