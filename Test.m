%简单的主程序，主要用于测试各类主加密解密程序的有效性
clear all
clc
close all

m=imread('lenna256.bmp');

encrypt=@(m)enc_pixel_level_modadd(m);
decrypt=@(m)dec_pixel_level_modadd(m);


c=encrypt(m);
d=decrypt(c);

dd=double(d)-double(m);

