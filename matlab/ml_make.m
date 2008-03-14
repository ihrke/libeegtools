function ml_make;
%
% ml_make.m
% 
% Made by (Matthias Ihrke)
% Login   <mihrke@localhost>
% 
% Started on  Sun Jul 15 19:24:20 2007 Matthias Ihrke
%

unix('make')
mex -v ml_denoise.c denoising.o helper.o
