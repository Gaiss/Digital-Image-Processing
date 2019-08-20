clear all
clc

lena = imread('lena.png');
lena = double(lena);
[m,n] = size(lena);
%% dft的图像信息
dft = DFT(lena);
lena_dft = dft*lena*dft/n;
lena_dft = shift(lena_dft);
lena_r = real(lena_dft);
lena_i = imag(lena_dft);
lena_p = angle(lena_dft)*180/pi;
figure(1)
subplot(2,2,1)
imshow(log(abs(lena_dft) + 1),[])
title('dft频谱图')
subplot(2,2,2)
imshow(log(abs(lena_r) + 1),[])
title('dft实部图')
subplot(2,2,3)
imshow(log(abs(lena_i) + 1),[])
title('dft虚部图')
subplot(2,2,4)
imshow(lena_p,[])
title('dft相位图')
%% dft的低频频谱还原
figure(2)
[m,n] = size(lena_dft);
fre1 = [8,16,32,64];
for i = 1:4
    lena_dft1 = zeros(m,n);
    t = (n/2-fre1(i)/2 + 1) : (n/2+fre1(i)/2);
    lena_dft1(t,t) = lena_dft(t,t);
    lena1 = inv(dft)*lena_dft1*inv(dft);
    subplot(2,2,i)
    imshow(abs(lena1),[])
    title(['低频',num2str(fre1(i)),'*',num2str(fre1(i)),'窗口'])
end
%% dft的高频频谱还原
figure(3)
fre2 = [256,128,64,32];
for i = 1:4
    lena_dft2 = lena_dft;
    t = (n/2-fre2(i)/2 + 1) : (n/2+fre2(i)/2);
    lena_dft2(t,t) = 0;
    lena2 = inv(dft)*lena_dft2*inv(dft);
    subplot(2,2,i)
    imshow(abs(lena2),[])
    title(['高频：挖去',num2str(fre1(i)),'*',num2str(fre1(i)),'窗口'])
end
%% dft的中频频谱还原
figure(4)
lena_dft3 = zeros(m,n);
for i =1:3
     t1 = (n/2-fre1(i)/2 + 1) : (n/2+fre1(i)/2);
     t2 = (n/2-fre2(i)/2 + 1) : (n/2+fre2(i)/2);
     lena_dft3(t2,t2) = lena_dft(t2,t2);
     lena_dft3(t1,t1) = 0;
     lena3 = inv(dft)*lena_dft3*inv(dft);
     subplot(2,2,i)
     imshow(abs(lena3),[])
     title(['中频范围为',num2str(fre1(i)),'到',num2str(fre2(i))])
end
%% dct图像信息
% dct = DCT(lena);
% a = zeros(1,n);
% a(1) = sqrt(1/n);
% for i = 2:n
%     a(i) = sqrt(2/n);
% end
% A = a'* a;
% lena_dct = A .* (dct*lena*dct);
% lena_dct = shift(lena_dct);
% lena_r = real(lena_dct);
% lena_i = imag(lena_dct);
% lena_p = angle(lena_dct)*180/pi;
% figure()
% subplot(2,2,1)
% imshow(log(lena_dct + 1),[])
% title('dct频谱图')
% subplot(2,2,2)
% imshow(log(abs(lena_r) + 1),[])
% title('dct实部图')
% subplot(2,2,3)
% imshow(log(abs(lena_i) + 1),[])
% title('dct虚部图')
% subplot(2,2,4)
% imshow(lena_p,[])
% title('dct相位图')
% %% dct的低频频谱还原
% figure()
% [m,n] = size(lena_dct);
% fre1 = [8,16,32,64];
% for i = 1:4
%     lena_dct1 = zeros(m,n);
%     t = (n/2-fre1(i)/2 + 1) : (n/2+fre1(i)/2);
%     lena_dct1(t,t) = lena_dct(t,t);
%     lena1 = inv(dct)*lena_dct1*inv(dct);
%     subplot(2,2,i)
%     imshow(abs(lena1),[])
%     title(['低频',num2str(fre1(i)),'*',num2str(fre1(i)),'窗口'])
% end
% %% dft的高频频谱还原
% figure()
% fre2 = [256,128,64,32];
% for i = 1:4
%     lena_dct2 = lena_dct;
%     t = (n/2-fre2(i)/2 + 1) : (n/2+fre2(i)/2);
%     lena_dct2(t,t) = 0;
%     lena2 = inv(dct)*lena_dct2*inv(dct);
%     subplot(2,2,i)
%     imshow(abs(lena2),[])
%     title(['高频：挖去',num2str(fre1(i)),'*',num2str(fre1(i)),'窗口'])
% end
% %% dct的中频频谱还原
% figure()
% lena_dct3 = zeros(m,n);
% for i =1:3
%      t1 = (n/2-fre1(i)/2 + 1) : (n/2+fre1(i)/2);
%      t2 = (n/2-fre2(i)/2 + 1) : (n/2+fre2(i)/2);
%      lena_dct3(t2,t2) = lena_dct(t2,t2);
%      lena_dct3(t1,t1) = 0;
%      lena3 = inv(dct)*lena_dct3*inv(dct);
%      subplot(2,2,i)
%      imshow(abs(lena3),[])
%      title(['中频范围为',num2str(fre1(i)),'到',num2str(fre2(i))])
% end