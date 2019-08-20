%% 直方图均衡化
clear all
A = imread('plate.jpg');
% B = histeq(A);
% H = zeros(numel(A),1);
% hp = 0;
% for i = 0:255
%     index = find(A == i);
%     hp = hp + length(index);
%     H(index) = hp/numel(A);
% end
% Range = max(max(A))-min(min(A));
% H = uint8(round(double(Range) * H+double(min(min(A)))+0.5));
% H = reshape(H,size(A));
% 
% figure(1)
% subplot(2,3,1),imshow(A,[])
% title('原图像')
% subplot(2,3,2),imshow(H,[])
% title('自编函数均衡化后图')
% subplot(2,3,3),imshow(B,[])
% title('histeq均衡化后图')
% subplot(2,3,4),imhist(A,[])
% title('原图像直方图')
% subplot(2,3,5),imhist(H,[])
% title('自编函数均衡化后直方图')
% subplot(2,3,6),imhist(B,[])
% title('histeq均衡化后直方图')
%% 使用Matlab内部函数edge实现边缘检测算法
ED1 = edge(A,'sobel');  % 用Sobel算子进行边缘检测
ED2 = edge(A,'prewitt');  % 用Prewitt算子进行边缘检测
ED3 = edge(A,'roberts');  % 用Roberts算子进行边缘检测
ED4 = edge(A,'log');  % 用LOG算子进行边缘检测
ED5 = edge(A,'canny');  % 用Canny算子进行边缘检测
% figure(2)
% subplot(2,3,1),imshow(A,[])
% title('Original')
% subplot(2,3,2),imshow(ED1,[])
% title('Sobel')
% subplot(2,3,3),imshow(ED2,[])
% title('Prewitt')
% subplot(2,3,4),imshow(ED3,[])
% title('Roberts')
% subplot(2,3,5),imshow(ED4,[])
% title('LOG')
% subplot(2,3,6),imshow(ED5,[])
% title('Canny')
%% Sobel
s_h = [-1 -2 -1;0 0 0;1 2 1];% Sobel梯度模板
s_gx = conv2(A,s_h,'same'); % 计算图像的Sobel垂直梯度 
s_gy = conv2(A,s_h','same'); % 计算图像的sobel水平梯度
s_g = abs(s_gx) + abs(s_gy);  %得到图像的sobel梯度
sobel = A;
sobel(find(s_g > 250)) = 255; % 图像二值化
sobel(find(s_g <= 250)) = 0;
s_se = strel([0,0,0;0,1,0;0,0,1]);
sobel = imerode(sobel,s_se); % 腐蚀
figure()
subplot(1,2,1),imshow(ED1,[])
title('MATLAB自带Sobel算子')
subplot(1,2,2),imshow(sobel,[])
title('自编Sobel算子')
%% Prewitt
p_h = [-1 -1 -1;0 0 0;1 1 1];% Prewitt梯度模板
p_gx = conv2(A,p_h,'same'); % 计算图像的Prewitt垂直梯度 
p_gy = conv2(A,p_h','same'); % 计算图像的Prewitt水平梯度
p_g = abs(p_gx) + abs(p_gy);  %得到图像的Prewitt梯度
prewitt = A;
prewitt(find(p_g > 250)) = 255; % 图像二值化
prewitt(find(p_g <= 250)) = 0;
figure()
subplot(1,2,1),imshow(ED2,[])
title('MATLAB自带Prewitt算子')
subplot(1,2,2),imshow(prewitt,[])
title('自编Prewitt算子')
%% Roberts
r_h = [-1 0;0 1];% Roberts梯度模板
r_gx = conv2(A,r_h,'same'); % 计算图像的Roberts垂直梯度 
r_gy = conv2(A,r_h','same'); % 计算图像的Roberts水平梯度
r_g = abs(r_gx) + abs(r_gy);  %得到图像的Roberts梯度
roberts = A;
roberts(find(r_g > 90)) = 255; % 图像二值化
roberts(find(r_g <= 90)) = 0;
figure()
subplot(1,2,1),imshow(ED3,[])
title('MATLAB自带Roberts算子')
subplot(1,2,2),imshow(roberts,[])
title('自编Roberts算子')
%% LOG
l_h = [0 0 1 0 0;0 1 2 1 0;1 2 -16 2 1;0 1 2 1 0;0 0 1 0 0];% Roberts梯度模板
l_gx = conv2(A,l_h,'same'); % 计算图像的Roberts垂直梯度 
l_gy = conv2(A,l_h','same'); % 计算图像的Roberts水平梯度
l_g = abs(l_gx) + abs(l_gy);  %得到图像的Roberts梯度
LOG = A;
LOG(find(l_g > 250)) = 255; % 图像二值化
LOG(find(l_g <= 250)) = 0;
l_se = strel([0,0,1;0,1,0;0,1,0]);
LOG = imerode(LOG,l_se); % 腐蚀
figure()
subplot(1,2,1),imshow(ED4,[])
title('MATLAB自带LOG算子')
subplot(1,2,2),imshow(LOG,[])
title('自编LOG算子')
%% Canny
[m n] = size(A);
A = double(A);
% 高斯滤波
c_h = fspecial('gaussian',[5 5]);
canny = imfilter(A,c_h,'replicate');
% sobel边缘检测
c_h = fspecial('sobel');
c_gx = imfilter(canny,c_h,'replicate');   % 求横边缘
c_gy = imfilter(canny,c_h','replicate');    % 求竖边缘
canny = sqrt(c_gx.^2+c_gy.^2); 
c_gmax = max(canny(:));
canny = canny/c_gmax; % 归一化
temp = canny;
counts = imhist(canny,64); 
high_threshold = find(cumsum(counts) > 0.7*m*n,1,'first') / 64; % 设定非边缘点比例确定高门限
low_threshold = 0.4 * high_threshold;  % 设定低门限为高门限乘以比例因子
% 非极大抑制 & 双阈值
for q = 2 : m-1 
    for p = 2 : n-1 
        if((c_gx(q,p) <= 0 && c_gx(q,p) > -c_gy(q,p)) ...
                || (c_gy(q,p) >= 0 && c_gx(q,p) < -c_gy(q,p)))  % 0-45度
            d = abs(c_gy(q,p) / c_gx(q,p));  
            gradmag = canny(q,p); 
            gradmag1 = canny(q,p+1)*(1-d) + canny(q-1,p+1)*d;  % 插值
            gradmag2 = canny(q,p-1)*(1-d) + canny(q+1,p-1)*d; 
        elseif((c_gx(q,p) > 0 && -c_gy(q,p) >= c_gx(q,p)) ...
                || (c_gx(q,p) < 0 && -c_gy(q,p) <= c_gx(q,p)))  % 45-90度
            d = abs(c_gx(q,p)/c_gy(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q-1,p)*(1-d) + canny(q-1,p+1)*d;  
            gradmag2 = canny(q+1,p)*(1-d) + canny(q+1,p-1)*d; 
        elseif((c_gx(q,p) <= 0 && c_gx(q,p) > c_gy(q,p)) ...
                || (c_gx(q,p) >= 0 && c_gx(q,p) < c_gy(q,p)))   % 90-135度
            d = abs(c_gx(q,p)/c_gy(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q-1,p)*(1-d) + canny(q-1,p-1)*d;  
            gradmag2 = canny(q+1,p)*(1-d) + canny(q+1,p+1)*d; 
        elseif((c_gy(q,p) < 0 && c_gx(q,p) <= c_gy(q,p)) ...
                || (c_gy(q,p) > 0 && c_gx(q,p) >= c_gy(q,p)))   % 135-180度
            d = abs(c_gy(q,p)/c_gx(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q,p-1)*(1-d) + canny(q-1,p-1)*d;  
            gradmag2 = canny(q,p+1)*(1-d) + canny(q+1,p+1)*d; 
        end 
        if(gradmag >= gradmag1 && gradmag >= gradmag2) % 局部极大值点，保留
            if(gradmag >= high_threshold)   % 大于阈值上界的保留，得到强边界
                temp(q,p) = 255; 
            elseif(gradmag >= low_threshold) % 大于阈值下界的保留，得到弱边界
                temp(q,p) = 125; 
            else 
                temp(q,p) = 0;  % 其余丢弃
            end 
        else 
            temp(q,p) = 0;  % 非局部极大值点，抑制
        end 
    end  
end 
num = 0;      % 当前堆栈个数
flag = zeros(80000,2);       % 堆栈
temp_flag = zeros(80000,2);  % 临时堆栈
for q = 2 : m-1      % 高阈值判别，检查高阈值邻域8个方向范围内是否存在低阈值
    for p = 2 : n-1
        if(temp(q,p) == 255)
            canny(q,p) = 255;
            if(temp(q-1,p-1) == 125)
                temp(q-1,p-1) = 255;
                canny(q-1,p-1) = 255;
                if((q-1 > 1) && (p-1 > 1))
                    num = num + 1;
                    flag(num,1) = q-1;
                    flag(num,2) = p-1;
                end
            end
            if(temp(q-1,p) == 125)
                temp(q-1,p) = 255;
                canny(q-1,p) = 255;
                if(q-1 > 1)
                    num = num + 1;
                    flag(num,1) = q-1;
                    flag(num,2) = p;
                end
            end
            if(temp(q-1,p+1) == 125)
                temp(q-1,p+1) = 255;
                canny(q-1,p+1) = 255;
                if((q-1 > 1) && (p+1 < n))
                    num = num + 1;
                    flag(num,1) = q-1;
                    flag(num,2) = p+1;
                end
            end
            if(temp(q,p-1) == 125)
                temp(q,p-1) = 255;
                canny(q,p-1) = 255;
                if(p-1 > 1)
                    num = num + 1;
                    flag(num,1) = q;
                    flag(num,2) = p-1;
                end
            end
            if(temp(q,p+1) == 125)
                temp(q,p+1) = 255;
                canny(q,p+1) = 255;
                if(p+1 < n)
                    num = num + 1;
                    flag(num,1) = q;
                    flag(num,2) = p+1;
                end
            end
            if(temp(q+1,p-1) == 125)
                temp(q+1,p-1) = 255;
                canny(q+1,p-1) = 255;
                if((q+1 < m) && (p-1 > 1))
                    num = num + 1;
                    flag(num,1) = q+1;
                    flag(num,2) = p-1;
                end
            end
            if(temp(q+1,p) == 125)
                temp(q+1,p) = 255;
                canny(q+1,p) = 255;
                if(q+1 < m)
                    num = num + 1;
                    flag(num,1) = q+1;
                    flag(num,2) = p;
                end
            end
            if(temp(q+1,p+1) == 125)
                temp(q+1,p+1) = 255;
                canny(q+1,p+1) = 255;
                if((q+1 < m) && (p+1 < n))
                    num = num + 1;
                    flag(num,1) = q+1;
                    flag(num,2) = p+1;
                end
            end
        end
    end
end
done = num;  % 完成标志,等于0表示当前连线已完成
while done ~= 0
    num = 0;
    for temp_num = 1:done
        q = flag(temp_num,1);
        p = flag(temp_num,2);
        if(temp(q-1,p-1) == 125)
            temp(q-1,p-1) = 255;
            canny(q-1,p-1) = 255;
            if((q-1 > 1) && (p-1 > 1))
                num = num + 1;
                temp_flag(num,1) = q-1;
                temp_flag(num,2) = p-1;
            end
        end
        if(temp(q-1,p) == 125)
            temp(q-1,p) = 255;
            canny(q-1,p) = 255;
            if(q-1 > 1)
                num = num + 1;
                temp_flag(num,1) = q-1;
                temp_flag(num,2) = p;
            end
        end
        if(temp(q-1,p+1) == 125)
            temp(q-1,p+1) = 255;
            canny(q-1,p+1) = 255;
            if((q-1 > 1) && (p+1 < n))
                num = num + 1;
                temp_flag(num,1) = q-1;
                temp_flag(num,2) = p+1;
            end
        end
        if(temp(q,p-1) == 125)
            temp(q,p-1) = 255;
            canny(q,p-1) = 255;
            if(p-1 > 1)
                num = num + 1;
                temp_flag(num,1) = q;
                temp_flag(num,2) = p-1;
            end
        end
        if(temp(q,p+1) == 125)
            temp(q,p+1) = 255;
            canny(q,p+1) = 255;
            if(p+1 < n)
                num = num + 1;
                temp_flag(num,1) = q;
                temp_flag(num,2) = p+1;
            end
        end
        if(temp(q+1,p-1) == 125)
            temp(q+1,p-1) = 255;
            canny(q+1,p-1) = 255;
            if((q+1 < m) && (p-1 > 1))
                num = num + 1;
                temp_flag(num,1) = q+1;
                temp_flag(num,2) = p-1;
            end
        end
        if(temp(q+1,p) == 125)
            temp(q+1,p) = 255;
            canny(q+1,p) = 255;
            if(q+1 < m)
                num = num + 1;
                temp_flag(num,1) = q+1;
                temp_flag(num,2) = p;
            end
        end
        if(temp(q+1,p+1) == 125)
            temp(q+1,p+1) = 255;
            canny(q+1,p+1) = 255;
            if((q+1 < m) && (p+1 < n))
                num = num + 1;
                temp_flag(num,1) = q+1;
                temp_flag(num,2) = p+1;
            end
        end
    end
    done = num;
    flag = temp_flag;
end
figure()
subplot(1,2,1),imshow(ED5,[])
title('MATLAB自带Canny算子')
subplot(1,2,2),imshow(canny,[])
title('自编Canny算子')
%% hough变换提取圆
% figure(3)
% imshow(A)
% c = []; r = [];
% for i = 1:2:9
%     [centers, radii] = imfindcircles(ED5,[10*i 10*i+79]); % 寻找指定半径范围内的圆
%     c = [c;centers]; % 找出的圆的中心
%     r = [r;radii]; % 找出的圆的半径
% end
% viscircles(c,r,'EdgeColor','b'); % 将找出的圆用蓝色线画出