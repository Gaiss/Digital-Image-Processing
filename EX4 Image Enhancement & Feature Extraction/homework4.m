%% ֱ��ͼ���⻯
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
% title('ԭͼ��')
% subplot(2,3,2),imshow(H,[])
% title('�Աຯ�����⻯��ͼ')
% subplot(2,3,3),imshow(B,[])
% title('histeq���⻯��ͼ')
% subplot(2,3,4),imhist(A,[])
% title('ԭͼ��ֱ��ͼ')
% subplot(2,3,5),imhist(H,[])
% title('�Աຯ�����⻯��ֱ��ͼ')
% subplot(2,3,6),imhist(B,[])
% title('histeq���⻯��ֱ��ͼ')
%% ʹ��Matlab�ڲ�����edgeʵ�ֱ�Ե����㷨
ED1 = edge(A,'sobel');  % ��Sobel���ӽ��б�Ե���
ED2 = edge(A,'prewitt');  % ��Prewitt���ӽ��б�Ե���
ED3 = edge(A,'roberts');  % ��Roberts���ӽ��б�Ե���
ED4 = edge(A,'log');  % ��LOG���ӽ��б�Ե���
ED5 = edge(A,'canny');  % ��Canny���ӽ��б�Ե���
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
s_h = [-1 -2 -1;0 0 0;1 2 1];% Sobel�ݶ�ģ��
s_gx = conv2(A,s_h,'same'); % ����ͼ���Sobel��ֱ�ݶ� 
s_gy = conv2(A,s_h','same'); % ����ͼ���sobelˮƽ�ݶ�
s_g = abs(s_gx) + abs(s_gy);  %�õ�ͼ���sobel�ݶ�
sobel = A;
sobel(find(s_g > 250)) = 255; % ͼ���ֵ��
sobel(find(s_g <= 250)) = 0;
s_se = strel([0,0,0;0,1,0;0,0,1]);
sobel = imerode(sobel,s_se); % ��ʴ
figure()
subplot(1,2,1),imshow(ED1,[])
title('MATLAB�Դ�Sobel����')
subplot(1,2,2),imshow(sobel,[])
title('�Ա�Sobel����')
%% Prewitt
p_h = [-1 -1 -1;0 0 0;1 1 1];% Prewitt�ݶ�ģ��
p_gx = conv2(A,p_h,'same'); % ����ͼ���Prewitt��ֱ�ݶ� 
p_gy = conv2(A,p_h','same'); % ����ͼ���Prewittˮƽ�ݶ�
p_g = abs(p_gx) + abs(p_gy);  %�õ�ͼ���Prewitt�ݶ�
prewitt = A;
prewitt(find(p_g > 250)) = 255; % ͼ���ֵ��
prewitt(find(p_g <= 250)) = 0;
figure()
subplot(1,2,1),imshow(ED2,[])
title('MATLAB�Դ�Prewitt����')
subplot(1,2,2),imshow(prewitt,[])
title('�Ա�Prewitt����')
%% Roberts
r_h = [-1 0;0 1];% Roberts�ݶ�ģ��
r_gx = conv2(A,r_h,'same'); % ����ͼ���Roberts��ֱ�ݶ� 
r_gy = conv2(A,r_h','same'); % ����ͼ���Robertsˮƽ�ݶ�
r_g = abs(r_gx) + abs(r_gy);  %�õ�ͼ���Roberts�ݶ�
roberts = A;
roberts(find(r_g > 90)) = 255; % ͼ���ֵ��
roberts(find(r_g <= 90)) = 0;
figure()
subplot(1,2,1),imshow(ED3,[])
title('MATLAB�Դ�Roberts����')
subplot(1,2,2),imshow(roberts,[])
title('�Ա�Roberts����')
%% LOG
l_h = [0 0 1 0 0;0 1 2 1 0;1 2 -16 2 1;0 1 2 1 0;0 0 1 0 0];% Roberts�ݶ�ģ��
l_gx = conv2(A,l_h,'same'); % ����ͼ���Roberts��ֱ�ݶ� 
l_gy = conv2(A,l_h','same'); % ����ͼ���Robertsˮƽ�ݶ�
l_g = abs(l_gx) + abs(l_gy);  %�õ�ͼ���Roberts�ݶ�
LOG = A;
LOG(find(l_g > 250)) = 255; % ͼ���ֵ��
LOG(find(l_g <= 250)) = 0;
l_se = strel([0,0,1;0,1,0;0,1,0]);
LOG = imerode(LOG,l_se); % ��ʴ
figure()
subplot(1,2,1),imshow(ED4,[])
title('MATLAB�Դ�LOG����')
subplot(1,2,2),imshow(LOG,[])
title('�Ա�LOG����')
%% Canny
[m n] = size(A);
A = double(A);
% ��˹�˲�
c_h = fspecial('gaussian',[5 5]);
canny = imfilter(A,c_h,'replicate');
% sobel��Ե���
c_h = fspecial('sobel');
c_gx = imfilter(canny,c_h,'replicate');   % ����Ե
c_gy = imfilter(canny,c_h','replicate');    % ������Ե
canny = sqrt(c_gx.^2+c_gy.^2); 
c_gmax = max(canny(:));
canny = canny/c_gmax; % ��һ��
temp = canny;
counts = imhist(canny,64); 
high_threshold = find(cumsum(counts) > 0.7*m*n,1,'first') / 64; % �趨�Ǳ�Ե�����ȷ��������
low_threshold = 0.4 * high_threshold;  % �趨������Ϊ�����޳��Ա�������
% �Ǽ������� & ˫��ֵ
for q = 2 : m-1 
    for p = 2 : n-1 
        if((c_gx(q,p) <= 0 && c_gx(q,p) > -c_gy(q,p)) ...
                || (c_gy(q,p) >= 0 && c_gx(q,p) < -c_gy(q,p)))  % 0-45��
            d = abs(c_gy(q,p) / c_gx(q,p));  
            gradmag = canny(q,p); 
            gradmag1 = canny(q,p+1)*(1-d) + canny(q-1,p+1)*d;  % ��ֵ
            gradmag2 = canny(q,p-1)*(1-d) + canny(q+1,p-1)*d; 
        elseif((c_gx(q,p) > 0 && -c_gy(q,p) >= c_gx(q,p)) ...
                || (c_gx(q,p) < 0 && -c_gy(q,p) <= c_gx(q,p)))  % 45-90��
            d = abs(c_gx(q,p)/c_gy(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q-1,p)*(1-d) + canny(q-1,p+1)*d;  
            gradmag2 = canny(q+1,p)*(1-d) + canny(q+1,p-1)*d; 
        elseif((c_gx(q,p) <= 0 && c_gx(q,p) > c_gy(q,p)) ...
                || (c_gx(q,p) >= 0 && c_gx(q,p) < c_gy(q,p)))   % 90-135��
            d = abs(c_gx(q,p)/c_gy(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q-1,p)*(1-d) + canny(q-1,p-1)*d;  
            gradmag2 = canny(q+1,p)*(1-d) + canny(q+1,p+1)*d; 
        elseif((c_gy(q,p) < 0 && c_gx(q,p) <= c_gy(q,p)) ...
                || (c_gy(q,p) > 0 && c_gx(q,p) >= c_gy(q,p)))   % 135-180��
            d = abs(c_gy(q,p)/c_gx(q,p)); 
            gradmag = canny(q,p); 
            gradmag1 = canny(q,p-1)*(1-d) + canny(q-1,p-1)*d;  
            gradmag2 = canny(q,p+1)*(1-d) + canny(q+1,p+1)*d; 
        end 
        if(gradmag >= gradmag1 && gradmag >= gradmag2) % �ֲ�����ֵ�㣬����
            if(gradmag >= high_threshold)   % ������ֵ�Ͻ�ı������õ�ǿ�߽�
                temp(q,p) = 255; 
            elseif(gradmag >= low_threshold) % ������ֵ�½�ı������õ����߽�
                temp(q,p) = 125; 
            else 
                temp(q,p) = 0;  % ���ඪ��
            end 
        else 
            temp(q,p) = 0;  % �Ǿֲ�����ֵ�㣬����
        end 
    end  
end 
num = 0;      % ��ǰ��ջ����
flag = zeros(80000,2);       % ��ջ
temp_flag = zeros(80000,2);  % ��ʱ��ջ
for q = 2 : m-1      % ����ֵ�б𣬼�����ֵ����8������Χ���Ƿ���ڵ���ֵ
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
done = num;  % ��ɱ�־,����0��ʾ��ǰ���������
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
title('MATLAB�Դ�Canny����')
subplot(1,2,2),imshow(canny,[])
title('�Ա�Canny����')
%% hough�任��ȡԲ
% figure(3)
% imshow(A)
% c = []; r = [];
% for i = 1:2:9
%     [centers, radii] = imfindcircles(ED5,[10*i 10*i+79]); % Ѱ��ָ���뾶��Χ�ڵ�Բ
%     c = [c;centers]; % �ҳ���Բ������
%     r = [r;radii]; % �ҳ���Բ�İ뾶
% end
% viscircles(c,r,'EdgeColor','b'); % ���ҳ���Բ����ɫ�߻���