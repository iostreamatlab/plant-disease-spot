function varargout = drsvmgui(varargin)
% DRSVMGUI MATLAB code for drsvmgui.fig
%      DRSVMGUI, by itself, creates a new DRSVMGUI or raises the existing
%      singleton*.
%
%      H = DRSVMGUI returns the handle to a new DRSVMGUI or the handle to
%      the existing singleton*.
%
%      DRSVMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRSVMGUI.M with the given input arguments.
%
%      DRSVMGUI('Property','Value',...) creates a new DRSVMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drsvmgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drsvmgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drsvmgui
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drsvmgui_OpeningFcn, ...
                   'gui_OutputFcn',  @drsvmgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before drsvmgui is made visible.
function drsvmgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drsvmgui (see VARARGIN)

% Choose default command line output for drsvmgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drsvmgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drsvmgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im   %定义一个全局变量im
[filename,pathname]=uigetfile({'*.*';'*.bmp';'*.tif';'*.png'},'请选择多肉植物病害图片');  %选择图片路径
str=[pathname filename];  %合成路径+文件名
im=imread(str);   %读取图片
% if isrgb(im)
%     im=rgb2gray(im);
% end
axes(handles.axes1);  %使用第一个axes
imshow(im);  %显示图片
title('原始多肉植物病害图像');

tuxiang;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning off;
global im;
I=im;

% 色彩转换
cform = makecform('srgb2lab');
lab_he = applycform(I,cform);

% 构建颜色矩阵
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 3;
[cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);

segmented_images = cell(1,3);
% 构造色彩标签
rgb_label = repmat(pixel_labels,[1,1,3]);

for k = 1:nColors
    colors = I;
    colors(rgb_label ~= k) = 0;
    segmented_images{k} = colors;
end
% 显示分类内容
axes(handles.axes2);
imshow(segmented_images{2});title('多肉植物病斑分割');
% 特征解析

i = 3;

seg_img = segmented_images{i};

% 色彩转换
if ndims(seg_img) == 3
   img = rgb2gray(seg_img);
end

% 构造色阶
glcms = graycomatrix(img);

% 构造特征数据结构
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Mean = mean2(seg_img);
Standard_Deviation = std2(seg_img);
Entropy = entropy(seg_img);
RMS = mean2(rms(seg_img));
Variance = mean2(var(double(seg_img)));
a = sum(double(seg_img(:)));
Smoothness = 1-(1/(1+a));
Kurtosis = kurtosis(double(seg_img(:)));
Skewness = skewness(double(seg_img(:)));



m = size(seg_img,1);
n = size(seg_img,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = seg_img(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);

% 构造特征数组
%feat_disease = [Contrast,Correlation,Energy,Homogeneity, Mean(:,:,2), Standard_Deviation, Entropy, RMS, Variance, Smoothness, Kurtosis, Skewness, IDM];
%2014a版本



feat_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, RMS, Variance, Smoothness, Kurtosis, Skewness, IDM];  %2017b按本
% libsvm分类

% 读取训练矩阵
load Diseaseset.mat
svmStructDisease = svmtrain(diseasefeat,diseasetype);
% 得到训练结果
species_disease = svmclassify(svmStructDisease,feat_disease);
% 展示结果

set(handles.edit1,'string',species_disease);

if feat_disease(5)>90
    set(handles.edit2,'string','未感染');
elseif feat_disease(5)>50
    set(handles.edit2,'string','不严重（轻微）');
elseif feat_disease(5)>10
    set(handles.edit2,'string','较严重');
else 
    set(handles.edit2,'string','严重');
end

    









function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im;
% 坐标轴：
% 1。imshow图像是按矩阵显示（第一个下标为行，第二个下标为列）
% 2。plot本来是以左下角为坐标原点，右向为x轴正向，上方为y轴正向的。
% 3。如果在imshow出来的图形上用plot，则坐标轴的右方为x轴正向，下方为y轴正向，刚好分别对应和图像的列和行。

% 读图像与设置初始参数
t=im;
I = im2double(t);%fly_gray.jpg
%I = im2double(imread('mammo.bmp'));%fly_gray.jpg
axes(handles.axes3);
imshow(I);

hold on;
axis on;

% 图像的大小
sizeX = size(I,2);
sizeY = size(I,1);

% 图像的梯度
[gIx,gIy]=gradient(I);
I_grad = zeros(sizeX,sizeY);  % 梯度的模，先一次性计算出来，为了计算能量所用，而且调整了坐标，使得和plot函数的坐标一致
for xi=1:sizeX
    for yi=1:sizeY
         I_grad(xi,yi) = sqrt(gIx(yi,xi)^2 + gIy(yi,xi)^2);
    end
end

VCount = 200; % 轮廓上点的个数
V_init = zeros(VCount+2,2); % V是要找的轮廓，共有VCount个点，（V(i,1),V(i,2)）是轮廓上第i个点的坐标
% 初始化轮廓为一个圆（顺时针的顺序）
R =120;%150
for vi=1:VCount
    V_init(vi+1,1) = floor(130 + R * cos(vi * 2*3.14/VCount));
    V_init(vi+1,2) = floor(140 + R * sin(vi * 2*3.14/VCount));
end
% 超出边界的部分赋值为边界
V_init(find(V_init(:,1)>sizeX),1) = sizeX-4;
V_init(find(V_init(:,1)<1),1) = 2;
V_init(find(V_init(:,2)>sizeY),2) = sizeY-4;
V_init(find(V_init(:,2)<1),2) = 2;

V_init(1,:) = V_init(VCount+1,:);
V_init(VCount+2,:) = V_init(2,:);


hold on;
plot(V_init(:,1),V_init(:,2),'b');
plot(V_init(:,1),V_init(:,2),'b.');
plot(V_init(2,1),V_init(2,2),'b*');
plot(V_init(3,1),V_init(3,2),'b*');
%plot(V_init(:,1),V_init(:,2),'--rs','LineWidth',2, 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
% -----------------------------------------------------------------------------------------------------------
% 根据Snake算法寻找轮廓

repeatK = 500; % 迭代次数
minE = 10000; % 初始化的最小能量
V_min = V_init;  % 当前最小能量对应的轮廓

[I_gradSort,I_gradIndex] = sort(I_grad(:));
avgGrad = mean(I_gradSort((length(I_gradSort)-floor(length(I_gradSort)/50)):length(I_gradSort))); % 轮廓上大概的平均梯度，梯度最大的1/50的点的平均

for k = 1:repeatK
    
    vi=2; avgCon=0;
    while vi <= VCount
        viDist = norm(V_min(vi+1,:)-V_min(vi,:));
        % 如果两个点距离的很近，则合并
        if viDist <= 3
            for vii = (vi+1):(VCount+1)
                V_min(vii,:) = V_min(vii+1,:);
            end
            V_min(VCount+2,:) = [];
            VCount = VCount - 1;
        elseif viDist >= 7  % 如果两个点距离的很远，则中间插入一点
            for vii = (VCount+1):-1:(vi+1)
                V_min(vii+1,:) = V_min(vii,:);
            end
            V_min(vi+1,:) = floor((V_min(vi,:) + V_min(vi+2,:))/2);
            
            VCount = VCount + 1;

            V_min(1,:) = V_min(VCount+1,:);
            V_min(VCount+2,:) = V_min(2,:);
        else
            avgCon = avgCon + viDist; 
        end
        vi = vi + 1;
    end
    avgCon = avgCon / VCount;  % 所有点之间的平均距离，规范化所用 
    
%     avgCon = 0;
%     for vi = 2:(VCount+1)
%         avgCon = avgCon + norm(V_min(vi+1,:)-V_min(vi,:)); % 所有点之间的平均距离，规范化所用 
%     end
%     avgCon = avgCon / VCount;
    
    % 每次循环把所有点都向周围8个点移动一次，根据局部能量调整轮廓
    for vi = 2:(VCount+1)

        vi_01 = (V_min(vi,:) - V_min(vi-1,:)) / norm(V_min(vi,:) - V_min(vi-1,:));
        vi_12 = (V_min(vi+1,:) - V_min(vi,:)) / norm(V_min(vi+1,:) - V_min(vi,:));
        theta = acos( sum(vi_01 .* vi_12) );
        % 如果两个向量的夹角超过90度，则直接取中间的点
        if theta > 3.14/2 
            V_min(vi,:) = floor((V_min(vi-1,:) + V_min(vi+1,:))/2);
 
            % 防止重复的点
            if V_min(vi,1)==V_min(vi-1,1) && V_min(vi,2)==V_min(vi-1,2)
                V_min(vi,1) = V_min(vi,1) - 1;
            end

            if V_min(vi,1)==V_min(vi+1,1) && V_min(vi,2)==V_min(vi+1,2)
                V_min(vi,2) = V_min(vi,2) - 1;
            end
        
            V_min(1,:) = V_min(VCount+1,:);
            V_min(VCount+2,:) = V_min(2,:);
            continue;
        end
        
        % 该点的切线方向（顺时针方向，右边是轮廓内部）
        t_vi = (V_min(vi,:) - V_min(vi-1,:)) / norm(V_min(vi,:) - V_min(vi-1,:)) + (V_min(vi+1,:) - V_min(vi,:)) / norm(V_min(vi+1,:) - V_min(vi,:));

        % 法线方向（把切线逆时针旋转90度）
        n_vi = [t_vi(2),-t_vi(1)];
        n_vi = n_vi / norm(n_vi);

        Nei_vi = snake_Neighbor(sizeX,sizeY,V_min(vi,1),V_min(vi,2),V_min); % 得到轮廓上某点周围的点的坐标（8个）
        
        Nei_Egrad = zeros(1,size(Nei_vi,1));  % 第vi个点的领域的所有点的梯度（外部能量）
        Nei_Econ = zeros(1,size(Nei_vi,1)); % 平滑性约束（内部能量）
        Nei_Ebal = zeros(1,size(Nei_vi,1)); % 气球约束（内部能量）
        
        for ni=1:size(Nei_vi,1)
            
            Nei_Egrad(ni) = -I_grad(Nei_vi(ni,1),Nei_vi(ni,2));  % 外部能量（梯度）
            Nei_Econ(ni) = norm(Nei_vi(ni,:) - (V_min(vi-1,:) + V_min(vi+1,:))/2); % 平滑性约束
            
            % 邻域的点和原来的点在法线方向上的增长度量（前面的系数为正，要使得能量尽量小，则有往轮廓外扩张的趋势）
            Nei_Ebal(ni) = - sum( n_vi .* (V_min(vi,:)-Nei_vi(ni,:)) );
        end

        Nei_Egrad_norm = Nei_Egrad/avgGrad;
        Nei_Econ_norm = Nei_Econ/avgCon;
        
%         % 所有的能量都规范化成[0,1]，这样才可比较
%         if (max(Nei_Egrad)-min(Nei_Egrad))~=0
%             Nei_Egrad_norm = (Nei_Egrad-min(Nei_Egrad)) / (max(Nei_Egrad)-min(Nei_Egrad));
%         else
%             Nei_Egrad_norm = 0;
%         end
%         if (max(Nei_Econ)-min(Nei_Econ))~=0
%             Nei_Econ_norm = (Nei_Econ-min(Nei_Econ)) / (max(Nei_Econ)-min(Nei_Econ));
%         else
%             Nei_Econ_norm = 0;
%         end
        if (max(Nei_Ebal)-min(Nei_Ebal))~=0
            vin_grad = (I_grad(V_min(vi,1),V_min(vi,2)) + I_grad(V_min(vi-1,1),V_min(vi-1,2)) + I_grad(V_min(vi+1,1),V_min(vi+1,2))) / 3;
            Nei_Ebal_norm = (1 - vin_grad/avgGrad) * (Nei_Ebal-min(Nei_Ebal)) / (max(Nei_Ebal)-min(Nei_Ebal));
        else
            Nei_Ebal_norm = 0;
        end
        
%         Nei_E =  Nei_Econ_norm + Nei_Egrad_norm;
%         Nei_E =  Nei_Ebal_norm + Nei_Egrad_norm;
%         Nei_E =  Nei_Econ_norm + Nei_Ebal_norm;
        Nei_E = Nei_Econ_norm + Nei_Ebal_norm +15*Nei_Egrad_norm;
        
        % 取具有最小能量的点
        [Nei_ESort,Nei_EIndex] = sort(Nei_E);
        V_min(vi,:) = Nei_vi(Nei_EIndex(1),:);
        
        V_min(1,:) = V_min(VCount+1,:);
        V_min(VCount+2,:) = V_min(2,:);
    end
    
    if mod(k,10)==0
        plot(V_min(:,1),V_min(:,2),'r');
        plot(V_min(:,1),V_min(:,2),'r.');
    end
end

% -----------------------------------------------------------------------------------------------------------
% 画出结果轮廓
plot(V_min(:,1),V_min(:,2),'y');
plot(V_min(:,1),V_min(:,2),'y.');
title('处理过程');

% ------------------------------------------------------------------------------------------------------------


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
irisData = dlmread('drzw.data');


plotAxis = [min(irisData(:, 3)) max(irisData(:, 3)) min(irisData(:, 4)) max(irisData(:, 4))];
kernel = 'rbf';

c = [1 10 100 1000 10000];
sigma = 1 : 0.5: 3;
[C, Sigma] = meshgrid(c, sigma);
C = C(:);
Sigma = Sigma(:);
Err = zeros(size(C));
n = length(Err);
Xtest = irisData(1 : 150, 3 : 4);



X = irisData(1 : 100, 3 : 4);
Y = irisData(1 : 100, 5);
Y(find(Y==1)) = 1;
Y(find(Y==2)) = -1;
Xtrain = X;

for i = 1 : n
    [alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
    [YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));
  
end
[mErr, i] = min(Err);

[alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
[YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));

[Y12, Z] = yxcSVMclassifer(Xtrain, Xtest, Y, alphaStar, bStar, kernel, Sigma(i));
Y12(find(Y12==1)) = 1;
Y12(find(Y12==-1)) = 2;

X = irisData(51 : 150, 3 : 4);
Y = irisData(51 : 150, 5);
Y(find(Y==2)) = 1;
Y(find(Y==3)) = -1;
Xtrain = X;

for i = 1 : n
    [alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
    [YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));
  
end
[mErr, i] = min(Err);


[alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
[YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));

[Y23, Z] = yxcSVMclassifer(Xtrain, Xtest, Y, alphaStar, bStar, kernel, Sigma(i));
Y23(find(Y23==1)) = 2;
Y23(find(Y23==-1)) = 3;


X = irisData([1:50 101 : 150], 3 : 4);
Y = irisData([1:50 101 : 150], 5);
Y(find(Y==1)) = 1;
Y(find(Y==3)) = -1;
Xtrain = X;

for i = 1 : n
    [alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
    [YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));
 
end
[mErr, i] = min(Err);


[alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C(i), kernel, Sigma(i));
[YClassified, Z, Err(i)] = yxcSVMclassifer(Xtrain, Xtrain, Y, alphaStar, bStar, kernel, Sigma(i));

 
[Y13, Z] = yxcSVMclassifer(Xtrain, Xtest, Y, alphaStar, bStar, kernel, Sigma(i));
Y13(find(Y13==1)) = 1;
Y13(find(Y13==-1)) = 3;

Y123 = sort([Y12 Y23 Y13], 2);
Y123 = Y123(:, 2);


iY1 = find(Y123 == 1);
iY2 = find(Y123 == 2);
iY3 = find(Y123 == 3);
x = irisData(:, 3);
y = irisData(:, 4);
axes(handles.axes4);
plot(x(iY1), y(iY1), 'g+');
axis(plotAxis);
hold on;
plot(x(iY2), y(iY2), 'r*');
hold on;
plot(x(iY3), y(iY3), 'c.');
hold on;

% 错误点寻找
Y = irisData(:, 5);
iYwrong = find(Y123 ~= Y);
cwl=length(iYwrong)./length(Y);
zql=(1-cwl)*100;




plot(x(iYwrong), y(iYwrong), 'ko');
title('样本识别分布（黑色点为识别错误）')
hold off;
set(handles.edit3,'string',zql);



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
