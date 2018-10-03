function varargout = tuxiang(varargin)
% TUXIANG M-file for tuxiang.fig
%      TUXIANG, by itself, creates a new TUXIANG or raises the existing
%      singleton*.
%
%      H = TUXIANG returns the handle to a new TUXIANG or the handle to
%      the existing singleton*.
%
%      TUXIANG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TUXIANG.M with the given input arguments.
%
%      TUXIANG('Property','Value',...) creates a new TUXIANG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tuxiang_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tuxiang_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tuxiang_OpeningFcn, ...
                   'gui_OutputFcn',  @tuxiang_OutputFcn, ...
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


% --- Executes just before tuxiang is made visible.
function tuxiang_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tuxiang (see VARARGIN)

% Choose default command line output for tuxiang
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tuxiang wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tuxiang_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --------------------------------------------------------------------
function AC_Callback(hObject, eventdata, handles)
% hObject    handle to AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)



% --------------------------------------------------------------------
function AA_Callback(hObject, eventdata, handles)
% hObject    handle to AA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
[filename, pathname] = uigetfile( ...
       {'*.bmp;*.jpg;*.gif', 'All Files (*.bmp,*.jpg,*.gif)'; 
        '*.*','All Files (*.*)'}, ...
        '��������')
str=[pathname filename];
tu=imread(str);
subplot(1,2,1);
imshow(tu);
title('ԭͼ');








% --------------------------------------------------------------------
function ba_Callback(hObject, eventdata, handles)
% hObject    handle to ba (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%������
global tu
a=rgb2gray(tu);
I=imhist(a);   %����
I=histeq(a);   %������
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(I),title('ֱ��ͼ���⻯');


% --------------------------------------------------------------------
function bb_Callback(hObject, eventdata, handles)
% hObject    handle to bb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
global tu
a=rgb2gray(tu);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('�Ҷ�ͼ');


% --------------------------------------------------------------------
function bc_Callback(hObject, eventdata, handles)
% hObject    handle to bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
global tu
background=imopen(tu,strel('disk',80));
i=imsubtract(tu,background);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(i),title('ȥ����');


% --------------------------------------------------------------------
function bd_Callback(hObject, eventdata, handles)
% hObject    handle to bd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%��
global tu
a=graythresh(tu);
b=im2bw(tu,a);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(b),title('��ֵ��');

% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function ca_Callback(hObject, eventdata, handles)
% hObject    handle to ca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imadd(tu,50);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��������');






% --------------------------------------------------------------------
function cb_Callback(hObject, eventdata, handles)
% hObject    handle to cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imsubtract(tu,50);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('���ȼ���');




% --------------------------------------------------------------------
function cc_Callback(hObject, eventdata, handles)
% hObject    handle to cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=immultiply(tu,1.2);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('����(��һ������)');



% --------------------------------------------------------------------
function cd_Callback(hObject, eventdata, handles)
% hObject    handle to cd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imdivide(tu,1.5);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('���ʱ任(��һ������)');




% --------------------------------------------------------------------
function cfa_Callback(hObject, eventdata, handles)
% hObject    handle to cfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cfb_Callback(hObject, eventdata, handles)
% hObject    handle to cfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imrotate(tu,90,'bilinear');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('ͼ����ת90��');

% --------------------------------------------------------------------
function cf_Callback(hObject, eventdata, handles)
% hObject    handle to cf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cg_Callback(hObject, eventdata, handles)
% hObject    handle to cg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
figure;
imcontour(tu)


% --------------------------------------------------------------------
function ch_Callback(hObject, eventdata, handles)
% hObject    handle to ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function cfaa_Callback(hObject, eventdata, handles)
% hObject    handle to cfaa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imresize(tu,2,'nearest');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('����ڲ�ֵ');

% --------------------------------------------------------------------
function cfab_Callback(hObject, eventdata, handles)
% hObject    handle to cfab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imresize(tu,2,'bilinear');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('˫���Բ�ֵ');

% --------------------------------------------------------------------
function cfac_Callback(hObject, eventdata, handles)
% hObject    handle to cfac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imresize(tu,2,'bicubic');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('˫������ֵ');



% --------------------------------------------------------------------
function cfc_Callback(hObject, eventdata, handles)
% hObject    handle to cfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imcrop;
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('ͼ�����');

% --------------------------------------------------------------------
function cfd_Callback(hObject, eventdata, handles)
% hObject    handle to cfd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function daa_Callback(hObject, eventdata, handles)
% hObject    handle to daa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=edge(rgb2gray(tu),'sobel');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��Ե���sobel');

% --------------------------------------------------------------------
function dab_Callback(hObject, eventdata, handles)
% hObject    handle to dab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=edge(rgb2gray(tu),'prewitt');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��Ե���prewitt');


% --------------------------------------------------------------------
function dac_Callback(hObject, eventdata, handles)
% hObject    handle to dac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=edge(rgb2gray(tu),'roberts');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��Ե���roberts');


% --------------------------------------------------------------------
function dad_Callback(hObject, eventdata, handles)
% hObject    handle to dad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=edge(rgb2gray(tu),'log');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��Ե���log');


% --------------------------------------------------------------------
function dae_Callback(hObject, eventdata, handles)
% hObject    handle to dae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=edge(rgb2gray(tu),'zerocross');
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('��Ե���zerocross');


% --------------------------------------------------------------------
function da_Callback(hObject, eventdata, handles)
% hObject    handle to da (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function db_Callback(hObject, eventdata, handles)
% hObject    handle to db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dc_Callback(hObject, eventdata, handles)
% hObject    handle to dc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dd_Callback(hObject, eventdata, handles)
% hObject    handle to dd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function de_Callback(hObject, eventdata, handles)
% hObject    handle to de (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function d_Callback(hObject, eventdata, handles)
% hObject    handle to d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function dba_Callback(hObject, eventdata, handles)
% hObject    handle to dba (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imnoise(tu,'gaussian',0,0.1);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('gaussianz����');

% --------------------------------------------------------------------
function dbb_Callback(hObject, eventdata, handles)
% hObject    handle to dbb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imnoise(tu,'salt & pepper',0.05);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('salt & pepper����');

% --------------------------------------------------------------------
function dbc_Callback(hObject, eventdata, handles)
% hObject    handle to dbc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=imnoise(tu,'speckle',0.05);
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('speckle����');



% --------------------------------------------------------------------
function ea_Callback(hObject, eventdata, handles)
% hObject    handle to ea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=roipoly(tu)
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('�����ѡ��');

% --------------------------------------------------------------------
function Untitled_41_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function e_Callback(hObject, eventdata, handles)
% hObject    handle to e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function eb_Callback(hObject, eventdata, handles)
% hObject    handle to eb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
a=roicolor(rgb2gray(tu),155,255)
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(a),title('�ָ�ͼ��(�Ҷ�ֵ��155��255֮��)');

% --------------------------------------------------------------------
function Untitled_42_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Untitled_45_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_43_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --------------------------------------------------------------------
function fa_Callback(hObject, eventdata, handles)
% hObject    handle to fa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
ttu=rgb2gray(tu);
a=fftshift(fft2(ttu));
subplot(1,2,1),imshow(tu),title('ԭͼ');
subplot(1,2,2),imshow(log(abs(a)),[]),colormap(jet(64)),colorbar;
title('��ά����Ҷ�任ͼ');

% --------------------------------------------------------------------
function fb_Callback(hObject, eventdata, handles)
% hObject    handle to fb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
ttu=rgb2gray(tu);
a=dct2(ttu);
subplot(2,2,1),imshow(tu),title('ԭͼ');
subplot(2,2,3),imshow(log(abs(a)),[]),colormap(jet(64)),colorbar;
title('���ұ任ͼ');
a(abs(a)<10)=0;
k=idct2(a)/255;subplot(2,2,4),imshow(k),title('���ҷ��任�ָ�ͼ��');

% --------------------------------------------------------------------
function fc_Callback(hObject, eventdata, handles)
% hObject    handle to fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fd_Callback(hObject, eventdata, handles)
% hObject    handle to fd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function fca_Callback(hObject, eventdata, handles)
% hObject    handle to fca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tu
X=double(rgb2gray(tu))/255;
dwtmode('zpd');
lev=3,[C,S]=wavedec2(X,lev,'sym4');
a1=wrcoef2('a',C,S,'sym4',lev);
dwtmode('spd');
[C,S]=wavedec2(X,lev,'sym4');
a3=wrcoef2('a',C,S,'sym4',lev);
subplot(2,2,1),imshow(tu),title('ԭͼ');
subplot(2,2,2),imshow(a1),title('��sym4С���任�����ع����ͼ��');
subplot(2,2,3),imshow(a3),title('���⻬���С���任���ͼ��');


% --------------------------------------------------------------------
function fcb_Callback(hObject, eventdata, handles)
% hObject    handle to fcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function bf_Callback(hObject, eventdata, handles)
% hObject    handle to bf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


