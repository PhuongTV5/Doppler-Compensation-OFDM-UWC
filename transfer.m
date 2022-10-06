function varargout = transfer(varargin)
% Author Nguyen Quoc Khuong
% Ha noi University of Sience and Technology
% Radio electronic and Teclecommunication School
% mail:  khuong.nguyenquoc@mail.hust.vn

% This version insert 2 successive pilot in order to synchronization


%audioplayer
% TRANSFER MATLAB code for transfer.fig
%      TRANSFER, by itself, creates a new TRANSFER or raises the existing
%      singleton*.
%
%      H = TRANSFER returns the handle to a new TRANSFER or the handle to
%      the existing singleton*.
%
%      TRANSFER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRANSFER.M with the given input arguments.
%
%      TRANSFER('Property','Value',...) creates a new TRANSFER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before transfer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to transfer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help transfer

% Last Modified by GUIDE v2.5 05-Jul-2022 11:42:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @transfer_OpeningFcn, ...
                   'gui_OutputFcn',  @transfer_OutputFcn, ...
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


% --- Executes just before transfer is made visible.
function transfer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to transfer (see VARARGIN)

% Choose default command line output for transfer
handles.output = hObject;

% Continous_Tx handles structure
guidata(hObject, handles);

% UIWAIT makes transfer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = transfer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*txt'});
path=[pathname,filename];
if (filename==0)
    set(handles.load,'String','Cannot load');
    set(handles.transfer,'String','Cannot transfer')
    set(handles.transfer,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else
    set(handles.load,'String','Loaded');
    fid=fopen(path,'r');
    a=fileread(path); 
    fclose(fid);
    set(handles.disp_txt,'string',a);
    % Enable the Transfer button with its original name
    set(handles.transfer,'String','Transfer')
    set(handles.transfer,'Enable','on')
    savefile = 'pqfile.mat';
    save(savefile,'path');
end


function load2_Callback(hObject, eventdata, handles)
% hObject    handle to load2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in transfer.
[filename2,pathname2]=uigetfile({'*txt'});
path=[pathname2,filename2];
if (filename2==0)
    set(handles.load2,'String','Cannot load');
    set(handles.transfer,'String','Cannot transfer')
    set(handles.transfer,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else
    set(handles.load2,'String','Loaded');
    fid2=fopen(path,'r');
    a2=fileread(path); 
    fclose(fid2);
    set(handles.disp_txt2,'string',a2);
    % Enable the Transfer button with its original name
    set(handles.transfer,'String','Transfer')
    set(handles.transfer,'Enable','on')
    savefile = 'pqfile2.mat';
    save(savefile,'path2');
end



function transfer_Callback(hObject, eventdata, handles)
% hObject    handle to transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This  function use to transfer file through Sound  card 
% method can be   Bipolar,  BPSK,QPSK, OFDM,...

set(handles.stop,'Value',0); 


flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
 
 
f0=get(handles.edit_F0,'String');   
f0=str2double(f0);
 
  
df=get(handles.edit_df,'String');
df=str2double(df); 
 

tp=get(handles.time_period,'String');
tp=str2double(tp); 
  
bitrates=get(handles.edit_bitrates,'String');
bitrates=str2double(bitrates); 
 

Mod_Style=get(handles.popupmenu1,'String');       % Number of QASK modulation level
val = get(handles.popupmenu1,'Value');

M_ary=get(handles.ofdm_mod,'String');
M_ary=str2double(M_ary); 
        
payload_size=12;

y=get(handles.disp_txt,'String');
y2=get(handles.disp_txt2,'String');
y=y';
y2=y2';
y=y(:);
y2=y2(:);
y=y';
y2=y2';
L=length(y);
L2=length(y2);

  if L<L2
              for i=1:L2-L      
              y=[y 122];
              end %them bot ve do dai
  else
              for i=1:L-L2      
              y2=[y2 122];
              end  %them bot ve do dai
            end
y=y';
y2=y2';
L=length(y);
if mod(L,payload_size)~=0
    y=[y; ones(mod(L,payload_size),1)];
    y2=[y2; ones(mod(L,payload_size),1)];
end

 L=length(y); 
 z=[];
 z2=[];
 for k=0:L/payload_size-1
    yb=[];
    yb2=[];
    for i=1:payload_size
       yb=[yb de2bi(double(y(payload_size*k+i)),8,'left-msb')];
       yb2=[yb2 de2bi(double(y2(payload_size*k+i)),8,'left-msb')];
    end
    z=[z yb];
    z2=[z2 yb2];
 end
%  L2=length(y2);
%  z2=[];
%   for k2=0:L2/payload_size-1
%     yb2=[];
%     for i2=1:payload_size
%        yb2=[yb2 de2bi(double(y2(payload_size*k2+i2)),8,'left-msb')];
%     end
%     z2=[z2 yb2];
%  end

switch val
    case 1 %MIMO-OFDM  Luu day la phien ban co the dung truyen text
        %Chu y dieu che tham so 3000Hz-38000Hz voi 1024 hay 256 muc GI=128 NFFT=512 F=5 co the truyen
        %hoan hao  - tuy nhien khi thu len 3000-88000hz thi van truyen duoc
        %voi 256 muc cho truong hop NFFT=1024 GI=256- co le van de la phai dieu chinh tin hieu phat nam
        %duoc muc 1 de khong bi meo khi phat di (de GI=256 NFFT=1024 cung
        %duoc)
%      y=randi(256,3000,1)-1;
%      y2=randi(256,3000,1)-1;
       
        ofdm_f1=get(handles.ofdm_F1,'String');   
        ofdm_f1=str2double(ofdm_f1);

        ofdm_f2=get(handles.ofdm_F2,'String');   
        ofdm_f2=str2double(ofdm_f2);

        GI=get(handles.ofdm_GI,'String');
        GI=str2double(GI); 

        NFFT=get(handles.NFFT,'String');
        NFFT=str2double(NFFT); 
        
        D_t=get(handles.frame_size,'String');
        D_t=str2double(D_t);
        
        SubL=floor(2*ofdm_f1*NFFT/flm);
        SubL=SubL-mod(SubL,2);
        SubH=ceil(2*ofdm_f2*NFFT/flm);
        SubH=SubH+mod(SubH,2);
        set(handles.ofdm_SubL,'String',num2str(SubL)); 
        set(handles.ofdm_SubH,'String',num2str(SubH));  

        N_D= SubH- SubL; % so song mang con mang data
      
        PO1=[ones(1,N_D/2); zeros(1,N_D/2)];
        PO1=PO1(:);
        PO2=[zeros(1,N_D/2); ones(1,N_D/2)];
        PO2=PO2(:);
        
        pilot_mask=[zeros(1,SubL) PO1' zeros(1,NFFT-SubH)];
        
        pilot_mask2=[zeros(1,SubL) PO2' zeros(1,NFFT-SubH)];
        
        P_A = sqrt(M_ary);     % Amplitude of pilot symbol
        D_f = 2;
        K=NFFT*2+GI+1;

        switch NFFT
            case 64
        Pilot=[1 1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 0 -1 1 -1 1 -1 1 1 1 -1 1 -1 -1 1 -1 1 1 -1 -1   -1 1 -1 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
            case 256
               Pilot=  [ 1     1    -1     1     1     1    -1    -1    -1     1     1    -1     1     1    -1     1    -1    -1    -1     1    -1,...
     -1    -1    -1    -1     1     1    -1     1     1     1     1    -1    -1     1     1    -1     1    -1    -1     1     1,...
     1    -1     1     1     1     1     1     1    -1    -1     1    -1     1     1    -1    -1     1    -1     1    -1    -1,...
    -1    -1    -1     1     1    -1     1    -1    -1    -1     1     1     1    -1     1    -1    -1     1     1     1     1,...
     1    -1     1     1    -1     1     1    -1     1    -1    -1    -1     1    -1     1    -1    -1    -1     1    -1     1,...
     1     1     1     1     1     1    -1     1    -1    -1    -1     1     1     1    -1    -1     1    -1    -1    -1    -1,...
     1     1    -1    -1     1     1     1     1     1     1    -1     1    -1    -1    -1    -1    -1     1    -1     1    -1,...
    -1     1     1     1    -1    -1     1    -1    -1     1    -1     1    -1     1    -1     1    -1     1     1     1    -1,...
     1     1     1    -1     1    -1     1    -1     1    -1     1     1     1    -1    -1    -1     1    -1    -1    -1    -1,...
    -1     1    -1     1     1    -1     1     1    -1    -1     1    -1     1    -1    -1     1    -1    -1    -1    -1     1,...
     1     1    -1    -1     1    -1    -1    -1    -1    -1    -1    -1     1     1     1     1     1     1     1    -1    -1,...
    -1    -1    -1     1     1    -1     1    -1    -1    -1     1     1    -1    -1    -1     1     1    -1     1     1    -1,...
     1     1    -1     1];
       case 512
            load('ofdm_pilot_512.mat');
            Pilot=P;
        case 1024
            load('ofdm_pilot_1024.mat');
            Pilot=P;
        case 2048
            load('ofdm_pilot_2048.mat');
            Pilot=P;
        case 8192
            load('ofdm_pilot_8192.mat');
            Pilot=P;
        end
        size(pilot_mask2)
        Pilot2=sqrt(M_ary)*Pilot.*pilot_mask2;
        Pilot=sqrt(M_ary)*Pilot.*pilot_mask;
       
     
        L=length(z);
        L2=length(z2);
        
        Nbps=log2(M_ary);

%         z=[z zeros(1,Nbps-mod(length(z),Nbps))];
%         z2=[z2 zeros(1,Nbps-mod(length(z2),Nbps))];
          
        s=flipud(reshape(z,Nbps,length(z)/Nbps));
        s2=flipud(reshape(z2,Nbps,length(z2)/Nbps));
        y=bi2de(s');
        y2=bi2de(s2');
        y=y';
        y2=y2';
%         L=length(y);
%         L2=length(y2);  
%             if L<L2
%             y=[y randi(M_ary-1,1,L2-L)];   %them bot ve do dai
%             else
%             y2=[y2 randi(M_ary-1,1,L-L2)]; %them bot ve do dai
%             end
        d=mod(length(y),N_D);
        if d~=0
            zz=zeros(1,N_D-d);
            y=[y zz];
  
            y2=[y2 zz];
        end
        dataluu=[y; y2];
        save('mimo_data.mat','dataluu');
         
        y=qammod(y,M_ary);
        y2=qammod(y2,M_ary);
       
        L=length(y)/N_D;

        data=[];
        data2=[];
        for i=1:L
            data=[data;y((i-1)*N_D+1:i*N_D)];
            data2=[data2;y2((i-1)*N_D+1:i*N_D)];
        end
        [m,n]=size(data)


        %----------Zeros Insert  -----------------
         
        zL=zeros(m,SubL);
        zH=zeros(m,NFFT-SubH);

        data =[zL data zH];
        data2=[zL data2 zH];
 
        %------------------------
        % Preparing pilot pattern 
        %------------------------
        %%%%%%%%%%%%
        dataP=[];
        dataP2=[];
        
        for i=1:m
            if mod(i,D_t)==1 
                dataP=[dataP;Pilot];
                dataP=[dataP;data(i,:)];
                
                dataP2=[dataP2;Pilot2];
                dataP2=[dataP2;data2(i,:)];
            else
                dataP=[dataP;data(i,:)];
                dataP2=[dataP2;data2(i,:)];
            end
        end

        [mp,np]=size(dataP)
        %------------------------
        % calculate  IFFT for OFDM Symbol combination with  Pulse shaping after that
        %------------------------

        dataIFFT=[];
        dataIFFT2=[];
        for i=1:mp
        %    dataIFFT=[dataIFFT; ifft(dataP(i,:))];
            dataIFFT =[dataIFFT;  ifft([0 dataP(i,:)  fliplr(conj(dataP(i,:)))])];
            dataIFFT2=[dataIFFT2; ifft([0 dataP2(i,:) fliplr(conj(dataP2(i,:)))])];
        end
 
        [mp,np]=size(dataIFFT);

        %------------------------
        % Copy GI
        %------------------------

        dataGI=[dataIFFT(:,np-GI+1:np) dataIFFT];
        dataGI2=[dataIFFT2(:,np-GI+1:np) dataIFFT2];
%  insert more pilot in the header  for sys
       dataGI=[dataGI(1,:);dataGI];
       dataGI2=[dataGI2(1,:);dataGI2];
       [mm nn]=size(dataGI)
        set(handles.edit37,'String',num2str(mm));
        dataGI=dataGI.';
        dataGI2=dataGI2.';

        dd=dataGI(:);
        yy=real(dd);
        
        dd2=dataGI2(:);
        yy2=real(dd2);
      
%         f1=(ofdm_f1+ofdm_f2)/2;
%         tttt=0:2*pi/(flm/f1):100000*pi;
%         xxxx=sin(tttt)'/2;
%         yy=[xxxx(1:2*(NFFT+GI)); NFFT*yy/N_D];
        yy=2*NFFT*yy/N_D;  % chia 6 cho 1024 muc-chia 4 cho 256 muc - chia 2 cho 64 muc
        yy2=2*NFFT*yy2/N_D;
        yy=yy(1:length(yy)-K);
        yy2=yy2(1:length(yy2)-K);
         %%%% tao kenh mo phong
%          L_channel=16;
%          h11=(randn(1,L_channel)+j*randn(1,L_channel))/sqrt(2);
%          h12=(randn(1,L_channel)+j*randn(1,L_channel))/sqrt(2);
%          h21=(randn(1,L_channel)+j*randn(1,L_channel))/sqrt(2);
%          h22=(randn(1,L_channel)+j*randn(1,L_channel))/sqrt(2);
%          yy_sim=[zeros(20000,1); conv(yy,h11)+conv(yy2,h12); zeros(20000,1)];
%          yy_sim2=[zeros(20000,1); conv(yy,h21)+conv(yy2,h22); zeros(20000,1)];
%          n1=randn(length(yy_sim),1)+j*randn(length(yy_sim),1);
%          n2=randn(length(yy_sim),1)+j*randn(length(yy_sim),1);
%          yy_sim=yy_sim+n1/10;
%          yy_sim2=yy_sim2+n2/10;
% 
%          wavwrite([yy_sim yy_sim2],flm,16,'mimo_ofdm_channel_sim.wav');
        wavwrite([yy yy2],flm,16,'mimo_ofdm_BER_12_15_4_256_1024.wav');
        Khuech_dai=get(handles.edit39,'String');
        Khuech_dai=str2double(Khuech_dai);
        yy=Khuech_dai*yy;
        yy2=Khuech_dai*yy2;
        sound([yy yy2],flm);
        
         size([yy yy2])
 
        axes(handles.axes1);
        plot(yy);
        axes(handles.axes10);
        plot(yy2);
        set(handles.edit38,'String',num2str(length(yy)/flm));
        
         axes(handles.axes9);
         yy=abs(fft(yy(1:length(yy))));
         plot(yy(1:length(yy)/2));
         
         axes(handles.axes11);
         yy2=abs(fft(yy2(1:length(yy2))));
         plot(yy2(1:length(yy2)/2));
 
end

function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop,'Value',0); 
set(handles.disp_txt,'String','Show Text');
set(handles.edit_FS,'String',num2str(192000));
set(handles.edit_F0,'String',num2str(14300));
set(handles.edit_df,'String',num2str(300));
set(handles.edit_bitrates,'String',num2str(100));
set(handles.load,'String','Load');
path=0;
savefile='pqfile.mat';
save(savefile,'path');
cla(handles.axes1);

function disp_txt_Callback(hObject, eventdata, handles)

function disp_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disp_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_FS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FS as text
%        str2double(get(hObject,'String')) returns contents of edit_FS as a double
nummod = str2double(get(hObject,'String'));
if isnan(nummod) || ~isreal(nummod)  
    % isdouble returns NaN for non-numbers and nummod cannot be complex
    % Disable the Transfer button and change its string to say why
    set(handles.transfer,'String','Cannot transfer')
    set(handles.transfer,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Transfer button with its original name
    set(handles.transfer,'String','Transfer')
    set(handles.transfer,'Enable','on')
end

function edit_FS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_df_Callback(hObject, eventdata, handles)
% hObject    handle to edit_df (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_df as text
%        str2double(get(hObject,'String')) returns contents of edit_df as a double
numpro = str2double(get(hObject,'String'));
if isnan(numpro) || ~isreal(numpro)  
    % isdouble returns NaN for non-numbers and numpro cannot be complex
    % Disable the Transfer button and change its string to say why
    set(handles.transfer,'String','Cannot transfer')
    set(handles.transfer,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Transfer button with its original name
    set(handles.transfer,'String','Transfer')
    set(handles.transfer,'Enable','on')
end

function edit_df_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_bitrates_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bitrates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bitrates as text
%        str2double(get(hObject,'String')) returns contents of edit_bitrates as a double
numfft = str2double(get(hObject,'String'));
if isnan(numfft) || ~isreal(numfft)
    % isdouble returns NaN for non-numbers and numfft cannot be complex
    % Disable the Transfer button and change its string to say why
    set(handles.transfer,'String','Cannot transfer')
    set(handles.transfer,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Transfer button with its original name
    set(handles.transfer,'String','Transfer')
    set(handles.transfer,'Enable','on')
end

% --- Executes during object creation, after setting all properties.
function edit_bitrates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Continous_Tx.
function Continous_Tx_Callback(hObject, eventdata, handles)
% hObject    handle to Continous_Tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
val = get(handles.popupmenu1,'Value');

set(handles.stop,'Value',1); 

l= get(handles.axes1, 'Children'); 
y=get(l,'Ydata');

l2= get(handles.axes10, 'Children'); 
y2=get(l2,'Ydata');

pl=get(handles.stop,'Value');

yy=[];
yy2=[];

while pl 
play(p);
pause(.01);
pl=get(handles.stop,'Value');
end
stop(p)
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

 
    

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in repeat_transfer.
function repeat_transfer_Callback(hObject, eventdata, handles)
% hObject    handle to repeat_transfer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      % Number of QASK modulation level
flm=get(handles.edit_FS,'String');       % Number of QASK modulation level
flm=str2double(flm); 
disp (flm);
 
l= get(handles.axes1, 'Children'); 
y=get(l,'Ydata');
 
 
l2= get(handles.axes10, 'Children'); 
y2=get(l2,'Ydata');
 
 sound([y' y2'],flm);
 
 
 

function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stop,'Value',0);
 



function time_period_Callback(hObject, eventdata, handles)
% hObject    handle to time_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_period as text
%        str2double(get(hObject,'String')) returns contents of time_period as a double


% --- Executes during object creation, after setting all properties.
function time_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text34.
function text34_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function command_repeat_Callback(hObject, eventdata, handles)
% hObject    handle to command_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of command_repeat as text
%        str2double(get(hObject,'String')) returns contents of command_repeat as a double


% --- Executes during object creation, after setting all properties.
function command_repeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to command_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_FS_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_FS as text
%        str2double(get(hObject,'String')) returns contents of ofdm_FS as a double


% --- Executes during object creation, after setting all properties.
function ofdm_FS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_F1_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_F1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_F1 as text
%        str2double(get(hObject,'String')) returns contents of ofdm_F1 as a double


% --- Executes during object creation, after setting all properties.
function ofdm_F1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_F1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_F2_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_F2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_F2 as text
%        str2double(get(hObject,'String')) returns contents of ofdm_F2 as a double


% --- Executes during object creation, after setting all properties.
function ofdm_F2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_F2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_mod_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_mod as text
%        str2double(get(hObject,'String')) returns contents of ofdm_mod as a double


% --- Executes during object creation, after setting all properties.
function ofdm_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_GI_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_GI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_GI as text
%        str2double(get(hObject,'String')) returns contents of ofdm_GI as a double


% --- Executes during object creation, after setting all properties.
function ofdm_GI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_GI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NFFT_Callback(hObject, eventdata, handles)
% hObject    handle to NFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NFFT as text
%        str2double(get(hObject,'String')) returns contents of NFFT as a double


% --- Executes during object creation, after setting all properties.
function NFFT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_SubL_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_SubL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_SubL as text
%        str2double(get(hObject,'String')) returns contents of ofdm_SubL as a double


% --- Executes during object creation, after setting all properties.
function ofdm_SubL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_SubL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ofdm_SubH_Callback(hObject, eventdata, handles)
% hObject    handle to ofdm_SubH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ofdm_SubH as text
%        str2double(get(hObject,'String')) returns contents of ofdm_SubH as a double


% --- Executes during object creation, after setting all properties.
function ofdm_SubH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ofdm_SubH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_size_Callback(hObject, eventdata, handles)
% hObject    handle to frame_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_size as text
%        str2double(get(hObject,'String')) returns contents of frame_size as a double


% --- Executes during object creation, after setting all properties.
function frame_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse_width_Callback(hObject, eventdata, handles)
% hObject    handle to pulse_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse_width as text
%        str2double(get(hObject,'String')) returns contents of pulse_width as a double


% --- Executes during object creation, after setting all properties.
function pulse_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function disp_txt2_Callback(hObject, eventdata, handles)
% hObject    handle to disp_txt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of disp_txt2 as text
%        str2double(get(hObject,'String')) returns contents of disp_txt2 as a double


% --- Executes during object creation, after setting all properties.
function disp_txt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disp_txt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
