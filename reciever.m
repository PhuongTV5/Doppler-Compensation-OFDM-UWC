function varargout = reciever(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reciever_OpeningFcn, ...
                   'gui_OutputFcn',  @reciever_OutputFcn, ...
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

function reciever_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);

function varargout = reciever_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function load_data_Callback(hObject, eventdata, handles)
[filename,pathname]=uigetfile({'*wav'});
path=[pathname,filename];
if (filename==0)
    set(handles.load_data,'String','Cannot load');
    pause(1);
    set(handles.load_data,'String','Load wav file');
else
    y=wavread(path);

axes(handles.axes3);
plot(y);  
  
    kk=strfind(filename,'_');

    flm=filename(5:kk(1)-1);
    set(handles.edit_FS,'String',flm); 

    f1=filename(kk(1)+1:kk(2)-1);
    set(handles.ofdm_f1,'String',f1); 
    
    f2=filename(kk(2)+1:kk(3)-1);
    set(handles.ofdm_f2,'String',f2); 

    M_ary=filename(kk(3)+1:kk(4)-1);
    set(handles.ofdm_mod,'String',M_ary); 
    
   GI=filename(kk(4)+1:kk(5)-1);
    set(handles.ofdm_GI,'String',GI); 
    
    NFFT=filename(kk(5)+1:kk(6)-1);
    set(handles.NFFT,'String',NFFT);
  
end

function reciever_Callback(hObject, eventdata, handles)

        flm=get(handles.FS,'String');      
        flm=str2double(flm);
        
        rt=get(handles.edit_time,'String');       % Number of QASK modulation level
        rt=str2double(rt);

        ofdm_f1=get(handles.ofdm_F1,'String');   
        ofdm_f1=str2double(ofdm_f1);    

        ofdm_f2=get(handles.ofdm_F2,'String');   
        ofdm_f2=str2double(ofdm_f2);

        GI=get(handles.ofdm_GI,'String');
        GI=str2double(GI);   % guard interval

        NFFT=get(handles.NFFT,'String');
        NFFT=str2double(NFFT); 

        M_ary=get(handles.ofdm_mod,'String');
        M_ary=str2double(M_ary); % Number of QAM modulation level
        
        D_t=get(handles.frame_size,'String');
        D_t=str2double(D_t);  %Distance in time domain for Pilot insertion
        
        SubL=floor(2*ofdm_f1*NFFT/flm);
        SubL=SubL-mod(SubL,2);
        SubH=ceil(2*ofdm_f2*NFFT/flm);
        SubH=SubH+mod(SubH,2);
        
        set(handles.ofdm_SubL,'String',num2str(SubL)); 
        set(handles.ofdm_SubH,'String',num2str(SubH));  
        
        N_D = SubH- SubL;  % The number  of subcarrier data
        PO1=[ones(1,N_D/2); zeros(1,N_D/2)];
        PO1=PO1(:);
        PO2=[zeros(1,N_D/2); ones(1,N_D/2)];
        PO2=PO2(:);
         
        pilot_mask1=[zeros(1,SubL) PO1' zeros(1,NFFT-SubH)];
        pilot_mask2=[zeros(1,SubL) PO2' zeros(1,NFFT-SubH)];
%record  data
    rt=rt+.4;
    r=audiorecorder(flm,16,2);
    record(r);
    pause(rt);
    z=getaudiodata(r,'double');
    y=z(:,1);
    y2=z(:,2);
    

    %bo loc
   [B A]=butter(3,[1.8*ofdm_f1/flm 2.2*ofdm_f2/flm]);
    y=filter(B,A,y);
    axes(handles.axes3);
    plot(y);
    y2=filter(B,A,y2); 
    axes(handles.axes13);
    plot(y2);

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
        
        Pilot2=sqrt(M_ary)*Pilot.*pilot_mask2;
        Pilot1=sqrt(M_ary)*Pilot.*pilot_mask1;
 
        Nbps=log2(M_ary);    

%%%%%%%%%%%%%%%%%%
 
 ss=find(y>max(y)/2);
 K=2*NFFT+GI+1; 
 Ldata=ss(length(ss))-ss(1);
 if ss(1)<K 
     ss(1)=K+1;
 end

    yy=y(ss(1)-K:ss(1)+7*K);
          syn=[];
         L=length(yy)
        for i=1:L-2*NFFT-GI-2
               syn=[syn sum(abs(yy(i:i+GI)-yy(i+K-GI:i+K)))];
        end
            syn=max(syn)-syn;
%          axes(handles.axes15);
%          plot(syn)
         syn2=[];
           for i=1:L-2*NFFT-GI-2
              syn2=[syn2 abs(yy(i:i+GI)'*yy(i+K-GI:i+K))];
           end
           
         axes(handles.axes15);
         syn2=syn.*syn2;
            %%%%%%%%%%%%%%%%%%%% Schmild method use to compare
              PD=[];
        for i=1:L-2*NFFT-GI-2
               PD=[PD sum(yy(i:i+GI).*yy(i+K-GI:i+K))];
        end
          RD=[];
        for i=1:L-2*NFFT-GI-2
               RD=[RD sum(yy(i+K-GI:i+K).^2)];
        end
        S_syn=PD./RD;
        Header_text=1:1:length(syn2);
        plot(Header_text,syn2/max(syn2),'-',Header_text,S_syn,'-.')
%%%%%%%%%%%%%%%%%%%%
             
             
     bb=find(syn2>max(syn2)/3);
     [ii dd]=max(syn2(1:bb(1)+K/2))
     
       Noise=y(ss(length(ss))+2*K:end);
       NL=length(Noise)
       Noise_Power=mean(Noise.^2);
       set(handles.Noise_Power,'String',num2str(log10(Noise_Power)));
       start_point=dd;
       y=y(start_point+ss(1)-K:start_point+ss(1)-K+Ldata);
       y2=y2(start_point+ss(1)-K:start_point+ss(1)-K+Ldata);
       Signal_Power=mean(y.^2);
       set(handles.CNR,'String',num2str(log10(Signal_Power/Noise_Power)));

       L=length(y)/K;
       Header_text=1:1:length(syn2);
       plot(Header_text,syn2/max(syn2),'-',Header_text,S_syn,'-.')

y=y.';
y2=y2.';
dataGI=[];
dataGI2=[];
for i=1:L
    dataGI=[dataGI;y((i-1)*K+1:i*K)];
    dataGI2=[dataGI2;y2((i-1)*K+1:i*K)];
end
%-----------Remove GI
[m,n]=size(dataGI)
%%%%%%%%%%%%%%%%%%
dataP=dataGI(:,GI+1:n);
dataP=dataP(2:m,:);
%%%%%%%%%%%%%%%%%%
dataP2=dataGI2(:,GI+1:n);
dataP2=dataP2(2:m,:);
%------------
[mp,np]=size(dataP);

dataFFT=[];
for i=1:mp
    if mod(i,D_t+1)==1 
        Pr1= dataP(i,:);
        Pr2= dataP2(i,:) ;
    H11=fft(Pr1)./([0 Pilot1 fliplr(conj(Pilot1))]);
    H11=H11(SubL+1:SubH);
    H11=H11(2:2:length(H11));
    H11=interp(H11,2);
    
    H12=fft(Pr2)./([0 Pilot1 fliplr(conj(Pilot1))]);
    H12=H12(SubL+1:SubH);
    H12=H12(2:2:length(H12));
    H12=interp(H12,2);
    
    H21=fft(Pr1)./([0 Pilot2 fliplr(conj(Pilot2))]);
    H21=H21(SubL+2:SubH+1);
    H21=H21(2:2:length(H21));
    H21=interp(H21,2);
    
    H22=fft(Pr2)./([0 Pilot2 fliplr(conj(Pilot2))]);
    H22=H22(SubL+2:SubH+1);
    H22=H22(2:2:length(H22));
    H22=interp(H22,2);
    else
         RFFT=fft(dataP(i,:));
         RFFT=RFFT(SubL+2:SubH+1);
         R2FFT=fft(dataP2(i,:));
         R2FFT=R2FFT(SubL+2:SubH+1);
         R_data=[];
         for k=1:N_D
             H=[H11(k) H21(k);H12(k) H22(k)];
             R_data=[R_data inv(H)*[RFFT(k);R2FFT(k)]];
         end
        dataFFT=[dataFFT R_data];
    end
end

    yy=dataFFT(1,:);
    yy2=dataFFT(2,:);
%------------------------------------
    yy=yy.';
    yy=yy(:);
    yy2=yy2.';
    yy2=yy2(:);
    goc=atan(imag(yy)./real(yy));
    goc_duong=goc>0;
    goc=sum(goc.*goc_duong)/sum(goc_duong)
    goc2=atan(imag(yy2)./real(yy2));
    goc_duong2=goc2>0;
    goc2=sum(goc2.*goc_duong2)/sum(goc_duong2)
    yy=yy*exp(-j*pi/25);
    yy2=yy2*exp(-j*pi/25);
    axes(handles.axes2)
    scatter(real(yy),imag(yy),'.');
    axes(handles.axes14)
    scatter(real(yy2),imag(yy2),'.');
    
%------Demodulation
 yy=qamdemod(yy,M_ary);
 yy2=qamdemod(yy2,M_ary);

 load('mimo_data.mat');

%-----------------
datathu=[yy yy2]';
       if length(datathu)>length(dataluu) 
           datathu=datathu(:,1:length(dataluu));
            
       else
           dataluu=dataluu(:,1:length(datathu));
            
       end
       [nss1 err1]=symerr(dataluu,datathu)
      set(handles.edit_SER,'String',num2str(err1));
 
%--------------Convert to ASCII
L=length(yy);
y=[];
y2=[];
for i=1:L
    y=[y dec2bin(yy(i),Nbps)];
    y2=[y2 dec2bin(yy2(i),Nbps)];
end

t=num2str(y');
t2=num2str(y2');
st=[];
st2=[];
for i=1:length(t)/8
    st=[st char(bin2dec(t((i-1)*8+1:i*8)'))];
    st2=[st2 char(bin2dec(t2((i-1)*8+1:i*8)'))];
end
File=st;
File2=st2;
set(handles.disp_txt1,'String',File);
set(handles.disp_txt2,'String',File2);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes13);
cla(handles.axes14);
cla(handles.axes2);
cla(handles.axes3);
clc;
set(handles.disp_txt1,'String','Show text');
set(handles.disp_txt2,'String','Show text');
set(handles.edit_mod,'String',num2str(16));
set(handles.edit_pro,'String',num2str(16));
set(handles.edit_fft,'String',num2str(64));


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
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



function FS_Callback(hObject, eventdata, handles)
% hObject    handle to FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FS as text
%        str2double(get(hObject,'String')) returns contents of FS as a double


% --- Executes during object creation, after setting all properties.
function FS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time as text
%        str2double(get(hObject,'String')) returns contents of edit_time as a double


% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Noise_Power_Callback(hObject, eventdata, handles)
% hObject    handle to Noise_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Noise_Power as text
%        str2double(get(hObject,'String')) returns contents of Noise_Power as a double


% --- Executes during object creation, after setting all properties.
function Noise_Power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Noise_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CNR_Callback(hObject, eventdata, handles)
% hObject    handle to CNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CNR as text
%        str2double(get(hObject,'String')) returns contents of CNR as a double


% --- Executes during object creation, after setting all properties.
function CNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SER_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SER as text
%        str2double(get(hObject,'String')) returns contents of edit_SER as a double


% --- Executes during object creation, after setting all properties.
function edit_SER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edit_VBLAST.
function edit_VBLAST_Callback(hObject, eventdata, handles)
% hObject    handle to edit_VBLAST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edit_VBLAST
