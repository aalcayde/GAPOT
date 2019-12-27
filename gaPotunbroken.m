function varargout = gaPotunbroken(varargin)
% GAPOTUNBROKEN MATLAB code for gaPotunbroken.fig
%      GAPOTUNBROKEN, by itself, creates a new GAPOTUNBROKEN or raises the existing
%      singleton*.
%
%      H = GAPOTUNBROKEN returns the handle to a new GAPOTUNBROKEN or the handle to
%      the existing singleton*.
%
%      GAPOTUNBROKEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAPOTUNBROKEN.M with the given input arguments.
%
%      GAPOTUNBROKEN('Property','Value',...) creates a new GAPOTUNBROKEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gaPotunbroken_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gaPotunbroken_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gaPotunbroken

% Last Modified by GUIDE v2.5 17-Dec-2019 14:50:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gaPotunbroken_OpeningFcn, ...
                   'gui_OutputFcn',  @gaPotunbroken_OutputFcn, ...
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


% --- Executes just before gaPotunbroken is made visible.
function gaPotunbroken_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gaPotunbroken (see VARARGIN)

% Choose default command line output for gaPotunbroken
handles.output = hObject;
ptd = PowerTimeDomain;
handles.ptd = ptd;
handles.ColumnName = {'Values'};  
set(handles.uitResults,'data',[],'ColumnName',{'Values'});
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gaPotunbroken wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gaPotunbroken_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function VoltageVector_Callback(hObject, eventdata, handles)
% hObject    handle to VoltageVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VoltageVector as text
%        str2double(get(hObject,'String')) returns contents of VoltageVector as a double


% --- Executes during object creation, after setting all properties.
function VoltageVector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoltageVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrentVector_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentVector as text
%        str2double(get(hObject,'String')) returns contents of CurrentVector as a double


% --- Executes during object creation, after setting all properties.
function CurrentVector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SimulationTime_Callback(hObject, eventdata, handles)
% hObject    handle to SimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SimulationTime as text
%        str2double(get(hObject,'String')) returns contents of SimulationTime as a double


% --- Executes during object creation, after setting all properties.
function SimulationTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppmShowPlot.
function ppmShowPlot_Callback(hObject, eventdata, handles)
% hObject    handle to ppmShowPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppmShowPlot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppmShowPlot
ptd = handles.ptd;
contents = get(handles.ppmShowPlot,'value');
Fs = ptd.Fs;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(ptd.SimulationTime);             % Length of signal
switch contents
    case 1 
        
        plot(handles.axes1,ptd.SimulationTime,ptd.Va,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Voltage (V)')
        fftVa = fft(ptd.Va);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Va(f)|')
        set(handles.axes2, 'XScale', 'log')
    case 2
        
        plot(handles.axes1,ptd.SimulationTime,ptd.HVa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Voltage (V)')
        fftVa = fft(ptd.Va);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Va(f)|')
        set(handles.axes2, 'XScale', 'log')
    case 3
        disp('zero')
    case 7
        
        plot(handles.axes1,ptd.SimulationTime,ptd.Ia,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Ia);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ia(f)|')
        set(handles.axes2, 'XScale', 'log')
     case 8
        
        plot(handles.axes1,ptd.SimulationTime,ptd.HIa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.HIa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ia(f)|')
        set(handles.axes2, 'XScale', 'log')
      case 9       
        plot(handles.axes1,ptd.SimulationTime,ptd.Ipa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Ipa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ipa(f)|')
        set(handles.axes2, 'XScale', 'log')
      case 10
        plot(handles.axes1,ptd.SimulationTime,ptd.Iqa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Iqa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Iqa(f)|')
        set(handles.axes2, 'XScale', 'log')
      case 11
        plot(handles.axes1,ptd.SimulationTime,ptd.Ifa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Ifa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ifa(f)|')
        set(handles.axes2, 'XScale', 'log')
       case 12
        plot(handles.axes1,ptd.SimulationTime,ptd.Ixa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Ixa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ixa(f)|')
        set(handles.axes2, 'XScale', 'log')
       case 13
        plot(handles.axes1,ptd.SimulationTime,ptd.Ibra,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Ibra);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Ibra(f)|')
        set(handles.axes2, 'XScale', 'log')
       case 19
        plot(handles.axes1,ptd.SimulationTime,ptd.Mpa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Mpa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Mpa(f)|')
        set(handles.axes2, 'XScale', 'log')
       case 20
        plot(handles.axes1,ptd.SimulationTime,ptd.Mqa,'LineWidth',2);
        xlabel(handles.axes1,'Time (s)')
        ylabel(handles.axes1,'Current (A)')
        fftVa = fft(ptd.Mqa);
        P2 = abs(fftVa/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 20*log10(2*P1(2:end-1));
        f = Fs*(0:(L/2))/L;
        plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
%         title('Single-Sided Amplitude Spectrum of Va(t)')
        xlabel(handles.axes2,'f (Hz)')
        ylabel(handles.axes2,'|Mqa(f)|')
        set(handles.axes2, 'XScale', 'log')
       case 25
        plot(handles.axes1,ptd.Mpa,ptd.Mqa,'LineWidth',2);
        xlabel(handles.axes1,' Mp ')
        ylabel(handles.axes1,' Mq ')
%         fftVa = fft(ptd.Mqa);
%         P2 = abs(fftVa/L);
%         P1 = P2(1:L/2+1);
%         P1(2:end-1) = 20*log10(2*P1(2:end-1));
%         f = Fs*(0:(L/2))/L;
%         plot(handles.axes2,f(1:end-2),P1(1:end-2)) 
% %         title('Single-Sided Amplitude Spectrum of Va(t)')
%         xlabel(handles.axes2,'f (Hz)')
%         ylabel(handles.axes2,'|Mqa(f)|')
%         set(handles.axes2, 'XScale', 'log')
        case 26
            plot(handles.axes1,ptd.Ipa,ptd.Iqa,'LineWidth',2);
            xlabel(handles.axes1,' Ipa ')
            ylabel(handles.axes1,' Iqa ')
        case 27
            plot(handles.axes1,ptd.Ifa,ptd.Ibra,'LineWidth',2);
            xlabel(handles.axes1,' Ifa ')
            ylabel(handles.axes1,' Ibra ')
    otherwise
        disp('other value')
end

% --- Executes during object creation, after setting all properties.
function ppmShowPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppmShowPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tbHoldOnOff.
function tbHoldOnOff_Callback(hObject, eventdata, handles)
% hObject    handle to tbHoldOnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbHoldOnOff

if ( get(hObject,'Value') == true)
    hold( handles.axes1, 'on' )
    hold( handles.axes2, 'on' )
else
    hold( handles.axes1, 'off' )
    hold( handles.axes2, 'off' )
end

% --- Executes on button press in pbCalculate.
function pbCalculate_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ptd = handles.ptd;
clifford_signature(6,0)
ColumnName = handles.ColumnName;
if ptd. LoadedDataSinglePhase == true            
            ptd.Varms = rms(ptd.Va);
            ptd.HVa = imag(hilbert(ptd.Va));
            ptd.HIa = imag(hilbert(ptd.Ia));
            ptd.gaVa =  ptd.Va*e1+ptd.HVa*e2;
            ptd.gaV2a = ptd.gaVa.*ptd.gaVa;
            ptd.gaIa =  ptd.Ia*e1+ptd.HIa*e2;
            ptd.gaMa = ptd.gaVa.*ptd.gaIa;
            ptd.Mpa = grade(ptd.gaMa,0);
            ptd.Mqa = ptd.gaMa-ptd.Mpa;
            ptd.Ipa =  ptd.gaVa./ ptd.gaV2a.*ptd.Mpa;
            ptd.Iqa = ptd.gaIa-ptd.Ipa;
            ptd.Ifa = part(mean(ptd.Mpa)/2,1)/(ptd.Varms^2).*ptd.gaVa;
            ptd.Ixa = ptd.Ipa-ptd.Ifa;
            ptd.Ibra = part(mean(ptd.Mqa)/2,8)/(ptd.Varms^2).*ptd.gaVa;
            ptd.Mpa = part(ptd.Mpa,1);
            ptd.Mqa = part(ptd.Mqa,8);
            ptd.Ipa=part(ptd.Ipa,2);
            ptd.Iqa=part(ptd.Iqa,2);
            ptd.Ifa=part(ptd.Ifa,2);
            ptd.Ixa=part(ptd.Ixa,2);
            ptd.Ibra=part(ptd.Ibra,3);
            ptd.mMpa = round(mean(ptd.Mpa),2);
            ptd.mMqa = round(mean(ptd.Mqa),2);
            ptd.P = ptd.mMpa/2;
            ptd.Q = ptd.mMqa/2;
             ptd.Iarms = round(rms(ptd.Ia),2);
            ptd.Iparms = round(rms(ptd.Ipa),2);
            ptd.Iqarms = round(rms(ptd.Iqa),2);
            ptd.Ifarms = round(rms(ptd.Ifa),2);
            ptd.Ixarms = round(rms(ptd.Ixa),2);
            ptd.Ibrarms = round(rms(ptd.Ibra),2);
            ptd.LoadedDataSinglePhase = false;
            handles.ptd = ptd;

            Values = [ptd.Varms; ptd.Iarms;0; ptd.Iparms;ptd.Iqarms; ptd.Ifarms; ptd.Ixarms; ptd.Ibrarms; ptd.P;ptd.Q];
            set(handles.uitResults,'data',Values,'ColumnName',ColumnName);
            set(handles.ppmShowPlot,'Enable','on');
            assignin('base','simulation',ptd)
else
     
    
    paramTime = get(handles.SimulationTime,'String');
    frequency = str2num(get(handles.edFrequency,'String'));
    ptd.SimulationTime = str2num(paramTime);
    ptd.W = 2*pi*frequency;
    Fs = strsplit(paramTime,':');
    ptd.Fs = 1/str2num(Fs{2});

    A = get(handles.VoltageVector,'String');
    B = strsplit(A,';');
    paramVa = str2num(strcat(B{1},']'));

    C = get(handles.CurrentVector,'String');
    D = strsplit(C,';');
    paramIa = str2num(strcat(D{1},']'));

    E = get(handles.edVoltagePhase,'String');
    F = strsplit(E,';');
    paramPhVa = str2num(strcat(F{1},']'));

    G = get(handles.edCurrentPhase,'String');
    H = strsplit(G,';');
    paramPhIa = str2num(strcat(H{1},']'));

    try
      paramVb = str2num(strcat('[',strcat(B{2},']')));
      paramVc = str2num(strcat('[',B{3}));
      paramPhVc = str2num(strcat('[',strcat(F{2},']')));
      paramPhVd = str2num(strcat('[',F{3}));
      paramIb = str2num(strcat('[',strcat(D{2},']')));
      paramIc = str2num(strcat('[',D{3}));
      paramPhIc = str2num(strcat('[',strcat(G{2},']')));
      paramPhId = str2num(strcat('[',G{3}));
      ptd.SinglePhase = false;
    catch ME

    %    if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
    %       msg = ['Dimension mismatch occurred: First argument has ', ...
    %             num2str(size(A,2)),' columns while second has ', ...
    %             num2str(size(B,2)),' columns.'];
    %         causeException = MException('MATLAB:myCode:dimensions',msg);
    %         ME = addCause(ME,causeException);
    %    end
    %    rethrow(ME)

    ptd.SinglePhase = true;

    A = get(handles.VoltageVector,'String');
    B = strsplit(A,';');
    paramVa = str2num(B{1});

    C = get(handles.CurrentVector,'String');
    D = strsplit(C,';');
    paramIa = str2num(D{1});

    E = get(handles.edVoltagePhase,'String');
    F = strsplit(E,';');
    paramPhVa = str2num(F{1});

    G = get(handles.edCurrentPhase,'String');
    H = strsplit(G,';');
    paramPhIa = str2num(H{1});

    end 

    switch ptd.SinglePhase
        case true
            Va = zeros(length(ptd.SimulationTime),1);
            Ia = zeros(length(ptd.SimulationTime),1);
            
%             for i = 1:length(ptd.SimulationTime)
%                for j = 1:length(paramVa)
%                     Va(i) = Va(i) + paramVa(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhVa(j));           
%                end
%                for j = 1:length(paramIa)
%                     Ia(i) = Ia(i) + paramIa(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhIa(j));           
%                end
%             end
            
            for i = 1:length(ptd.SimulationTime)
               for j = 1:length(paramVa)/2
                    Va(i) = Va(i) + paramVa(2*j-1)*sin(ptd.W*(j)*ptd.SimulationTime(i))+ paramVa(2*j)*cos(ptd.W*(j)*ptd.SimulationTime(i));           
               end
               for j = 1:length(paramIa)/2
                    Ia(i) = Ia(i) + paramIa(2*j-1)*sin(ptd.W*(j)*ptd.SimulationTime(i))+ paramIa(2*j)*cos(ptd.W*(j)*ptd.SimulationTime(i));             
               end
            end
            
            
            ptd.Va = sqrt(2)*Va;
            ptd.Ia = sqrt(2)*Ia;
            ptd.Varms = rms(ptd.Va);
            ptd.HVa = imag(hilbert(ptd.Va));
            ptd.HIa = imag(hilbert(ptd.Ia));
            ptd.gaVa =  ptd.Va*e1+ptd.HVa*e2;
            ptd.gaV2a = ptd.gaVa.*ptd.gaVa;
            ptd.gaIa =  ptd.Ia*e1+ptd.HIa*e2;
            ptd.gaMa = ptd.gaVa.*ptd.gaIa;
            ptd.Mpa = grade(ptd.gaMa,0);
            ptd.Mqa = ptd.gaMa-ptd.Mpa;
            ptd.Ipa =  ptd.gaVa./ ptd.gaV2a.*ptd.Mpa;
            ptd.Iqa = ptd.gaIa-ptd.Ipa;
            ptd.Ifa = part(mean(ptd.Mpa)/2,1)/(ptd.Varms^2).*ptd.gaVa;
            ptd.Ixa = ptd.Ipa-ptd.Ifa;
            ptd.Ibra = part(mean(ptd.Mqa)/2,8)/(ptd.Varms^2).*ptd.gaVa;
            ptd.Mpa = part(ptd.Mpa,1);
            ptd.Mqa = part(ptd.Mqa,8);
            ptd.Ipa=part(ptd.Ipa,2);
            ptd.Iqa=part(ptd.Iqa,2);
            ptd.Ifa=part(ptd.Ifa,2);
            ptd.Ixa=part(ptd.Ixa,2);
            ptd.Ibra=part(ptd.Ibra,3);
            ptd.mMpa = round(mean(ptd.Mpa),2);
            ptd.mMqa = round(mean(ptd.Mqa),2);
            ptd.P = ptd.mMpa/2;
            ptd.Q = ptd.mMqa/2;
            ptd.Iarms = round(rms(ptd.Ia),2);
            ptd.Iparms = round(rms(ptd.Ipa),2);
            ptd.Iqarms = round(rms(ptd.Iqa),2);
            ptd.Ifarms = round(rms(ptd.Ifa),2);
            ptd.Ixarms = round(rms(ptd.Ixa),2);
            ptd.Ibrarms = round(rms(ptd.Ibra),2);
            handles.ptd = ptd;

            Values = [ptd.Varms; ptd.Iarms;0; ptd.Iparms;ptd.Iqarms; ptd.Ifarms; ptd.Ixarms; ptd.Ibrarms; ptd.P;ptd.Q];
            set(handles.uitResults,'data',Values,'ColumnName',ColumnName);
            set(handles.ppmShowPlot,'Enable','on');   
            assignin('base','simulation',ptd)
        otherwise
            Va = zeros(length(ptd.SimulationTime),1);
            Vb = zeros(length(ptd.SimulationTime),1);
            Vc = zeros(length(ptd.SimulationTime),1);
            Ia = zeros(length(ptd.SimulationTime),1);
            Ib = zeros(length(ptd.SimulationTime),1);
            Ic = zeros(length(ptd.SimulationTime),1);
            for i = 1:length(ptd.SimulationTime)
               for j = 1:length(paramVa)
                    Va(i) = Va(i) + paramVa(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhVa(j));           
               end
               for j = 1:length(paramVb)
                    Vb(i) = Vb(i) + paramVb(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhVb(j));           
               end
               for j = 1:length(paramVc)
                    Vc(i) = Vc(i) + paramVc(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhVc(j));           
               end
               for j = 1:length(paramIa)
                    Ia(i) = Ia(i) + paramIa(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhIa(j));           
               end
               for j = 1:length(paramIb)
                    Ib(i) = Ib(i) + paramIb(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhIb(j));           
               end
               for j = 1:length(paramIc)
                    Ic(i) = Ic(i) + paramIc(j)*sin(ptd.W*j*ptd.SimulationTime(i)+pi/180*paramPhIc(j));           
               end
            end
            ptd.Va = Va;
            ptd.Vb = Vb;
            ptd.Vc = Vc;
            ptd.Ia = Ia;
            ptd.Ib = Ib;
            ptd.Ic = Ic;
            handles.ptd = ptd;
    end  
end

guidata(hObject, handles);

% --- Executes on button press in pbClear.
function pbClear_Callback(hObject, eventdata, handles)
% hObject    handle to pbClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edFrequency_Callback(hObject, eventdata, handles)
% hObject    handle to edFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edFrequency as text
%        str2double(get(hObject,'String')) returns contents of edFrequency as a double


% --- Executes during object creation, after setting all properties.
function edFrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edFrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edVoltagePhase_Callback(hObject, eventdata, handles)
% hObject    handle to edVoltagePhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edVoltagePhase as text
%        str2double(get(hObject,'String')) returns contents of edVoltagePhase as a double


% --- Executes during object creation, after setting all properties.
function edVoltagePhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edVoltagePhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edCurrentPhase_Callback(hObject, eventdata, handles)
% hObject    handle to edCurrentPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCurrentPhase as text
%        str2double(get(hObject,'String')) returns contents of edCurrentPhase as a double


% --- Executes during object creation, after setting all properties.
function edCurrentPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCurrentPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Read_from_file_Callback(hObject, eventdata, handles)
% hObject    handle to Read_from_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Import the file
% Import the file
ptd = handles.ptd;
[file,path] = uigetfile('*.mat');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
data = load( file);
% Create new variables in the base workspace from those fields.
vars = fieldnames(data);
for i = 1:length(vars)
    assignin('base', vars{i}, data.(vars{i}));
end
Matrix = data.(vars{1});
[m,n] = size(Matrix);
if m == 3
    ptd.LoadedDataSinglePhase = true;
    ptd.Va = Matrix(1,:);
    ptd.Ia = Matrix(2,:);
    ptd.SimulationTime = Matrix(3,:);
    ptd.Fs = length(Matrix(3,:))/(Matrix(3,end)-Matrix(3,1));
     
else
  
end
handles.ptd = ptd;
guidata(hObject, handles); 


% --- Executes on button press in tbGrid.
function tbGrid_Callback(hObject, eventdata, handles)
% hObject    handle to tbGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB


% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbGrid

if ( get(hObject,'Value') == true)
    grid( handles.axes1, 'on' )
    grid( handles.axes2, 'on' )
else
    grid( handles.axes1, 'off' )
    grid( handles.axes2, 'off' )
end
