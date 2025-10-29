clc
clear
close all;

%% -------------------------QPSK & 16-QAM with No coding---------------------------%%

Number_of_Bits = 435200;     % 43520 , 435200 , 1044480
Random_Data = randi([0 1],1,Number_of_Bits);
Eb = 1;
Eb_QAM = 2.5 * Eb;
SNR_range_db = (-5:15);   

%---------------------Interleaver------------------%
for i=1:256:Number_of_Bits 
     QPSK_intrlv(1,i:i+255) = matintrlv(Random_Data(1,i:i+255),16,16);      % for QPSK Modulation
end

for i=1:512:Number_of_Bits 
     QAM_intrlv(1,i:i+511) = matintrlv(Random_Data(1,i:i+511),32,16);       % for 16-QAM Modulation
end

%------------------------Mapper----------------------%
QPSK_mapped = Mapper(QPSK_intrlv,2,Eb);            
QAM_mapped = Mapper(QAM_intrlv,4,Eb_QAM);

%-------------------------IFFT-----------------------%
QPSK_IFFT = zeros(1,length(QPSK_mapped));
for i=1:128:length(QPSK_mapped)
    QPSK_IFFT(1,i:i+127) = ifft(QPSK_mapped(1,i:i+127),128);
end

QAM_IFFT = zeros(1,length(QAM_mapped));
for i=1:128:length(QAM_mapped)
    QAM_IFFT(1,i:i+127) = ifft(QAM_mapped(1,i:i+127),128);
end

%---------------------Cyclic Extension-----------------%
count1 = 1; 
QPSK_CYC = [];              % Initialize the QPSK_CYC matrix
for i = 1:128:(Number_of_Bits/2) 
    % Add cyclic prefix with size 32 (cyclic prefix length) and append to QPSK_IFFT to produce 160 bit block QPSK_cyclic(128+32)
    QPSK_CYC((160 * (count1 - 1) + 1):(count1 * 160)) = [QPSK_IFFT(1, ((i + 96):(i + 127))), QPSK_IFFT(1, (i:(i + 127)))];
    count1 = count1 + 1; 
end

count2 = 1; 
QAM_CYC = [];              % Initialize the QAM_CYC matrix
for i = 1:128:(Number_of_Bits/4) 
    % Add cyclic prefix and append to QPSK_cyclic
    QAM_CYC((160 * (count2 - 1) + 1):(count2 * 160)) = [QAM_IFFT(1, ((i + 96):(i + 127))), QAM_IFFT(1, (i:(i + 127)))];
    count2 = count2 + 1; 
end

%-------------------------Channel---------------------%

% outputs from Rayleigh flat fading channel
QPSK_ch = Flat_channel(QPSK_CYC,SNR_range_db);
QAM_ch = Flat_channel(QAM_CYC,SNR_range_db);

% outputs from Frequency selective (FS) Fading channel
QPSK_ch_FS = FS_channel(QPSK_CYC,SNR_range_db);
QAM_ch_FS = FS_channel(QAM_CYC,SNR_range_db);

%=======================At Reciever=========================%

%--------------------De-Cyclic Extension--------------------%
QPSK_DeCYC = Decyclic(QPSK_ch,SNR_range_db);           % First case for QPSK [flat channel]
QPSK_DeCYC_FS = Decyclic(QPSK_ch_FS,SNR_range_db);     % Second case for QPSK [frequency selective channel]
QAM_DeCYC = Decyclic(QAM_ch,SNR_range_db);             % Third case for QAM [flat channel]
QAM_DeCYC_FS = Decyclic(QAM_ch_FS,SNR_range_db);       % Fourth case for QAM [frequency selective channel]

%--------------------------FFT---------------------------%
QPSK_FFT = FFT(QPSK_DeCYC,SNR_range_db);
QPSK_FFT_FS = FFT(QPSK_DeCYC_FS,SNR_range_db);
QAM_FFT = FFT(QAM_DeCYC,SNR_range_db);
QAM_FFT_FS = FFT(QAM_DeCYC_FS,SNR_range_db);

%-----------------------Demapper-------------------------%
QPSK_demapped = DeMapper(QPSK_FFT,SNR_range_db,2,Eb);
QPSK_demapped_FS = DeMapper(QPSK_FFT_FS,SNR_range_db,2,Eb);
QAM_demapped = DeMapper(QAM_FFT,SNR_range_db,4,Eb_QAM);
QAM_demapped_FS = DeMapper(QAM_FFT_FS,SNR_range_db,4,Eb_QAM);

%--------------------Deinterleaver-----------------------%
QPSK_Deintrlv = Deinterleaver(QPSK_demapped,SNR_range_db,2);
QPSK_FS_Deintrlv = Deinterleaver(QPSK_demapped_FS,SNR_range_db,2);
QAM_Deintrlv = Deinterleaver(QAM_demapped,SNR_range_db,4);
QAM_FS_Deintrlv = Deinterleaver(QAM_demapped_FS,SNR_range_db,4);

%-------------------BER Calculation-----------------------%
BER_QPSK = compute_BER(QPSK_Deintrlv , Random_Data,SNR_range_db);
BER_QPSK_FS = compute_BER(QPSK_FS_Deintrlv , Random_Data,SNR_range_db);
BER_QAM = compute_BER(QAM_Deintrlv , Random_Data,SNR_range_db);
BER_QAM_FS = compute_BER(QAM_FS_Deintrlv , Random_Data,SNR_range_db);


%% -------------------------QPSK & 16-QAM with Repeatition coding---------------------%%

QPSK_rep = SymbolRepetitionEncoder(Random_Data,2);
QAM_rep = SymbolRepetitionEncoder(Random_Data,4);

%---------------------Interleaver (Repetition)------------------%
for i=1:255:length(QPSK_rep)
     intrlv_out = matintrlv([QPSK_rep(1,i:i+254),0],16,16);      % for QPSK Modulation with repeatiton
     if(i==1)
       QPSK_rep_intrlv = intrlv_out;
     else
       QPSK_rep_intrlv = [QPSK_rep_intrlv intrlv_out];
     end
end

for i=1:510:length(QAM_rep) 
     intrlv_out = matintrlv([QAM_rep(1,i:i+509),0,0],32,16);      % for 16-QAM Modulation with repeatition
     if(i==1)
        QAM_rep_intrlv = intrlv_out;
     else
        QAM_rep_intrlv = [QAM_rep_intrlv intrlv_out];
     end
end

%------------------------Mapper (Repetition)----------------------%
QPSK_rep_mapped = Mapper(QPSK_rep_intrlv,2,Eb/sqrt(3));         % Same energy per info
QAM_rep_mapped = Mapper(QAM_rep_intrlv,4,Eb_QAM/sqrt(3));       % Same energy per info

%-------------------------IFFT (Repetition)-----------------------%
QPSK_rep_IFFT = zeros(1,length(QPSK_rep_mapped));
for i=1:128:length(QPSK_rep_mapped)
    QPSK_rep_IFFT(1,i:i+127) = ifft(QPSK_rep_mapped(1,i:i+127),128);
end

QAM_IFFT = zeros(1,length(QAM_rep_mapped));
for i=1:128:length(QAM_rep_mapped)
    QAM_rep_IFFT(1,i:i+127) = ifft(QAM_rep_mapped(1,i:i+127),128);
end

%---------------------Cyclic Extension (Repetition)-----------------%
count3 = 1; 
QPSK_rep_CYC = [];              % Initialize the QPSK_CYC matrix
for i = 1:128:length(QPSK_rep_mapped) 
    QPSK_rep_CYC((160 * (count3 - 1) + 1):(count3 * 160)) = [QPSK_rep_IFFT(1,((i + 96):(i + 127))), QPSK_rep_IFFT(1,(i:(i + 127)))];
    count3 = count3 + 1; 
end

count4 = 1; 
QAM_rep_CYC = [];
for i = 1:128:length(QAM_rep_mapped)
    QAM_rep_CYC((160 * (count4 - 1) + 1):(count4 * 160)) = [QAM_rep_IFFT(1,((i + 96):(i + 127))), QAM_rep_IFFT(1,(i:(i + 127)))];
    count4 = count4 + 1; 
end

%-------------------------------Channel----------------------------------%
% outputs from Rayleigh flat fading channel
QPSK_rep_ch = Flat_channel(QPSK_rep_CYC,SNR_range_db);    % for QPSK
QAM_rep_ch = Flat_channel(QAM_rep_CYC,SNR_range_db);      % for QAM

% outputs from Frequency selective (FS) Fading channel
QPSK_rep_ch_FS = FS_channel(QPSK_rep_CYC,SNR_range_db);   % for QPSK
QAM_rep_ch_FS = FS_channel(QAM_rep_CYC,SNR_range_db);     % for QAM

%=============================At Reciever=============================%

%--------------------De-Cyclic Extension--------------------%
% First case for QPSK [flat channel]
QPSK_rep_DeCYC = Decyclic(QPSK_rep_ch,SNR_range_db);   
% Second case for QPSK [frequency selective channel]
QPSK_rep_DeCYC_FS = Decyclic(QPSK_rep_ch_FS,SNR_range_db);     
% Third case for QAM [flat channel]  
QAM_rep_DeCYC = Decyclic(QAM_rep_ch,SNR_range_db);             
% Fourth case for QAM [frequency selective channel]
QAM_rep_DeCYC_FS = Decyclic(QAM_rep_ch_FS,SNR_range_db);      

%--------------------------FFT (Repitition)---------------------------%
QPSK_rep_FFT = FFT(QPSK_rep_DeCYC,SNR_range_db);        % QPSK [Flat-Fading]
QPSK_rep_FFT_FS = FFT(QPSK_rep_DeCYC_FS,SNR_range_db);  % QPSK [Frequency Slective]

QAM_rep_FFT = FFT(QAM_rep_DeCYC,SNR_range_db);          % QAM [Flat-Fading]
QAM_rep_FFT_FS = FFT(QAM_rep_DeCYC_FS,SNR_range_db);    % QAM [Frequency Slective]

%-----------------------Demapper (Repitition)-------------------------%
QPSK_rep_demapped = DeMapper(QPSK_rep_FFT,SNR_range_db,2,Eb);              % QPSK [Flat-Fading]
QPSK_rep_demapped_FS = DeMapper(QPSK_rep_FFT_FS,SNR_range_db,2,Eb);        % QPSK [Frequency Slective]

QAM_rep_demapped = DeMapper(QAM_rep_FFT,SNR_range_db,4,Eb_QAM);            % QAM [Flat-Fading]
QAM_rep_demapped_FS = DeMapper(QAM_rep_FFT_FS,SNR_range_db,4,Eb_QAM);      % QAM [Frequency Slective]

%--------------------Deinterleaver(Repitition)-----------------------%
QPSK_rep_Deintrlv = Deinterleaver(QPSK_rep_demapped,SNR_range_db,2);
QPSK_rep_FS_Deintrlv = Deinterleaver(QPSK_rep_demapped_FS,SNR_range_db,2);
QAM_rep_Deintrlv = Deinterleaver(QAM_rep_demapped,SNR_range_db,4);
QAM_rep_FS_Deintrlv = Deinterleaver(QAM_rep_demapped_FS,SNR_range_db,4);

% Remove the padding done at the interleaving stage at the Tx 
QPSK_rep_Deintrlv = pad_removal(QPSK_rep_Deintrlv,2);
QPSK_rep_FS_Deintrlv = pad_removal(QPSK_rep_FS_Deintrlv,2);
QAM_rep_Deintrlv = pad_removal(QAM_rep_Deintrlv,4);
QAM_rep_FS_Deintrlv = pad_removal(QAM_rep_FS_Deintrlv,4);

%------------------------Majority Decoding---------------------------%
for i=1:length(SNR_range_db)  
    QPSK_rep_decoded(i,:) = majority_decoder(QPSK_rep_Deintrlv(i,:),2);
    QPSK_rep_FS_decoded(i,:) = majority_decoder(QPSK_rep_FS_Deintrlv(i,:),2);
    QAM_rep_decoded(i,:) = majority_decoder(QAM_rep_Deintrlv(i,:),4);
    QAM_rep_FS_decoded(i,:) = majority_decoder(QAM_rep_FS_Deintrlv(i,:),4);
end

%-------------------BER Calculation (Repitition)------------------------%
BER_QPSK_rep = compute_BER(QPSK_rep_decoded , Random_Data,SNR_range_db);
BER_QPSK_rep_FS = compute_BER(QPSK_rep_FS_decoded, Random_Data,SNR_range_db);
BER_QAM_rep = compute_BER(QAM_rep_decoded , Random_Data,SNR_range_db);
BER_QAM_rep_FS = compute_BER(QAM_rep_FS_decoded , Random_Data,SNR_range_db);


%--------------------Plotting BER vs SNR(Eb/No)------------------------%
figure;
semilogy(SNR_range_db , BER_QPSK ,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_QPSK_rep ,'r','linewidth',1.5);
xlabel('Eb/No (dB)');
grid on;
title('BER of Flat fading channel QPSK');
legend('without coding','with Repitition coding','location','southwest')

figure;
semilogy(SNR_range_db , BER_QPSK_FS ,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_QPSK_rep_FS ,'r','linewidth',1.5);
xlabel('Eb/No (dB)')
grid on;
title('BER of frequency selective channel QPSK');
legend('without coding','with Repitition coding','location','southwest')

figure;
semilogy(SNR_range_db , BER_QAM ,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_QAM_rep ,'r','linewidth',1.5);
xlabel('Eb/No (dB)');
grid on;
title('BER of Flat fading channel 16-QAM');
legend('without coding','with Repitition coding','location','southwest')

figure;
semilogy(SNR_range_db , BER_QAM_FS ,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_QAM_rep_FS ,'r','linewidth',1.5);
xlabel('Eb/No (dB)');
grid on;
title('BER of frequency selective channel 16-QAM');
legend('without coding','with Repitition coding','location','southwest')


%% ----------------------Repeatition encoder--------------------------- %%
function encodedStream = SymbolRepetitionEncoder(inputStream, symbolBits)   % Repeats each symbol by 3 times
    % Reshape the input stream into symbols
    numSymbols = length(inputStream) / symbolBits;
    symbols = reshape(inputStream, symbolBits, numSymbols)';

    % Repeat each symbol 3 times
    repeatedSymbols = repelem(symbols, 3, 1);

    % Flatten the repeated symbols back into a binary stream
    encodedStream = repeatedSymbols';
    encodedStream = encodedStream(:)';
end

%% ------------------------Mapper function-----------------------------%%
function mapped_constellations = Mapper(inputstream, symbolBits, Eb)            % Eb represents the info energy
Number_of_Bits = length(inputstream);
mapped_constellations = zeros(1,Number_of_Bits/symbolBits);

switch symbolBits
    case 2              % QPSK mapping
        for i = 1:symbolBits:Number_of_Bits
           if(inputstream(i) == 0 && inputstream(i+1) == 0 )
              mapped_constellations((i+1)/2) = sqrt(Eb) * complex(-1,-1);
           elseif(inputstream(i) == 0 && inputstream(i+1) == 1)
              mapped_constellations((i+1)/2) = sqrt(Eb) * complex(-1,1);
           elseif(inputstream(i) == 1 && inputstream(i+1) == 0)
              mapped_constellations((i+1)/2) = sqrt(Eb) * complex(1,-1);
           elseif(inputstream(i) == 1 && inputstream(i+1) == 1)
              mapped_constellations((i+1)/2) = sqrt(Eb) * complex(1,1);
           end
        end
    case 4              % 16-QAM mapping
        real_MQAM = zeros(1,length(inputstream)/symbolBits);
        img_MQAM = zeros(1,length(inputstream)/symbolBits);
        for i = 1:symbolBits:Number_of_Bits
            if(inputstream(i) == 0 && inputstream(i+1) == 0)            % First two bits control the real part of the MQAM signal
                real_MQAM((i+3)/4) = -3;
            elseif(inputstream(i) == 0 && inputstream(i+1) == 1)
                real_MQAM((i+3)/4) = -1;
            elseif(inputstream(i) == 1 && inputstream(i+1) == 1)
                real_MQAM((i+3)/4) = 1;
            elseif(inputstream(i) == 1 && inputstream(i+1) == 0)
                real_MQAM((i+3)/4) = 3;
            end
            if(inputstream(i+2) == 0 && inputstream(i+3) == 0)           % Second two bits control the imaginary part of the MQAM signal
                img_MQAM((i+3)/4) = -3;
            elseif(inputstream(i+2) == 0 && inputstream(i+3) == 1)
                img_MQAM((i+3)/4) = -1;
            elseif(inputstream(i+2) == 1 && inputstream(i+3) == 1)
                img_MQAM((i+3)/4) = 1;
            elseif(inputstream(i+2) == 1 && inputstream(i+3) == 0)
                img_MQAM((i+3)/4) = 3;
            end
        mapped_constellations((i+3)/4) = sqrt(Eb) * complex(real_MQAM((i+3)/4),img_MQAM((i+3)/4));
        end
end
        
end

%% -------------Rayleigh flat fading channel Function------------------%%
function Received_data = Flat_channel(cyclic_input,SNR_range_db)

SNR_range_linear = zeros(1,length(SNR_range_db));
Received_data = zeros(length(SNR_range_db),length(cyclic_input));
    for i = 1:length(SNR_range_db)
        SNR_range_linear(i) = 10 ^ (SNR_range_db(i)/10);
        for j = 1 : 160 : length(cyclic_input)
            Eavg = mean(abs(cyclic_input(j:j+159)).^2);
            h_channel = (1/sqrt(2)) .* complex(randn(1,160) , randn(1,160));
            noise = sqrt(Eavg/(2*SNR_range_linear(i))) .* complex(randn(1,160) , randn(1,160));
            Received_data(i,j:j+159) = (cyclic_input(j:j+159) .* h_channel + noise) ./ h_channel;
        end
    end
end

%% -------------Frequency selective fading channel Function-------------%%
function Received_data = FS_channel(cyclic_input,SNR_range_db)  
Cyclic_symbols = 160;
Number_of_subchannels = 2;
channel_length = Cyclic_symbols/Number_of_subchannels;  
gain = [0.5,0.5];  % FS gain of the subchannels (2 subchannels) (equal channel assignment)
% Each Rayleigh subchannel is scaled with different gain
ch1_time = gain(1)*(1/sqrt(2)) * complex(randn(1,channel_length),randn(1,channel_length));
ch2_time = gain(2)*(1/sqrt(2)) * complex(randn(1,channel_length),randn(1,channel_length));
data_due_channel = zeros(1,(length(cyclic_input)/channel_length)*(2*Cyclic_symbols-2)); 

indx = 1;
SNR_range_linear = zeros(1,length(SNR_range_db));
% Received_data = zeros(length(SNR_range_db),length(cyclic_input));
    for i= 1:length(SNR_range_db)
        SNR_range_linear(i) = 10 ^ (SNR_range_db(i)/10);
        for j = 1:channel_length:length(cyclic_input)
            if(mod(indx,2)==1) % Assign to channel 1
                data_due_channel(j:j+channel_length-1) = cyclic_input(j:j+channel_length-1) .* ch1_time;
                Eavg = mean(abs(cyclic_input(j:j+channel_length-1)).^2);
                noise = sqrt(Eavg/(2*SNR_range_linear(i))) .* complex(randn(1,channel_length) , randn(1,channel_length));
                Received_data(i,j:j+channel_length-1) = (data_due_channel(j:j+channel_length-1) + noise) ./ ch1_time;
            else               % Assign to channel 2 
                data_due_channel(j:j+channel_length-1) = cyclic_input(j:j+channel_length-1) .* ch2_time;
                Eavg = mean(abs(cyclic_input(j:j+channel_length-1)).^2);
                noise = sqrt(Eavg/(2*SNR_range_linear(i))) .* complex(randn(1,channel_length) , randn(1,channel_length));
                Received_data(i,j:j+channel_length-1) = (data_due_channel(j:j+channel_length-1) + noise) ./ ch2_time;
            end
            indx = indx +1;
        end
    end  
end

%% -----------------Decyclic extension Function------------------------%%
function decyclic_out = Decyclic(input_data,SNR_range_db)

%decyclic_out = zeros(length(SNR_range_db),(Number_of_Bits/symbolBits));
% decyclic_out = zeros(length(SNR_range_db),size(input_data,2));
    for j=1:length(SNR_range_db)
        for i=160:160:size(input_data,2)
            decyclic_out(j,((i/160)*128)-127:(i/160)*128)=input_data(j,i-127:i);
        end
    end
end

%% -----------------------FFT Function (At Rx)------------------------%%
function out_data_fft = FFT(input_data,SNR_range_db)

out_data_fft = zeros(length(SNR_range_db), size(input_data,2));        % Initialize output array

    for j = 1 : length(SNR_range_db)
        for i = 1:128:size(input_data,2)
            out_data_fft(j, i:i+127) = fft(input_data(j, i:i+127), 128);
        end
    end
end

%% ------------------ DeMapper Function----------------------%%
function demapped_data_binary = DeMapper (received_data ,SNR_range_db,symbolBits,Eb)

demapped_data_binary = zeros(length(SNR_range_db),size(received_data,2)*symbolBits);
demapped_data_decimal = zeros(length(SNR_range_db),size(received_data,2));

    for i = 1 : length(SNR_range_db)
        switch symbolBits
            case 2       % QPSK case
                QPSK_table = sqrt(Eb) * [complex(-1,-1), complex(-1,1), complex(1,-1), complex(1,1)];
                for j = 1 : size(received_data,2)
                    [~ , Min_index] = min(abs(received_data(i,j) - QPSK_table));
                    demapped_data_decimal(i,j) = Min_index - 1;
                end
                demapped_data_binary(i,:) = reshape(de2bi(demapped_data_decimal(i,:),2,'left-msb')',1,[]);
            case 4       % 16-QAM case
                MQAM_table = sqrt(Eb)*[complex(-3,-3), complex(-3,-1), complex(-3,3), complex(-3,1), complex(-1,-3), complex(-1,-1), complex(-1,3), complex(-1,1), complex(3,-3), complex(3,-1), complex(3,3), complex(3,1),complex(1,-3), complex(1,-1), complex(1,3), complex(1,1)];                     
                for j = 1 : size(received_data,2)
                    [~ , Min_index] = min(abs(received_data(i,j) - MQAM_table));
                    demapped_data_decimal(i,j) = Min_index - 1;
                end
                demapped_data_binary(i,:) = reshape(de2bi(demapped_data_decimal(i,:),4,'left-msb')',1,[]);        
        end
    end
end

%% -----------------------Deinterleaver Function---------------------%%
function deintrlv_out = Deinterleaver(demapped_data,SNR_range_db,symbolBits)

deintrlv_out = zeros(length(SNR_range_db), size(demapped_data,2));       % Initialize the output

    for j = 1 : length(SNR_range_db)
        for i = 1:(128*symbolBits):size(demapped_data,2)
            deintrlv_out(j, i:i+(128*symbolBits)-1) = matdeintrlv(demapped_data(j, i:i+(128*symbolBits)-1), symbolBits*8, 16); 
        end
    end
end

%% ------------------------Pad Removal function (In case of Repition code)--------------------%%
function Deintrlv_final = pad_removal(Deintrlv_output,symbolBits)

    for i = 1:size(Deintrlv_output,1)
        k=1;
        if(symbolBits == 2) % QPSK ---> every 256 bits remove 1 bit from the end
            for j  = 1:256:size(Deintrlv_output,2)
                frame = Deintrlv_output(i,j:j+255);
                Deintrlv_final(i,k:k+254) = frame(1:end-1);
                k=k+255;
            end
        else % 16-QAM ---> every 512 bits remove 2 bits from the end
            for j  = 1:512:size(Deintrlv_output,2)
                frame = Deintrlv_output(i,j:j+511);
                Deintrlv_final(i,k:k+509) = frame(1:end-2);
                k=k+510;
            end
        end
    end
end

%% ---------------------- Decoding Function ------------------ %%
function [decoded_bitstream] = majority_decoder (encoded_bitstream , symbolBits )
    hash_table = zeros(1,2^symbolBits);
    iterations = 0;
    j=1;
    decoded_bitstream = zeros(1,length(encoded_bitstream)/3); %% Assume 3-repeatition code 
    for i = 1 : symbolBits : length(encoded_bitstream)
        iterations = iterations + 1 ;
        % Convert the symbol into the correspoding decimal into a hashtable
        % to see the no. of occurrences of such symbol for decision
        index = bin2dec(num2str((encoded_bitstream(i:i+symbolBits-1))))+1;
        hash_table(index) = hash_table(index) + 1;
        if(iterations==3)
             %decide based on the most no. of occurrences
             % then get the corrsponding binary symbol
             index_of_most_repeated_symbol=find(hash_table==max(hash_table),1);
             binaryStr = dec2bin(index_of_most_repeated_symbol-1, symbolBits);
             decoded_bitstream(j:j+symbolBits-1) = double(binaryStr)-'0';
             j = j + symbolBits;
             hash_table = zeros(1,2^symbolBits); % Reset the hash table for the next set of symbols
             iterations = 0;
        end
    end
end

%% ------------------------BER Calculation function--------------------%%
function BER_value = compute_BER(received_data,Random_data,SNR_range_db)

BER_value = zeros(1,length(SNR_range_db));
    for i = 1 : length(SNR_range_db)
        BER_value(1,i) = length(find(xor(received_data(i,:),Random_data) == 1))/length(Random_data);
    end
end
