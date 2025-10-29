clc
clear
close all;

%% ------------------------------ Generating the Random Data Sequence -----------------------------------%%

% we can choose any value of number of bits under conditions that it can be divided into number of bits in each symbol that is 1(for BPSK) & 2(for QPSK) & 4(for 16-QAM)
Number_of_Bits = 120000;
Eb = 1;                   % Eb represent the bit eneregy
Eb_QAM = 2.5;
SNR_range_db = (-5:15);   % where SNR_range = Eb/No
SNR_range_linear = 10 .^ (SNR_range_db/10);
No = Eb ./ SNR_range_linear;
Random_Data = randi([0 1] , 1 , Number_of_Bits);   % generate all the random bits
 
%% ------------------------------BPSK Modulation Scheme ---------------------%%
 
% ----------------Declare BPSK (Nocode)--------------%
BPSK = Mapper(Random_Data,1,Eb); % Symbol_Bits = 1 , Eb = 1 

% -----------BPSK Channel stage (No code)-----------%
BPSK_recieved = zeros(1,length(BPSK));
BER_BPSK = zeros(1,length(SNR_range_db));
h_BPSK = (1/sqrt(2)) .* complex(randn(1,length(BPSK)) , randn(1,length(BPSK)));  % channel model h(t) with zero mean & variance = 1/2
 
for i = 1 : length(SNR_range_db) 
    noise_BPSK = sqrt(No(i)/2).* complex(randn(1,length(BPSK)) , randn(1,length(BPSK))) ;   % sqrt(No/2) refer to standard deviation 
    BPSK_recieved = (BPSK .* h_BPSK + noise_BPSK)./ h_BPSK;        % we divide received signal over the channel for applying channel inversion
    BPSK_demapped = DeMapper(BPSK_recieved,1);         % calling the demapper function
    
    %---------BER calculation for the BPSK Scheme(No code)---------%
    error_bits_BPSK = 0;
    for j = 1 : length(Random_Data)
        if(BPSK_demapped(j) ~= Random_Data(j))
            error_bits_BPSK = error_bits_BPSK + 1;
        end
    end
    BER_BPSK(1,i) = error_bits_BPSK/Number_of_Bits;
end    
% ----------------Declare BPSK (Rep Code)--------------%
BPSK_rep_bstream = SymbolRepetitionEncoder(Random_Data,1); 
BPSK_rep = Mapper(BPSK_rep_bstream,1,Eb);

% -----------BPSK Channel stage (Rep code)-----------%
BPSK_rep_recieved = zeros(1,length(BPSK_rep));
BER_BPSK_rep = zeros(1,length(SNR_range_db));

h_BPSK_rep = (1/sqrt(2)) .* complex(randn(1,length(BPSK_rep)) , randn(1,length(BPSK_rep)));  % channel model h(t) with zero mean & variance = 1/2
 
for i = 1 : length(SNR_range_db) 
    noise_BPSK_rep = sqrt(No(i)/2).* complex(randn(1,length(BPSK_rep)) , randn(1,length(BPSK_rep))) ;   % sqrt(No/2) refer to standard deviation 
    BPSK_rep_recieved = (BPSK_rep .* h_BPSK_rep + noise_BPSK_rep)./ h_BPSK_rep;        % we divide received signal over the channel for applying channel inversion
    BPSK_rep_demapped = DeMapper(BPSK_rep_recieved,1);         % calling the demapper function

BPSK_rep_decoded = majority_decoder(BPSK_rep_demapped,1); % Apply the majority decoder function
%---------BER calculation for the BPSK Scheme(Rep code)---------%
    error_bits_BPSK = 0;
    for j = 1 : length(Random_Data)
        if(BPSK_rep_decoded(j) ~= Random_Data(j))
            error_bits_BPSK = error_bits_BPSK + 1;
        end
    end
    BER_BPSK_rep(1,i) = error_bits_BPSK/Number_of_Bits;
end    
  
% ------------------plotting the BER for BPSK--------------%
figure;
semilogy(SNR_range_db , BER_BPSK,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_BPSK_rep,'r','linewidth',1.5);
%hold off;
title('BER of BPSK Modulation Scheme');
xlabel('Eb/No (dB)');
grid on;
legend('BER without repetition coding','BER with repetition coding','Location','southwest')

%% ------------------------------QPSK Modulation Scheme ---------------------%%

% -----------------Declare QPSK (No code)-------------%
QPSK = Mapper(Random_Data,2,Eb);

% ------------------QPSK channel stage(No code)-----------------%
QPSK_recieved = zeros(1,length(QPSK));
BER_QPSK = zeros(1,length(SNR_range_db));
h_QPSK = (1/sqrt(2)) .* complex(randn(1,length(QPSK)) , randn(1,length(QPSK)));

for i = 1 : length(SNR_range_db)   
    noise_QPSK = sqrt(No(i)/2).* complex(randn(1,length(QPSK)) , randn(1,length(QPSK)));   
    QPSK_recieved = (QPSK .* h_QPSK + noise_QPSK) ./ h_QPSK;
    QPSK_demapped = DeMapper(QPSK_recieved,2);      % Calling demapper function
    
    %-----------------------BER calculation for QPSK (No code)--------------------%
    error_bits_QPSK = 0;
    inc_var = 1;                 % this variable specified for incrementing the index of Random_Data variable in the loop
    for j = 1 : length(Random_Data)/2
        for k = 1 : 2
            if(QPSK_demapped(j,k) ~= Random_Data(inc_var))
                error_bits_QPSK = error_bits_QPSK + 1;     
            end
            inc_var = inc_var + 1;
        end
    end  
    BER_QPSK(1,i) = error_bits_QPSK/Number_of_Bits;
end

% -----------------Declare QPSK (Rep code)-------------%
QPSK_rep_bstream =  SymbolRepetitionEncoder(Random_Data,2);
QPSK_rep = Mapper(QPSK_rep_bstream,2,Eb);

% ------------------QPSK channel stage(Rep code)-----------------%
QPSK_rep_recieved = zeros(1,length(QPSK_rep));
BER_QPSK_rep = zeros(1,length(SNR_range_db));
h_QPSK_rep = (1/sqrt(2)) .* complex(randn(1,length(QPSK_rep)) , randn(1,length(QPSK_rep)));

for i = 1 : length(SNR_range_db)   
    noise_QPSK_rep = sqrt(No(i)/2).* complex(randn(1,length(QPSK_rep)) , randn(1,length(QPSK_rep)));   
    QPSK_rep_recieved = (QPSK_rep .* h_QPSK_rep + noise_QPSK_rep) ./ h_QPSK_rep;
    QPSK_demapped_rep = DeMapper(QPSK_rep_recieved,2);      % Calling demapper function
    QPSK_decoded_rep = majority_decoder(reshape(QPSK_demapped_rep',1,[]),2);
    numSymbols = length(QPSK_decoded_rep) / 2;
    QPSK_decoded_rep = inverseReshape(QPSK_decoded_rep,[numSymbols,2]);
    
    %-------------BER calculation for QPSK (Rep code)-------------%
    error_bits_QPSK = 0;
    inc_var = 1;                 
    for j = 1 : length(Random_Data)/2
        for k = 1 : 2
            if(QPSK_decoded_rep(j,k) ~= Random_Data(inc_var))
                error_bits_QPSK = error_bits_QPSK + 1;     
            end
            inc_var = inc_var + 1;
        end
    end  
    BER_QPSK_rep(1,i) = error_bits_QPSK/Number_of_Bits;
end

% -----------------plotting the BER for QPSK------------------%
figure;
semilogy(SNR_range_db , BER_QPSK ,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_QPSK_rep ,'r','linewidth',1.5);
title('BER of QPSK Modulation Scheme');
xlabel('Eb/No (dB)');
grid on;
legend('BER without repetition coding','BER with repetition coding','Location','southwest')

%% ------------------------------16-QAM Modulation Scheme ---------------------%%

% ---------------------Declare 16-QAM(No code)-----------------------%
MQAM = Mapper(Random_Data,4,Eb_QAM);

%---------------16-QAM channel stage(No code)------------%
MQAM_recieved = zeros(1,length(MQAM));
BER_MQAM = zeros(1,length(SNR_range_db));
h_MQAM = (1/sqrt(2)) .* complex(randn(1,length(MQAM)) , randn(1,length(MQAM)));

for i = 1 : length(SNR_range_db)  
    noise_MQAM = sqrt(No(i)/2).* complex(randn(1,length(MQAM)) , randn(1,length(MQAM)));      
    MQAM_recieved = (MQAM .* h_MQAM + noise_MQAM) ./ h_MQAM;
    MQAM_demapped = DeMapper(MQAM_recieved,4);             % Calling demapper function            
    
    %--------------------BER calculation for 16-QAM scheme(No code)----------------%
    error_bits_MQAM = 0;
    inc_var = 1;                 
    for j = 1 : length(Random_Data)/4
        for k = 1 : 4
            if(MQAM_demapped(j,k) ~= Random_Data(inc_var))
                error_bits_MQAM = error_bits_MQAM + 1;
            end
            inc_var = inc_var + 1;
        end
    end  
    BER_MQAM(1,i) = error_bits_MQAM/Number_of_Bits;       
end

% ---------------------Declare 16-QAM(Rep Code)-----------------------%
MQAM_rep_bstream = SymbolRepetitionEncoder(Random_Data,4);
MQAM_rep = Mapper(MQAM_rep_bstream,4,Eb_QAM);

%---------------16-QAM channel stage(Rep code)------------%
MQAM_rep_recieved = zeros(1,length(MQAM_rep));
BER_MQAM_rep = zeros(1,length(SNR_range_db));
h_MQAM_rep = (1/sqrt(2)) .* complex(randn(1,length(MQAM_rep)) , randn(1,length(MQAM_rep)));

for i = 1 : length(SNR_range_db)  
    noise_MQAM_rep = sqrt(No(i)/2).* complex(randn(1,length(MQAM_rep)) , randn(1,length(MQAM_rep)));     
    MQAM_rep_recieved = (MQAM_rep .* h_MQAM_rep + noise_MQAM_rep) ./ h_MQAM_rep;
    MQAM_demapped_rep = DeMapper(MQAM_rep_recieved,4);  % Calling demapper function            
    MQAM_decoded_rep = majority_decoder(reshape(MQAM_demapped_rep',1,[]),4);
    numSymbols = length(MQAM_decoded_rep) / 4;
    MQAM_decoded_rep = inverseReshape(MQAM_decoded_rep,[numSymbols,4]);
    
    %--------------------BER calculation for 16-QAM scheme(Rep code)----------------%
    error_bits_MQAM = 0;
    inc_var = 1;                 
    for j = 1 : length(Random_Data)/4
        for k = 1 : 4
            if(MQAM_decoded_rep(j,k) ~= Random_Data(inc_var))
                error_bits_MQAM = error_bits_MQAM + 1;
            end
            inc_var = inc_var + 1;
        end
    end  
    BER_MQAM_rep(1,i) = error_bits_MQAM/Number_of_Bits;       
end

% -----------------plotting the BER for 16-QAM---------------%
figure;
semilogy(SNR_range_db , BER_MQAM,'k','linewidth',1.5);
hold on;
semilogy(SNR_range_db , BER_MQAM_rep,'r','linewidth',1.5);
title('BER of 16-QAM Modulation Scheme');
xlabel('Eb/No (dB)');
grid on;
legend('BER without repetition coding','BER with repetition coding','Location','southwest')


%% ------------------------Mapper function-----------------------------%%
function mapped_constellations = Mapper(inputstream, symbolBits,Eb)         % Eb represents the info energy
Number_of_Bits = length(inputstream);
mapped_constellations = zeros(1,Number_of_Bits/symbolBits);

switch symbolBits
    case 1 % BPSK mapping
        mapped_constellations = sqrt(Eb) * (2 * inputstream - 1);
    case 2 % QPSK mapping
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
    case 4 % 16-QAM mapping
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

%% ----------------------Repeatition encoder--------------------------- %%
function encodedStream = SymbolRepetitionEncoder(inputStream, symbolBits)       % Repeats each symbol by 3 times
    % Reshape the input stream into symbols
    numSymbols = length(inputStream) / symbolBits;
    symbols = reshape(inputStream, symbolBits, numSymbols)';
    repeatedSymbols = repelem(symbols, 3, 1);
    % Flatten the repeated symbols back into a binary stream
    encodedStream = repeatedSymbols';
    encodedStream = encodedStream(:)';
end

%% ------------------ DeMapper Function----------------------%%
function demapped_data = DeMapper (received_data ,symbolBits)
    Eb = 1; 
    Eb_QAM = 1;  
    demapped_data = zeros(1,length(received_data));
    if symbolBits == 1           % BPSK case
        BPSK_table = sqrt(Eb) * [complex(-1,0), complex(1,0)];
        for j = 1 : length(received_data)
            [~ , Min_index] = min(abs(received_data(j) - BPSK_table));
            demapped_data(j) = Min_index - 1;
        end
        demapped_data = de2bi(demapped_data,1,'left-msb');
    elseif symbolBits == 2       % QPSK case
        QPSK_table = sqrt(Eb) * [complex(-1,-1), complex(-1,1), complex(1,-1), complex(1,1)];
        for j = 1 : length(received_data)
            [~ , Min_index] = min(abs(received_data(j) - QPSK_table));
            demapped_data(j) = Min_index - 1;
        end
        demapped_data = de2bi(demapped_data,2,'left-msb');
    elseif symbolBits == 4       % 16-QAM case
        MQAM_table = sqrt(Eb_QAM)*[complex(-3,-3), complex(-3,-1), complex(-3,3), complex(-3,1), complex(-1,-3), complex(-1,-1), complex(-1,3), complex(-1,1), complex(3,-3), complex(3,-1), complex(3,3), complex(3,1),complex(1,-3), complex(1,-1), complex(1,3), complex(1,1)];                     
        for j = 1 : length(received_data)
            [~ , Min_index] = min(abs(received_data(j) - MQAM_table));
            demapped_data(j) = Min_index - 1;
        end
        demapped_data = de2bi(demapped_data,4,'left-msb');        
    end
end
%% ---------------------- Decoding Function ------------------ %%
function decoded_bitstream = majority_decoder (encoded_bitstream , symbolBits )
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

%% ---------------------- Inverse reshape Function ------------------ %%
function originalMatrix = inverseReshape(rowVector, originalSize)
    % Reshape the row vector into the transposed matrix shape
    transposedMatrix = reshape(rowVector, fliplr(originalSize));
    
    % Transpose the matrix to restore the original order
    originalMatrix = transposedMatrix.';
end
