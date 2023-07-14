clc;
clear all;
% close all;

dim_vec = [10,5,1]; % dimensions of subcodebooks 
n_dim = length(dim_vec)-1; % number of trellis stages
bit_vec = 6*ones(1,n_dim); % bits per subcodebook - Matlab code currently works only for equal number of bits per stream
CB_size_vec = 2.^(bit_vec); % subcodebook sizes
NN_rand = 1; % number of random codebook realizations simulated
BB = 500; % number of codeblocks per random codebook realization
KK = 2^(10); % number of info bits per block
% KK = 1120;
NN = 3*KK;  % number of coded bits per block (lets use rate 1/3; i.e., the base rate of the turbo code -- no rate matching required)
n_sym = NN/bit_vec(1);
n = dim_vec(1);
m = dim_vec(end);
SNR_vec_dB = [2.75:0.25:6]; % SNRs to simulate
% SNR_vec_dB = 6;
sigma_n2_vec = 10.^(-SNR_vec_dB/10); % noise variance
% sigma_n2_vec = sigma_n2_vec(end-1);
Nr = 1; % number of receive antennas
pruning_percentage = ones(1,n_dim); % trellis pruning (not used)
activate_interleaving = true; % interleaving of bits
NN_perm = 1e5; % number of random permutations considered in bitmapping
opt_codebook = true;
CB_best = {};
dim_string = [''];
for d_i = 1:length(dim_vec)
    dim_string = [dim_string num2str(dim_vec(d_i)) '-'];
end
dim_string = dim_string(1:end-1);
if opt_codebook    
    BB = BB*NN_rand;
    NN_rand = 1;    
    file_name = ['NC_CB' dim_string '_CB' num2str(sum(log2(CB_size_vec))) '_new.mat'];
%     file_name = ['NC_CB' num2str(dim_vec(1)) '-' num2str(dim_vec(1)-dim_vec(2)) '-' num2str(dim_vec(end)) '_CB' num2str(sum(log2(CB_size_vec))) '_new.mat'];    
    load(file_name,'CB_best');
end

disp('Notice: SIC is still based on hypergeometric LLRs --> only first stream is useful for LLR comparison')


if dim_vec(end) > 1
    error('Multidimensional case not implemented')
end

% adapt the MGFs below if dim_vec is changed
for d_i = 1:n_dim  
    alpha_b = dim_vec(end); % parameters of hypergeom
    beta_b = dim_vec(d_i+1)-dim_vec(end); % parameters of hypergeom

    syms X 
    MGF = hypergeom(alpha_b,alpha_b+beta_b,X);
     MGF    % evaluate this if you change dim_vec and update mgf below
    
end
if n_dim == 1
    mgf{1} = @(x) exp(x);
elseif n_dim == 2
    if prod(dim_vec == [10,5,1]) || prod(dim_vec == [12,5,1])
        mgf{1} = @(x)(-(24*x - 24*exp(x) + 12*x.^2 + 4*x.^3 + 24)./x.^4); % this is for dim_vec = [10,5,1]
    elseif prod(dim_vec == [10,4,1]) || prod(dim_vec == [12,4,1])
        mgf{1} = @(x)(-(6*x - 6*exp(x) + 3*x.^2 + 6)./x.^3); % this is for dim_vec = [10,4,1]
    elseif prod(dim_vec == [10,3,1]) || prod(dim_vec == [12,3,1])
        mgf{1} = @(x)(-(2*x - 2*exp(x) + 2)./x.^2); % this is for dim_vec = [10,3,1]
    elseif prod(dim_vec == [12,2,1])
        mgf{1} = @(x)(exp(x) - 1)./x;
    elseif prod(dim_vec == [10,6,1])
        mgf{1} = @(x)(-(120*x - 120*exp(x) + 60*x.^2 + 20*x.^3 + 5*x.^4 + 120)./x.^5);
    else
        error('Specify MGF')
    end
elseif n_dim == 3
    if prod(dim_vec == [10,5,2,1])
        mgf{1} = @(x)(-(24*x - 24*exp(x) + 12*x.^2 + 4*x.^3 + 24)./x.^4);
        mgf{2} = @(x)(exp(x) - 1)./x;
    elseif prod(dim_vec == [10,6,3,1]) || prod(dim_vec == [12,6,3,1])
        mgf{1} = @(x)(-(120*x - 120*exp(x) + 60*x.^2 + 20*x.^3 + 5*x.^4 + 120)./x.^5);
        mgf{2} = @(x)(-(2*x - 2*exp(x) + 2)./x.^2);
    elseif prod(dim_vec == [12,5,3,1])
        mgf{1} = @(x)(-(24*x - 24*exp(x) + 12*x.^2 + 4*x.^3 + 24)./x.^4);
        mgf{2} = @(x)(-(2*x - 2*exp(x) + 2)./x.^2);
    elseif prod(dim_vec == [12,7,3,1])
        mgf{1} = @(x)(-(720*x - 720*exp(x) + 360*x.^2 + 120*x.^3 + 30*x.^4 + 6*x.^5 + 720)./x.^6);
        mgf{2} = @(x)(-(2*x - 2*exp(x) + 2)./x.^2);
    else
        error('Specify MGF')
    end
elseif n_dim == 4
    if  prod(dim_vec == [12,8,5,3,1])
        mgf{1} = @(x)( -(5040*x - 5040*exp(x) + 2520*x.^2 + 840*x.^3 + 210*x.^4 + 42*x.^5 + 7*x.^6 + 5040)./x.^7);
        mgf{2} = @(x)(-(24*x - 24*exp(x) + 12*x.^2 + 4*x.^3 + 24)./x.^4);
        mgf{3} = @(x)(-(2*x - 2*exp(x) + 2)./x.^2);
    elseif prod(dim_vec == [12,8,4,2,1])
        mgf{1} = @(x)( -(5040*x - 5040*exp(x) + 2520*x.^2 + 840*x.^3 + 210*x.^4 + 42*x.^5 + 7*x.^6 + 5040)./x.^7);
        mgf{2} = @(x)(-(6*x - 6*exp(x) + 3*x.^2 + 6)./x.^3); 
        mgf{3} = @(x)(exp(x) - 1)./x;
    elseif prod(dim_vec == [12,7,4,2,1])
        mgf{1} = @(x)(-(720*x - 720*exp(x) + 360*x.^2 + 120*x.^3 + 30*x.^4 + 6*x.^5 + 720)./x.^6);
        mgf{2} = @(x)(-(6*x - 6*exp(x) + 3*x.^2 + 6)./x.^3); 
        mgf{3} = @(x)(exp(x) - 1)./x;
    else
        error('Specify MGF')
    end
end
mgf{n_dim} = @(x) exp(x); % for the last stage there is nothing to take an expectation over

bit_2_sym = cell(n_dim,1); 
same_bit_sym = cell(n_dim,1);
hamming_profile = cell(n_dim,1);
for d_i = 1:n_dim
    bit_2_sym{d_i} = (2.^(0:bit_vec(d_i)-1)).'; % used to convert from bits to symbol index (simple conversion from binary number to symbol index --> bit-mapping not optimized)
    bit_map = de2bi(0:CB_size_vec(d_i)-1); % all possible bit sequences    
    for b_i = 1:bit_vec(d_i)
        same_bit_sym{d_i}(b_i,:) = find(bit_map(:,b_i) == 1); % indices of all symbols for which b_i == 1
        for b_ii = 1:bit_vec(d_i)
            
        end
    end
    temp_profile = zeros(CB_size_vec(d_i),CB_size_vec(d_i));
    for b_i = 1:CB_size_vec(d_i)
        for b_ii = 1:CB_size_vec(d_i)
            temp_profile(b_i,b_ii) = sum(abs(bit_map(b_i,:) - bit_map(b_ii,:)));
        end
    end
    hamming_profile{d_i} = temp_profile/max(temp_profile(:));
end

% SER = zeros(length(sigma_n2_vec),n_dim); % symbol error ratio of trellis decoder
% SER_coded = zeros(length(sigma_n2_vec),n_dim); % symbol error ratio of trellis decoder
% BER = zeros(length(sigma_n2_vec),n_dim); % bit error ratio using hypergeometric function
% BER_uncoded = zeros(length(sigma_n2_vec),n_dim); % bit error ratio without coding (systematic bits)
SER_coded = true(n_sym,n_dim,BB,NN_rand,length(sigma_n2_vec)); % symbol error ratio of trellis decoder
BER = true(KK,n_dim,BB,NN_rand,length(sigma_n2_vec)); % bit error ratio using hypergeometric function
BER_jensen = true(KK,n_dim,BB,NN_rand,length(sigma_n2_vec)); % bit error ratio using hypergeometric function
BER_max = true(KK,n_dim,BB,NN_rand,length(sigma_n2_vec)); % bit error ratio using hypergeometric function
BER_uncoded = true(KK,n_dim,BB,NN_rand,length(sigma_n2_vec)); % bit error ratio without coding (systematic bits)
for snr_i = 1:length(sigma_n2_vec) % loop over SNR - can use parfor here
    snr_i
    seed = RandStream('mt19937ar','Seed',1); % random number seed (required for RANDOM_MIMO_CB) -- easier to debug since we can reproduce results
    sigma_n2 = sigma_n2_vec(snr_i); % noise variance
    rho = 1/sigma_n2; % SNR
    rho_term = (rho*n)/(1+rho*n); % SNR-term in the exponent; notice, there is actually a 1/sigma_2-term missing in our case, since we do not use normalized noise
%     sym_error_coded_sum = zeros(n_dim,1); % counter of symbol errors
%     bit_error_sum = zeros(n_dim,1); % counters of bit errors
%     bit_error_uncoded_sum = zeros(n_dim,1);

    sym_errors = true(n_sym,n_dim,BB,NN_rand);
    bit_errors = true(KK,n_dim,BB,NN_rand);
    bit_errors_jensen = true(KK,n_dim,BB,NN_rand);
    bit_errors_max = true(KK,n_dim,BB,NN_rand);
    bit_errors_uncoded = true(KK,n_dim,BB,NN_rand);

    for nn = 1:NN_rand % loop over random codebook realizations 
        if ~mod(nn-1,1)
            nn
        end
        CB = cell(n_dim,1);
        for d_i = 1:n_dim
            if opt_codebook
                CB{d_i,1} = CB_best{d_i,1};
            else
                CB{d_i,1} = RANDOM_MIMO_CB(dim_vec(d_i+1),dim_vec(d_i),CB_size_vec(d_i),seed,0,1,1); % random subcodebooks
            end
            dist_profile = zeros(CB_size_vec(d_i),CB_size_vec(d_i));
            for cb_i = 1:CB_size_vec(d_i)
                for cb_ii = 1:CB_size_vec(d_i)
                    dist_profile(cb_i,cb_ii) = abs(1 - 1/dim_vec(d_i+1)*norm(CB{d_i,1}(:,:,cb_i)'*CB{d_i,1}(:,:,cb_ii),'fro')^2);
                end
            end
            dist_ratio = hamming_profile{d_i}./dist_profile;
            min_max_ratio = max(dist_ratio(:));            
            best_perm = 1:CB_size_vec(d_i);
            for nn_r = 1:NN_perm
                rand_perm = randperm(seed,CB_size_vec(d_i));
                dist_prof_temp = dist_profile(:,rand_perm);
                dist_prof_temp = dist_prof_temp(rand_perm,:);
                dist_ratio = hamming_profile{d_i}./dist_prof_temp;
                temp_max_ratio = max(dist_ratio(:));
                if temp_max_ratio < min_max_ratio
                    min_max_ratio = temp_max_ratio;
                    best_perm = rand_perm;
                end
            end
            CB{d_i,1} = CB{d_i,1}(:,:,best_perm);
        end
        for bb = 1:BB % loop over transport blocks
            if ~mod(bb-1,5)
                bb
            end
            tx_bits = cell(n_dim,1); % information bits
            tx_coded_bits = cell(n_dim,1); % coded bits
            U_tx_sub = cell(n_dim,1); % selected subcodebook constellation point
            tx_sym = cell(n_dim,1); % indices of subcodebook constellation points
            n_prune = zeros(n_dim,1); % number of pruned bits from lteTurboEncode output (these bits correspond to tail phase of the encoder)
            interleave = cell(n_dim,1);
            for d_i = 1:n_dim % generate random transmit bits independently for each stream
                tx_bits{d_i} = int8(randi(seed,[0,1],1,KK)); 
                temp_bits = lteTurboEncode(tx_bits{d_i});
                tx_coded_bits{d_i} = temp_bits(1:NN); % throw away trailing bits
                n_prune(d_i) = length(temp_bits)-NN;
%                 rx_bits = lteTurboDecode(temp_bits); % just for debugging purposes
                if activate_interleaving
                    interleave{d_i} = randperm(seed,NN);
                else
                    interleave{d_i} = 1:NN;
                end
                bit_interleaved = tx_coded_bits{d_i}(interleave{d_i});
                bit_mat = reshape(bit_interleaved,bit_vec(d_i),NN/bit_vec(d_i)); 
                sym_inds = sum(bit_2_sym{d_i}.*double(bit_mat),1) + 1;
                tx_sym{d_i} = sym_inds;
                U_tx_sub{d_i} = CB{d_i}(:,:,sym_inds); % subcodebook constellation 
                if d_i > 1
                    U_tx_full = pagemtimes(U_tx_full,U_tx_sub{d_i}); % full transmit symbol - product over subcodebooks
                else
                    U_tx_full = U_tx_sub{d_i};
                end
            end    
%             n_sym = size(U_tx_full,3); % number of symbols per code block
            H = 1/sqrt(2)*(randn(seed,m,Nr,n_sym)+1i*randn(seed,m,Nr,n_sym)); % channel for each symbol (block fading per symbol)
            Y_rx = sqrt(n)/sqrt(m)*pagemtimes(U_tx_full,H) + sqrt(sigma_n2)/sqrt(2)*(randn(seed,n,Nr,n_sym)+1i*randn(seed,n,Nr,n_sym)); % received signal

            % calculate matrix YY before the loop
            YY = zeros(n,n,n_sym);
            for s_i = 1:n_sym % detection is symbol-per-symbol (a bit slow, but not possible otherwise)
                Y_rx_temp = Y_rx(:,:,s_i);
                YY(:,:,s_i) = (Y_rx_temp*Y_rx_temp'); % to avoid multiple calculations of this within the loop further below
            end

%             num_errors = zeros(n_dim,1);
%             uncoded_errors = zeros(n_dim,1);
%             sym_err_count_coded = zeros(n_dim,1);
            CBinds_SIC = cell(n_dim,1);
            % successive detection of streams to implement SIC
            for d_i = 1:n_dim 
                llr_full = []; % LLRs using hypergeometric function   
                llr_full_jensen = [];
                llr_full_max = [];
                
                YY_new = zeros(dim_vec(d_i),dim_vec(d_i),n_sym);
                for s_i = 1:n_sym
                                   
                    if d_i > 1 
                        CB_temp = CB{d_i-1,1}(:,:,CBinds_SIC{d_i-1}(s_i));
%                         CB_temp = CB{d_i-1,1}(:,:,tx_sym{d_i-1}(s_i));
                        YYs = CB_temp'*YY(:,:,s_i)*CB_temp;
                        YYs = 1/2*(YYs + YYs'); % this improves numerical stability as it ensures that YY is non-negative definite (real non-negative eigenvalues)
                        YY_new(:,:,s_i) = YYs;
                    else
                        YYs = YY(:,:,s_i);
                    end
                    n_val = CB_size_vec(d_i);
    
                    % calculate the likelihoods of all symbols
                    lambda_vec = zeros(n_val,1); % vector of eigenvalues for hypergeometric calculation
                    jensen_vec = zeros(n_val,1); % vector of norm-terms for Jensen LLRs
                    for n_i = 1:n_val
                        CB_temp = CB{d_i,1}(:,:,n_i);
                        mat_prod = CB_temp'*YYs*CB_temp;
                        mat_prod = 1/2*(mat_prod + mat_prod');
                        if dim_vec(d_i+1) > 1
                            lambda_vec(n_i) = abs(eigs(mat_prod,1)); % largest eigenvalue -- these should be real positive eigenvalues since they come from a quadratic form
                            jensen_vec(n_i) = 1/dim_vec(d_i+1)*trace(mat_prod);
                        else % no need for eig if it is already a scalar
                            lambda_vec(n_i) = abs(mat_prod);
                            jensen_vec(n_i) = lambda_vec(n_i); % the same for the last stage
                        end
                    end    
                    max_vec = lambda_vec; % largest singular value according to (12) is the same as the largest eigenvalue of the quadratic form
                    like_vec = mgf{d_i}(lambda_vec*rho_term); % vector of likelihoods for hypergeometric calculation   
                    like_vec_jensen = exp(rho_term*jensen_vec); % vector of likelihoods using Jensen
                    like_vec_max = exp(rho_term*max_vec); % vector of likelihoods using max-log
                    
                    llr_vec = zeros(1,bit_vec(d_i)); % LLRs of current symbol
                    max_log = true;
                    if max(like_vec)^(1/sigma_n2) < 1e10 % activate max-log instead of sum when numbers get too large - otherwise everything is Inf
                        max_log = false;
                    end

                    llr_vec_jensen = zeros(1,bit_vec(d_i)); % LLRs of current symbol
                    max_log_jensen = true;
                    if max(like_vec_jensen)^(1/sigma_n2) < 1e10 % activate max-log instead of sum when numbers get too large - otherwise everything is Inf
                        max_log_jensen = false;
                    end

                    llr_vec_max = zeros(1,bit_vec(d_i)); % LLRs of current symbol
                    max_log_max = true;
                    if max(like_vec_max)^(1/sigma_n2) < 1e10 % activate max-log instead of sum when numbers get too large - otherwise everything is Inf
                        max_log_max = false;
                    end
       
                    for b_i = 1:bit_vec(d_i) % LLR calculation for each bit -- here we consider only the subcodebook entries that give the largest trellis metric
                        % LLRs for hypergeometric calculation
                        A = 0; % numerator of LLR
                        B = 0; % denominator of LLR
                        Aj = 0; % numerator with Jensen
                        Bj = 0; % denominator with Jensen
                        Am = 0; % numerator with max-log
                        Bm = 0; % denominator with max-log
                        for n_i = 1:n_val
                            if sum(same_bit_sym{d_i}(b_i,:) == n_i) % current subcodebook entry contributes to 1-bit
                                A = [A;like_vec(n_i)];
                                Aj = [Aj;like_vec_jensen(n_i)];
                                Am = [Am;like_vec_max(n_i)];
                            else % current subcodebook entry contributes to 0-bit
                                B = [B;like_vec(n_i)];
                                Bj = [Bj;like_vec_jensen(n_i)];
                                Bm = [Bm;like_vec_max(n_i)];
                            end
                        end
                        if max_log % max-log approximation of the sum of LLRs
                            A = max(A);
                            B = max(B);
                        else % summing the LLRs
                            A = sum(A.^(1/sigma_n2)); % notice the additional exponentiation - this is necessary as we do not use normalized noise (as used in the Cube-Split paper)
                            B = sum(B.^(1/sigma_n2));
                        end
                        if max_log_jensen % max-log approximation of the sum of LLRs
                            Aj = max(Aj);
                            Bj = max(Bj);
                        else % summing the LLRs
                            Aj = sum(Aj.^(1/sigma_n2)); % notice the additional exponentiation - this is necessary as we do not use normalized noise (as used in the Cube-Split paper)
                            Bj = sum(Bj.^(1/sigma_n2));
                        end
                        if max_log_max % max-log approximation of the sum of LLRs
                            Am = max(Am);
                            Bm = max(Bm);
                        else % summing the LLRs
                            Am = sum(Am.^(1/sigma_n2)); % notice the additional exponentiation - this is necessary as we do not use normalized noise (as used in the Cube-Split paper)
                            Bm = sum(Bm.^(1/sigma_n2));
                        end
    
                        LLR = log(A)-log(B);
                        LLRj = log(Aj)-log(Bj);
                        LLRm = log(Am)-log(Bm);
                        llr_vec(b_i) = LLR;
                        llr_vec_jensen(b_i) = LLRj; 
                        llr_vec_max(b_i) = LLRm; 
                    end
                    llr_full = [llr_full,llr_vec]; % LLRs of all symbols
                    llr_full_jensen = [llr_full_jensen,llr_vec_jensen]; % LLRs of all symbols
                    llr_full_max = [llr_full_max,llr_vec_max]; % LLRs of all symbols
                end 
                if d_i > 1
                    YY = YY_new;
                end

                % hypergeometric LLRs
                llr_temp = llr_full;
                inf_ind = isinf(llr_temp);
                if sum(~inf_ind) % remove Inf LLRs (this happens if A or B is equal to zero - no evidence for 1 or 0 bits)
                    max_val = max(llr_temp(~inf_ind))*10;  % replace Inf with something larger than the largest non-Inf value
                else
                    max_val = 10; % this value is used if all values are +-Inf
                end
                llr_temp(inf_ind) = sign(llr_temp(inf_ind))*max_val; % replace Inf LLRs with something that we can handle numerically
                llr_temp_deinterleave = zeros(size(llr_temp));
                llr_temp_deinterleave(interleave{d_i}) = llr_temp;
                tx_decoded_bits = lteTurboDecode([llr_temp_deinterleave,zeros(1,n_prune(d_i))],10); % append pruned bits and decode
                bit_errors(:,d_i,bb,nn) = (tx_decoded_bits.' ~= tx_bits{d_i});

                % Jensen LLRs
                llr_temp = llr_full_jensen;
                inf_ind = isinf(llr_temp);
                if sum(~inf_ind) % remove Inf LLRs (this happens if A or B is equal to zero - no evidence for 1 or 0 bits)
                    max_val = max(llr_temp(~inf_ind))*10;  % replace Inf with something larger than the largest non-Inf value
                else
                    max_val = 10; % this value is used if all values are +-Inf
                end
                llr_temp(inf_ind) = sign(llr_temp(inf_ind))*max_val; % replace Inf LLRs with something that we can handle numerically
                llr_temp_deinterleave = zeros(size(llr_temp));
                llr_temp_deinterleave(interleave{d_i}) = llr_temp;
                tx_decoded_bits = lteTurboDecode([llr_temp_deinterleave,zeros(1,n_prune(d_i))],10); % append pruned bits and decode
                bit_errors_jensen(:,d_i,bb,nn) = (tx_decoded_bits.' ~= tx_bits{d_i});

                % max-log LLRs
                llr_temp = llr_full_max;
                inf_ind = isinf(llr_temp);
                if sum(~inf_ind) % remove Inf LLRs (this happens if A or B is equal to zero - no evidence for 1 or 0 bits)
                    max_val = max(llr_temp(~inf_ind))*10;  % replace Inf with something larger than the largest non-Inf value
                else
                    max_val = 10; % this value is used if all values are +-Inf
                end
                llr_temp(inf_ind) = sign(llr_temp(inf_ind))*max_val; % replace Inf LLRs with something that we can handle numerically
                llr_temp_deinterleave = zeros(size(llr_temp));
                llr_temp_deinterleave(interleave{d_i}) = llr_temp;
                tx_decoded_bits = lteTurboDecode([llr_temp_deinterleave,zeros(1,n_prune(d_i))],10); % append pruned bits and decode
                bit_errors_max(:,d_i,bb,nn) = (tx_decoded_bits.' ~= tx_bits{d_i});

                % uncoded transmission                
                uncoded_bits = llr_temp_deinterleave(1:KK) > 0; % extract systematic bits
                bit_errors_uncoded(:,d_i,bb,nn) = (uncoded_bits ~= tx_bits{d_i});
%                 uncoded_errors(d_i) = sum(abs(uncoded_bits - double(tx_bits{d_i}))); % number of errors without coding

                % re-encode for SIC
                temp_bits = lteTurboEncode(tx_decoded_bits);
                tx_reencoded_bits = temp_bits(1:NN); % throw away trailing bits
                tx_reencoded_bits = tx_reencoded_bits(interleave{d_i});
                bit_mat = reshape(tx_reencoded_bits,bit_vec(d_i),NN/bit_vec(d_i)); 
                sym_inds = sum(bit_2_sym{d_i}.*double(bit_mat),1) + 1;
                CBinds_SIC{d_i} = sym_inds;
                sym_errors(:,d_i,bb,nn) = (CBinds_SIC{d_i} ~= tx_sym{d_i});
%                 sym_err_count_coded(d_i) = sum(CBinds_SIC{d_i} ~= tx_sym{d_i});
            end
%             sym_error_sum = sym_error_sum + sym_err_count;
%             sym_error_coded_sum = sym_error_coded_sum + sym_err_count_coded;
%             bit_error_sum = bit_error_sum + num_errors;
%             bit_error_uncoded_sum = bit_error_uncoded_sum + uncoded_errors;
        end

    end
%     SER(snr_i,:) = sym_error_sum./(NN_rand*NN./bit_vec.'*BB);
%     SER_coded(snr_i,:) = sym_error_coded_sum./(NN_rand*NN./bit_vec.'*BB);
%     BER(snr_i,:) = bit_error_sum/(NN_rand*KK*BB);
%     BER_uncoded(snr_i,:) = bit_error_uncoded_sum/(NN_rand*KK*BB);

    SER_coded(:,:,:,:,snr_i) = sym_errors;
    BER(:,:,:,:,snr_i) = bit_errors;
    BER_jensen(:,:,:,:,snr_i) = bit_errors_jensen;
    BER_max(:,:,:,:,snr_i) = bit_errors_max;
    BER_uncoded(:,:,:,:,snr_i) = bit_errors_uncoded;

end
SER_stream = squeeze(mean(SER_coded,[1,3,4]));
BER_stream = squeeze(mean(BER,[1,3,4]));
BER_stream_jensen = squeeze(mean(BER_jensen,[1,3,4]));
BER_stream_max = squeeze(mean(BER_max,[1,3,4]));
BER_stream_uncoded = squeeze(mean(BER_uncoded,[1,3,4]));


figure(1) % plot of BERs
SNR_vec = -10*log10(sigma_n2_vec);
semilogy(SNR_vec,mean(BER_stream,1))
hold on
grid on
semilogy(SNR_vec,mean(BER_stream_max,1),'--')
semilogy(SNR_vec,mean(BER_stream_jensen,1),':')
% legend('Hypergeometric','Jensen','Max-log','Uncoded')
xlabel('SNR')
ylabel('BER')

% figure(2) % plot of SER
% semilogy(SNR_vec,mean(SER_stream,1))
% hold on
% grid on
% % semilogy(SNR_vec,SER_stream(1,:),'--')
% semilogy(SNR_vec,SER_stream(2,:),':')
% xlabel('SNR')
% ylabel('SER')
% % legend('Average SER of both streams','SER stream 1','SER stream 2')

file_name = ['SoftSIC_' dim_string '_CB' num2str(sum(log2(CB_size_vec))) '_LLRs.mat'];
save(fullfile(file_name));
