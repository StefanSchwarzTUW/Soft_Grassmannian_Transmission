function [max_val,trellis_CBinds,CB_searches,flops,Uq_store,trellis_per_stage_inds,per_stage_flops] = Grass_quant_trellis(U,bits_vec,CB3,pruning_percentage,rand_rot)
    U_stage = U;
    trellis_metric = 1;
    Nt = size(U,1);
    Nr = size(U,2);
    trellis_inds = {};
    CB_searches = 0;
    flops = 0;    
    per_stage_flops = [];
    trellis_per_stage_inds = [];
    for d_i = 1:size(CB3,1)
        CB_size_high = max(ceil(2^bits_vec(d_i)),1);
        CB_size_low = max(floor(2^bits_vec(d_i)),1);
        prob = 2^bits_vec(d_i) - CB_size_low;
        bern = rand(1) <= prob;
        CB_size = CB_size_low*(1-bern) + CB_size_high*bern; % in case the codebook size is not an integer, we randomize between the integer below (floor) and above (ceil)        
        prior_trellis_width = size(U_stage,3);
%         trellis_width = ceil(trellis_width*pruning_percentage);
        quante = zeros(CB_size,prior_trellis_width);
        m1 = size(U_stage,1);
        m2 = size(U_stage,2);        
        for u_i = 1:prior_trellis_width            
            if trellis_metric(u_i) > 0 % trellis 
%                 if d_i == 1
                    CB = CB3{d_i,1}(:,:,1:CB_size);
%                 else
%                     CB = CB3{d_i,u_i}(:,:,1:CB_size);
%                 end
                m3 = size(CB,2);
                U = U_stage(:,:,u_i);
                for b_i = 1:CB_size
                    temp = U'*rand_rot{d_i}*CB(:,:,b_i);
                    quante(b_i,u_i) = real(trace(temp*temp'))/Nr; % quantization metric    
                    flops = flops + 2*m2*m1*m3 + 2*m2*m3;
                end
                CB_searches = CB_searches + CB_size;
            end
        end
        quante = quante.*repmat(trellis_metric,CB_size,1);       
        [max_vals,max_inds] = max(quante,[],2);
        [sort_vals,sort_inds] = sort(max_vals,'descend');
        current_trellis_width = ceil(CB_size*pruning_percentage(d_i));
        last_val = sort_vals(current_trellis_width);
%         trellis_metric = max_vals.';
        trellis_inds{d_i} = max_inds;
        U_stage_new = zeros(size(CB,2),Nr,CB_size);
        trellis_metric = zeros(1,CB_size);
        for b_i = 1:CB_size
%             CB = CB3{d_i,max_inds(b_i)}(:,:,1:CB_size);
            if max_vals(b_i) >= last_val
                Uq = rand_rot{d_i}*CB(:,:,b_i);
                U = U_stage(:,:,max_inds(b_i));
                U_stage_new(:,:,b_i) = Uq'*U*(U'*(Uq*Uq')*U)^(-1/2); % this is the SQBC combiner
            	trellis_metric(b_i) = max_vals(b_i);
                flops = flops + 2*m2*m3*m2 + 2*m2^3 + 2*m3*m2^2; % notice U'*Uq has been calculated before during quantization --> can be reused
                % first term is multiplication inside brackets; second term
                % is SVD of matrix inside brackets to calculate the square
                % root and multiplication of the SVD matrices; last term is
                % outer multiplication
            end
        end
        U_stage = U_stage_new;
        [~,max_stage_ind] = max(max_vals);
        trellis_per_stage_inds = [trellis_per_stage_inds,max_stage_ind];
        if ~isempty(per_stage_flops)
            per_stage_flops = [per_stage_flops,flops-per_stage_flops(end-1)];
        end
        if d_i == 1
%             Uq_store = rand_rot{d_i}*CB(:,:,sort_inds(1:2));
           Uq_store = [];
        end
    end
   [max_val,max_ind] = max(max_vals);
   trellis_CBinds = [];
   for t_i = length(trellis_inds):-1:1
       trellis_CBinds = [max_ind,trellis_CBinds];
       max_ind = trellis_inds{t_i}(max_ind);
   end
end