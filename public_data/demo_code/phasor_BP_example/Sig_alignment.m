function Sig = Sig_alignment(Sig,max_bin)
[num_meas,num_bin] = size(Sig);
if num_bin > max_bin
    Sig = Sig(:,1:max_bin);
elseif num_bin < max_bin
    ext_bin = max_bin - num_bin;
    Sig = [Sig, zeros(num_meas,ext_bin)];
end
end

