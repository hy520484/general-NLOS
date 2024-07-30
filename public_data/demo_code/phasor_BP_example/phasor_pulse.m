function P = phasor_pulse(wavelength, cycle_times, c_delta_t, sub_Sig)
    pulse_size = round((cycle_times * wavelength) / c_delta_t);
    cycle_size = round(wavelength / c_delta_t);
    sigma = (cycle_times * wavelength) / 6;
    t = c_delta_t * ((1:pulse_size) - pulse_size/2);
    gaussian_pulse = exp(-(t .* t) / (2 * sigma * sigma));
    sin_wave = sin(2*pi*(1/cycle_size * (1:pulse_size)));
    cos_wave = cos(2*pi*(1/cycle_size * (1:pulse_size)));
    cos_pulse = cos_wave .* gaussian_pulse;
    sin_pulse = sin_wave .* gaussian_pulse;
    P_real = zeros(size(sub_Sig));
    P_imag = zeros(size(sub_Sig));
    for p = 1 : size(sub_Sig, 1)
            P_real(p,:) = conv(sub_Sig(p,:), cos_pulse, 'same');
            P_imag(p,:) = conv(sub_Sig(p,:), sin_pulse, 'same');
    end
    P = complex(P_real, P_imag);
end