
function TM = generate_tm(M, N)
    %generates a Transmission Matrix (TM) of a scattering medium
    % M: height of the TM
    % N: weight of the TM
    % Adapted from https://github.com/comediaLKB/NMF_PR/blob/master/generate_tm.m
      
    sz_grains = 4;   %size of a single speckle grain
    if mod(ceil(sqrt(M) / sz_grains), 2) == 0
        n_grains = ceil(sqrt(M) / sz_grains); %number of speckle grains
    else
        n_grains = ceil(sqrt(M) / sz_grains)+1; 
    end  
         
    [Y, X] = meshgrid(1:n_grains, 1:n_grains);
    pupil = (X-n_grains/2).^2 + (Y-n_grains/2).^2 < (n_grains/2)^2; %pupil definition 
    
    % A speckle is modeled by a random phase in the Fourier space
    bruit = exp(2*1i*pi * rand(n_grains, n_grains, N, 'single')); 
    bruit = bruit .* single(pupil);              %imaging system CTF/pupil func
    
    % Fourier transform to go in the object space with zero padding 
    TM = zeros(M, N,  'single');

    for j = 1:N
        Temp = fft2(fftshift(padarray(bruit(:, :, j), ...
            [sqrt(M)/2 - n_grains/2, sqrt(M)/2 - n_grains/2]))); 
        TM(:, j) = 1/M * Temp(:); 
    end
end
