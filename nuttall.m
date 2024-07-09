function [w] = nuttall(N, varargin)
    %NUTTALL Generates a Nuttall window.
    %
    % Syntax: [y] = nuttall(N, Name, Value)
    %
    % Inputs:
    %   N - The length of the window. Must be a positive integer.
    %
    % Name-Value Pair Arguments:
    %   'numTerms' - The number of terms to use in the window calculation. 
    %                Supported values are 3 and 4. Default is 4.
    %   'Order' - The order of the highest continuous derivative of the window, affecting the shape and mainlobe width. 
    %             Higher values result in a window with a narrower mainlobe and higher sidelobes.
    %             Default is 0, for the minimum window.
    %
    % Outputs:
    %   y - The generated Nuttall window of length N.
    %
    % Example:
    %   N = 100;
    %   w = nuttall(N, 'numTerms', 4, 'Order', 1);
    %   plot(w); title('100-point Nuttall Window');
    %
    % Notes:
    %   The Nuttall window is a sum of cosines that forms a low-pass filter.
    %   It is widely used in signal processing for minimizing spectral leakage.
    %   This function allows customization of the window by specifying the 
    %   number of terms and the sidelobe dropping rate, offering flexibility 
    %   for various signal processing needs.
    %
    %   The coefficients used in the window calculation are determined based 
    %   on the specified 'numTerms' and 'Order'. Unsupported values for 
    %   these parameters will result in an error. The coefficients are based on the
    %   following papers:
    %   
    %   A.H. Nuttall: ‘Some windows with very good sidelobe behaviour’, IEEE Trans. on
    %   Acoustics, Speech and Signal Processing, Vol ASSP-29 (1): 84–91 (1981).
    %
    %
    % Heinzel G. et al., “Spectrum and spectral density estimation by the Discrete Fourier
    % transform (DFT), including a comprehensive list of window functions and some new 
    % flat-top windows”, February 15, 2002 https://holometer.fnal.gov/GH_FFT.pdf
    % See also: hamming, hanning, blackman

    % Create an input parser
    p = inputParser;

    % Define default values for numTerms and sidelobeRate
    defaultNumTerms = 4; % Example default number of terms
    defaultOrder = 0; % Example default window order

    % Add name-value pair arguments for numTerms and sidelobeRate
    addParameter(p, 'numTerms', defaultNumTerms, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'Order', defaultOrder, @(x) isnumeric(x) && isscalar(x));
    
    % Parse the input
    parse(p, varargin{:});

    % Extract the parameters after parsing
    numTerms = p.Results.numTerms;
    Order = p.Results.Order;

    n = 0:1:(N-1);

    switch numTerms
        case 4
            switch Order
                case 0 % Case for minimum window
                    coefficients = [0.3635819, 0.4891775, 0.1365995, 0.0106411];
                case 1 % Case for Order 1
                    coefficients = [0.355768, 0.487396, 0.144232, 0.012604];
                case 3 % Case for Order 3
                    coefficients = [0.338946, 0.481973, 0.161054, 0.018027];
                case 5 % Case for Order 5
                    coefficients = [0.3125, 0.46875, 0.1875, 0.03125];
                otherwise
                    error('Unsupported order');
            end
        case 3
            switch Order
                case 0 % Case for minimum window
                    coefficients = [0.4243801, 0.4973406, 0.0782793, 0];
                case 1 % Case for Order 1
                    coefficients = [0.40897, 0.5, 0.09103, 0];
                case 3 % Case for Order 3
                    coefficients = [0.375, 0.5, 0.125, 0];
                otherwise
                    error('Unsupported order');
            end
        otherwise
            error('Unsupported number of terms');
    end

    w = coefficients(1);
    for k = 1:length(coefficients)-1
        w = w + (-1)^(k)*coefficients(k+1) * cos(2 * pi * k * n ./ (N-1));
    end
end