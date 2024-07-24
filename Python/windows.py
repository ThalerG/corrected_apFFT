import numpy as np

def nuttall(N, numTerms=4, Order=0):
    """
    Generates a Nuttall window.
    
    Parameters:
    - N: The length of the window.
    - numTerms: The number of terms in the window's sum. Default is 4.
    - Order: The order of the window. Default is 0.
    
    Returns:
    - w: The Nuttall window.
    """
    n = np.arange(0, N)

    if numTerms == 4:
        if Order == 0:
            coefficients = [0.3635819, 0.4891775, 0.1365995, 0.0106411]
        elif Order == 1:
            coefficients = [0.355768, 0.487396, 0.144232, 0.012604]
        elif Order == 3:
            coefficients = [0.338946, 0.481973, 0.161054, 0.018027]
        elif Order == 5:
            coefficients = [0.3125, 0.46875, 0.1875, 0.03125]
        else:
            raise ValueError('Unsupported order')
    elif numTerms == 3:
        if Order == 0:
            coefficients = [0.4243801, 0.4973406, 0.0782793, 0]
        elif Order == 1:
            coefficients = [0.40897, 0.5, 0.09103, 0]
        elif Order == 3:
            coefficients = [0.375, 0.5, 0.125, 0]
        else:
            raise ValueError('Unsupported order')
    else:
        raise ValueError('Unsupported number of terms')

    w = coefficients[0]
    for k in range(1, len(coefficients)):
        w += (-1)**k * coefficients[k] * np.cos(2 * np.pi * k * n / (N-1))

    return w