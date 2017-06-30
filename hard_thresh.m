    function b = hard_thresh(a,t)
    %
    %hard_thresh apply hard thresholding at "level" t to matrix a

    b = a .* (abs(a)>t);