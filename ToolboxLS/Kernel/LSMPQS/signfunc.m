    function k = signfunc(w,absgr,delx)
        k = w./((w.^2+(absgr.*delx).^2)).^0.5;
    end