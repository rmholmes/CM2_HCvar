function trends = lintrends(time,signal,window)
% Calculate linear trends in signal over an (odd) window length of
% window (which is in index units).
    
    half_window = (window-1)/2;
    if (half_window ~= round(half_window))
        'Window must be odd!!!';
        return;
    end
    
    trends = NaN*zeros(size(signal));
    
    for ti=half_window+1:(length(time)-half_window)
        [~,b] = linfit(time(ti-half_window:ti+half_window), ...
                       signal(ti-half_window:ti+half_window));
        trends(ti) = b(2);
    end

end
