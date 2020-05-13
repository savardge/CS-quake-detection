function writeHYPO( fidl,time,ssc, ts, stalst, sc, wt )
%WRITEHYPO Write detection to HYP2000 format in file given
% fidl: file identifier
% time: time string (year, month, day, hour, min of window)
% ssc: time in seconds (window)
% ts: time vector of delays in seconds (actual detection time delays)
% stalst: station list
% sc: Component ( N, E or B (both))
% wt: weight (0 to 2 )

nu = size(stalst,1);

% Print eventline data.
      fprintf(fidl,'%s%5.2f                                                                                                 \n', time,ssc);
 
% Print bogus P time with weight set to 3 (lowest weighting). Need to do this because hypoinverse 
% requires a P time in the file to get an initial starting location. Choose station with smallest S-time
% and subtract 3 seconds to get trial P time. Fred Klein indicates that this will allow hypoinverse
% to get an initial location and then result in P being downweighted to have no effect.
      [dum,itn]=min(ts);
      fprintf(fidl,'%-4s PO ZBHZ IPU3%s%5.2f                                                                             \n', stalst(itn,:),time,ssc+ts(itn)-3.0);
 
% Print station times.
      for it=1:nu
         fprintf(fidl,'%-4s PO ZBHZ     %s            %5.2f%sSU%d                                                              \n',stalst(it,:),time,ssc+ts(it),sc(it),wt(it));
      end

% Print new line to separate event
      fprintf(fidl,' \n');

end

