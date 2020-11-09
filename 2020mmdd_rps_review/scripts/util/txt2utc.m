function uu = txt2utc (tt)
        FMTLIST = {'yyyy-mmm-dd:HH:MM:SS', 'yyyy-mmm-dd:HH:MM', ...
                   'yyyy-mmm-dd:HH', 'yyyy-mmm-dd', 'yyyy-mmm'};
        tmp = [];
        for (ii = 1:length (FMTLIST))
                try
                        tmp = datenum (tt, FMTLIST{ii});
                        break;
                catch
                        continue;
                end;
        end;
        if isempty (tmp)
                error (['Could not interpret UTC time ' tt]);
        end;
        uu = tmp + 48987 - datenum ('Jan-1-1993') + 1;

