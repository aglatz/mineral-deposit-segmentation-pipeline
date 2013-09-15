function [] = save_xls(SubjectFile_traintest, fmt, data)
SubjectFile_traintest_csv = [SubjectFile_traintest '.csv'];
fd = fopen(SubjectFile_traintest_csv, 'w');
if fd < 0
    error('segment_us_cv:fileopen', 'could not open file');
end
for idx_sub = 1:length(data)
    N_fmt = length(fmt);
    for idx_fmt = 1:N_fmt
        fprintf(fd, [fmt{idx_fmt} ','], data{idx_sub, idx_fmt});
    end
    fprintf(fd, '\n');
end
fclose(fd);
system(sprintf('./text2xls.pl -i %s.csv -o %s.xls', ...
               SubjectFile_traintest, SubjectFile_traintest));
delete(SubjectFile_traintest_csv);