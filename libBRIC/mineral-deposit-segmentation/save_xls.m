function [] = save_xls(SubjectFile_traintest, Subjects_train)
SubjectFile_traintest_csv = [SubjectFile_traintest '.csv'];
fd = fopen(SubjectFile_traintest_csv, 'w');
if fd < 0
    error('segment_us_cv:fileopen', 'could not open file');
end
for idx_sub = 1:length(Subjects_train)
    fprintf(fd, '%s,\n', Subjects_train{idx_sub});
end
fclose(fd);
system(sprintf('./text2xls.pl -i %s.csv -o %s.xls', ...
               SubjectFile_traintest, SubjectFile_traintest));
delete(SubjectFile_traintest_csv);