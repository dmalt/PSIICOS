test_funciton_name = 'ps.PSIICOS';
test_path = which(test_funciton_name);

if ~isempty(test_path)
    disp('INSTALLATION SUCCESFUL');
else
    disp(['SOMETHING WENT WRONG (cant find ', test_funciton_name,  ' in path)']);
end


