addpath('../../MATLAB_Tools/MESHND');

ff = fopen ('meshndRes.m', 'w') ;
%%matlab format
fprintf(ff,'%% id nnzA ' );
fprintf(ff,'myErr umfErr logratio' );
fprintf(ff,' myElaps umfElaps ratio');
fprintf(ff,' mynnz umfnnz ratio');
fprintf(ff,' myflop umfflop ratio');
%fprintf(ff,' hardwareflop ratio/myflop');  %if COUNT_FLOP
fprintf(ff,' \n results = [');

for n = 2:512
    %generating the matrix
    A = meshsparse(meshnd(n,n));
    %or 
    n
    %A = meshsparse(meshnd(n,n,n));
    spy(A);
    [m nn] = size (A);
    pause(1e-9)

    %writing the matrix and run my code
    mmwrite('../Matrix/ParUTst/tmp.mtx', A);
    my_path = '../Demo/Res/';




    str = sprintf ('../Demo/umfout %d < ../Matrix/ParUTst/tmp.mtx', n );
    info_name = sprintf ('%d_info.txt',n);
    infofullname = strcat(my_path, info_name);


    for i=1:5
        intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
        str = strcat(intel, str);
        system(str);
        t_Info = load (infofullname);
        a_myElaps(i) = t_Info(1);
        a_fromCode_umf_Elaps(i) = t_Info(2);
    end

    a_myElaps = sort(a_myElaps);
    a_fromCode_umf_Elaps = sort(a_fromCode_umf_Elaps);
    a_myElaps = a_myElaps(2:4);
    a_fromCode_umf_Elaps = a_fromCode_umf_Elaps (2:4);
    myElaps = mean(a_myElaps);
    fromCode_umf_Elaps = mean(a_fromCode_umf_Elaps);



    %one time test
    % intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
    % str = strcat(intel, str);
    % system(str);
    % t_Info = load (infofullname);
    % myElaps = t_Info(1);
    % fromCode_umf_Elaps = t_Info(2);



    row_name = sprintf ('%d_row.txt',n);
    rowfullname = strcat(my_path, row_name);
    rowp = load (rowfullname);
    rowp = rowp+1;


    col_name = sprintf ('%d_col.txt',n);
    colfullname = strcat(my_path, col_name);
    colp = load (colfullname);
    colp = colp+1;


    LU_name = sprintf ('%d_LU.txt',n);
    LUfullname = strcat(my_path, LU_name);


    %reading the output
    [LU, paddingZ] = mmread (LUfullname);

    L=tril(LU,-1)+speye(size(LU));
    U=triu(LU); 


    umfStart= tic;
    %[l,u,p,q,D]=lu(A, 'vector');
    %[l, u, p, q, r, Info]= umfpack (A); % with scaling
    [l, u, p, q ]= umfpack (A);  %no scaling
    %umfElaps = toc(umfStart);
    umfElaps = fromCode_umf_Elaps;

    %ratio:
    myElaps/umfElaps

    myErr = lu_normest(A(rowp,colp),L,U)/norm(A,1);
    %umfErr = lu_normest(D(:,p)\A(:,q),l,u);
    umfErr = lu_normest(p*A*q,l,u)/norm(A,1);

    if(myErr <= 100*umfErr || myErr < err)
        fprintf('Pass\n')
    else
         fprintf('Fail\n')
         break;
    end


    umfpnnz= nnz(l)+nnz(u) - m;
    mynnz = nnz(LU); %+ nnz(paddingZ);



    myflop = luflop(L,U);
    umfflop = luflop(l,u);

    %%matlab format
    fprintf(ff,'%d %d ', n, nnz(A));
    fprintf(ff,'%g %g %g', myErr, umfErr, log10(myErr/umfErr));
    fprintf(ff,' %g %g %g', myElaps, umfElaps, myElaps/umfElaps);
    fprintf(ff,' %g %g %g', mynnz , umfpnnz, mynnz/umfpnnz );
    fprintf(ff,' %g %g %g', myflop, umfflop, myflop/umfflop);

    %fprintf(ff,' %g %g ', hardware_flp_cnt, hardware_flp_cnt/myflop);
    fprintf(ff,' \n');
    % cleaning the files because of the memory problem
    str = ['rm  ' my_path LU_name];    system(str);
    str = ['rm  ' my_path col_name];    system(str);
    str = ['rm  ' my_path row_name];    system(str);
    str = ['rm  ' my_path info_name];    system(str);
end

fprintf(ff,'];\n\n');
fprintf(ff,'id = results (:,1) ;\n');
fprintf(ff,'nnzA = results (:,2) ;\n');
fprintf(ff,'\n');

fprintf(ff,'myErr = results (:,3) ;\n');
fprintf(ff,'umfErr = results (:,4) ;\n');
fprintf(ff,'logratio = results (:,5) ;\n');
fprintf(ff,'\n');

fprintf(ff, ...
'noNanRatio = logratio(~any (isnan(logratio) | isinf(logratio),2),:); \n');


fprintf(ff,'myElaps = results (:,6) ;\n');
fprintf(ff,'umfElaps = results (:,7) ;\n');
fprintf(ff,'tratio = results (:,8) ;\n');
fprintf(ff,'\n');

fprintf(ff,'mynnz = results (:,9) ;\n');
fprintf(ff,'umfnnz = results (:,10) ;\n');
fprintf(ff,'nzratio = results (:,11) ;\n');
fprintf(ff,'\n');

fprintf(ff,'myflop = results (:,12) ;\n');
fprintf(ff,'umfflop = results (:,13) ;\n');
fprintf(ff,'flratio = results (:,14) ;\n');
fprintf(ff,'\n');

fprintf(ff,'intensity = umfflop ./ umfnnz ;\n');

fprintf(ff,'figure(1);\n');
fprintf(ff,'subplot (1,3,1) ;\n');
fprintf(ff,'loglog (intensity,  tratio, ''o'', ''MarkerSize'', 10) ;\n');
fprintf(ff,'ylabel (''time ratio'') \n');
fprintf(ff,'xlabel (''intensity'') \n');

fprintf(ff,'subplot (1,3,2) ;\n');
fprintf(ff,'loglog (intensity,  nzratio, ''o'', ''MarkerSize'', 10) ;\n');
fprintf(ff,'ylabel (''nnz ratio'') \n');
fprintf(ff,'xlabel (''intensity'') \n');


fprintf(ff,'subplot (1,3,3) ;\n');
fprintf(ff,'loglog (intensity,  flratio, ''o'', ''MarkerSize'', 10) ;\n');
fprintf(ff,'ylabel (''flops ratio'') \n');
fprintf(ff,'xlabel (''intensity'') \n');



fclose (ff) ;
