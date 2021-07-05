index = ssget ;

f = find (index.nrows == index.ncols & ...
index.sprank == index.ncols & ...
index.sprank == index.nrows & ...
index.isReal & ~index.isGraph);
%~index.posdef & ...
%index.numerical_symmetry <= .9 ) ;
%index.pattern_symmetry <= .6 ) ;

[ignore, i] = sort (index.nnz (f) + index.nzero (f)) ;

f = f (i) ;

nzlast = -1 ;
fnew = [ ] ;

nmat = length (f) ;
for k = 1:nmat
    id = f (k) ;

    thisnz = index.nnz (id) + index.nzero (id) ;
    if (thisnz ~= nzlast)
        fnew = [fnew id] ;
    end
    nzlast = thisnz ;
end

fnew = fnew' ;
nmat = length (fnew) ;


err = 1e-5;
%Don't scale the matrix if s ==0 scale otherwise
s = 0;

ff = fopen ('myRes.m', 'w') ;

% Headers

%%matlab format
fprintf(ff,'%% id nnzA ' );
fprintf(ff,'myErr umfErr logratio' );
fprintf(ff,' myElaps umfElaps ratio');
fprintf(ff,' mynnz umfnnz ratio');
fprintf(ff,' myflop umfflop ratio');
%fprintf(ff,' hardwareflop ratio/myflop');  %if COUNT_FLOP
fprintf(ff,' \n results = [');

%%csv format
%fprintf(ff,'id, nnzA, myErr, umfErr, logratio,' );
%fprintf(ff,' myElaps, umfElaps, ratio,');
%fprintf(ff,' mynnz, umfnnz, ratio,');
%fprintf(ff,' myflop, umfflop, ratio\n');

loop_cnt = 0;
NNZMat = 100000;

for k = 1:nmat
%for k = 605:nmat
    id = fnew (k); 
    % some problem in these matrice
    if ( id == 2056 || id == 2034 || id == 1867 || id == 2842 || ...
        id == 2843 ||    id == 2844 || id == 2845 || id == 1396 || ...
        id == 1397 || ... %ordering failed
        id == 1404 || id == 1297 || id == 788 || id == 1373 || id == 2265 || ...
        id == 274 ||  id == 273 || id == 2015 || id == 2104 || ...%newer tests singular
        id == 2384 || id == 2385 || id == 2835 || ... %memory!!
        id == 109 ||  id == 102 || id == 113 || id == 93 || id == 119 ||...  %smaller singular
        id == 95 ||  id == 16 || id == 17 || id == 116 || id == 117 ||...  %smaller singular
        id == 96 ||  id == 8 || id == 121 || id == 100 || id == 98 ||...  %smaller singular
        id == 99 ||  id == 123 || id == 114 || id == 128 || id == 125 ||...  %smaller singular
        id == 130 ||  id == 131 || id == 18 || id == 19 || id == 134 ||...  %smaller singular
        id == 20 ||  id == 21 || id == 105 || id == 135 || id == 104 ||...  %smaller singular
        id == 138 ||  id == 137 || id == 111 || id == 110 || id == 112 ||...  %smaller singular
        id == 1186 ||  id == 91 || id == 92 || id == 139 || id == 140 ||...  %smaller singular
        id == 22 ||  id == 250 || id == 105 || id == 135 || id == 104 ||...  %smaller singular
        id == 120||  id == 1185 || id == 1380 || id == 1303 || id == 104 ||...  %smaller singular
        id == 2267 || id == 2649 || id == 2847 || id == 2337 || id == 2841)  
        continue;
    end
    group = index.Group {id} ;
    name = index.Name {id} ;

    Prob = ssget(id);
    Aorig = Prob.A;

    if (nnz(Aorig) > NNZMat)
            continue;
    end

    [dp,dq,dr,ds,dcc,drr] = dmperm(Aorig);
    [m n] = size (Aorig);

    A = Aorig;

    if (size(dr) ~= 2 )
        %continue 
        if (norm(diff(dr)-diff(ds)) ~= 0 )
            sprintf('Unexpected')
            continue;
        end
        B = Aorig(dp,dq);
        [M,I] = max(diff(dr));
        p = dr(I):dr(I+1)-1;
        q = ds(I):ds(I+1)-1;
        %A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );
        A = B(p,p);


        [m n] = size (A);
        if ( m== 1)
            sprintf('not worth trying');
            continue;
        end
    end


    if (nnz(A) > NNZMat)
            continue;
    end

    loop_cnt = loop_cnt + 1;

    if (loop_cnt > 2150)
        break
    end

    id
    %max scaling
    A = spdiags (1./max (abs(A),[], 2), 0, size(A,1), size(A,2)) * A ;
    mmwrite('../Matrix/ParUTst/tmp.mtx', A);
    %    intel = sprintf('. /home/grads/a/aznaveh/intel/bin/compilervars.sh intel64;');
    %    intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
    str = sprintf ('../Demo/x_paru_demo %d < ../Matrix/ParUTst/tmp.mtx', id );
    mypath = '../Demo/Res/';
    %str = strcat(intel, str);

    %%%  {  five times test
    %info_name = sprintf ('%d_info.txt',id);
    %infofullname = strcat(mypath, info_name);

    %for i=1:5
    %    %str = strcat(intel, str);
    %    system(str);
    %    t_Info = load (infofullname);
    %    a_myElaps(i) = t_Info(1);
    %    a_fromCode_umf_Elaps(i) = t_Info(2);
    %end

    %a_myElaps = sort(a_myElaps);
    %a_fromCode_umf_Elaps = sort(a_fromCode_umf_Elaps);
    %a_myElaps = a_myElaps(2:4);
    %a_fromCode_umf_Elaps = a_fromCode_umf_Elaps (2:4);
    %myElaps = mean(a_myElaps);
    %fromCode_umf_Elaps = mean(a_fromCode_umf_Elaps);
    %%%  }  five times test

    %%% { one time test
    info_name = sprintf ('%d_info.txt',id);
    infofullname = strcat(mypath, info_name);
    system(str);
    t_Info = load (infofullname);
    myElaps = t_Info(1);
    fromCode_umf_Elaps = t_Info(2);
    %%% } one time test


    % Loading the results into Matlab
    %mypath = '/export/scratch/multifrontal/aznaveh/myResults/';

    row_name = sprintf ('%d_row.txt',id);
    rowfullname = strcat(mypath, row_name);
    rowp = load (rowfullname);
    rowp = rowp+1;

    col_name = sprintf ('%d_col.txt',id);
    colfullname = strcat(mypath, col_name);
    colp = load (colfullname);
    colp = colp+1;

    if(s == 1)
            s_name = sprintf ('%d_scale.txt', id);
            scalefullname = strcat(mypath, s_name);
            Rvec = load (scalefullname);
            R = spdiags (Rvec, 0, m, m);
    else 
        r = 1;
    end


    LU_name = sprintf ('%d_LU.txt',id);
    LUfullname = strcat(mypath, LU_name);

    %    flp_cnt_dgemm = t_Info(3);
    %    flp_cnt_trsm = t_Info(4);
    %    flp_cnt_dger = t_Info(5);
    %    hardware_flp_cnt = flp_cnt_dgemm + flp_cnt_trsm + flp_cnt_dger;



    [LU, paddingZ] = mmread (LUfullname);


    umfStart= tic;
    %[l,u,p,q,D]=lu(A, 'vector');
    %[l, u, p, q, r, Info]= umfpack (A); % with scaling
    [l, u, p, q ]= umfpack (A);  %no scaling
    %umfElaps = toc(umfStart);
    umfElaps = fromCode_umf_Elaps;


    L=tril(LU,-1)+speye(size(LU));
    U=triu(LU); 

    if (s == 1)
        sA = R \A;
    else
        sA = A;
    end
    myErr = lu_normest(sA(rowp,colp),L,U)/norm(A,1);
    %umfErr = lu_normest(D(:,p)\A(:,q),l,u);
    umfErr = lu_normest(p*(r\A)*q,l,u)/norm(A,1);

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
    fprintf(ff,'%d %d ', id, nnz(A));
    fprintf(ff,'%g %g %g', myErr, umfErr, log10(myErr/umfErr));
    fprintf(ff,' %g %g %g', myElaps, umfElaps, myElaps/umfElaps);
    fprintf(ff,' %g %g %g', mynnz , umfpnnz, mynnz/umfpnnz );
    fprintf(ff,' %g %g %g', myflop, umfflop, myflop/umfflop);

    %fprintf(ff,' %g %g ', hardware_flp_cnt, hardware_flp_cnt/myflop);
    fprintf(ff,' \n');

    %%csv format
    %    fprintf(ff,'%d %d %g %g %g', id, nnz(A), myErr, umfErr, ...
    %    fprintf(ff,'%d, %d, %g, %g, %g,', id, nnz(A), myErr, umfErr, ...
    %        log10(myErr/umfErr));
    %    fprintf(ff,' %g, %g, %g,', myElaps, umfElaps, myElaps/umfElaps);
    %    fprintf(ff,' %g, %g, %g,', mynnz , umfpnnz, mynnz/umfpnnz );
    %    fprintf(ff,' %g, %g, %g\n', myflop, umfflop, myflop/umfflop);



    % cleaning the files because of the memory problem
    str = ['rm  ' mypath LU_name];    system(str);
    str = ['rm  ' mypath col_name];    system(str);
    str = ['rm  ' mypath row_name];    system(str);
    str = ['rm  ' mypath info_name];    system(str);
    if (s == 1)
        str = ['rm  ' mypath s_name];    system(str);
    end
    %id
    %pause;
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
