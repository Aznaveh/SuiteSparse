clear

index = ssget ;

f = find (index.nrows == index.ncols & ...
index.sprank == index.ncols & ...
~index.posdef & ...
index.isReal & ~index.isGraph  & ...
index.pattern_symmetry <= .6 ) ;

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
fprintf(ff,'%% id nnzA myErr umfErr logratio' );
fprintf(ff,' myElaps umfElaps ratio');
fprintf(ff,' mynnz umfnnz ratio');
fprintf(ff,' myflop umfflop ratio\n results = [');

%%csv format
%fprintf(ff,'id, nnzA, myErr, umfErr, logratio,' );
%fprintf(ff,' myElaps, umfElaps, ratio,');
%fprintf(ff,' mynnz, umfnnz, ratio,');
%fprintf(ff,' myflop, umfflop, ratio\n');

loop_cnt = 0;


%for k = 1:100
for k = 1:nmat
    id = fnew (k); 
    % some problem in these matrice
    if ( id == 2056 || id == 2034 || id == 1867 || id == 2842 || ...
        id == 2843 ||    id == 2844 || id == 2845  ...
        || id == 893) % || id == 823 ||...
%        id == 2232 || id == 826  || id == 1212)  

        continue;
    end
    group = index.Group {id} ;
    name = index.Name {id} ;

    Prob = ssget(id);
    A = Prob.A;


    [dp,dq,dr,ds,dcc,drr] = dmperm(A);

    [m n] = size (A);
    if (size(dr) ~= 2 )
        %continue 
        if (norm(diff(dr)-diff(ds)) ~= 0 )
            sprintf('Unexpected')
            continue;
        end
        B = A(dp,dq);
        [M,I] = max(diff(dr));
        A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );


        [m n] = size (A);
        if ( m== 1)
            sprintf('not worth trying');
            continue;
        end
    end

    if (nnz(A) < 500)
            continue;
    end


    loop_cnt = loop_cnt + 1;
    id

    if (loop_cnt > 10)
        break;
    end

    %max scaling
    A = spdiags (1./max (A,[], 2), 0, size(A,1), size(A,2)) * A ;
    mmwrite('../Matrix/ParUTst/tmp.mtx', A);
    intel = sprintf('. /home/grads/a/aznaveh/intel/bin/compilervars.sh intel64;');
    intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
    str = sprintf ('../Demo/umfout %d < ../Matrix/ParUTst/tmp.mtx', id );
    str = strcat(intel, str);
    system(str);


    % Loading the results into Matlab
    %path = '/export/scratch/multifrontal/aznaveh/myResults/';
    path = '../Demo/Res/';

    row_name = sprintf ('%d_row.txt',id);
    rowfullname = strcat(path, row_name);
    rowp = load (rowfullname);
    rowp = rowp+1;

    col_name = sprintf ('%d_col.txt',id);
    colfullname = strcat(path, col_name);
    colp = load (colfullname);
    colp = colp+1;

    if(s == 1)
            s_name = sprintf ('%d_scale.txt', id);
            scalefullname = strcat(path, s_name);
            Rvec = load (scalefullname);
            R = spdiags (Rvec, 0, m, m);
    else 
        r = 1;
    end


    LU_name = sprintf ('%d_LU.txt',id);
    LUfullname = strcat(path, LU_name);

    info_name = sprintf ('%d_info.txt',id);
    infofullname = strcat(path, info_name);
    t_Info = load (infofullname);
    myElaps = t_Info(1);
    fromCode_umf_Elaps = t_Info(2);

%%    flp_cnt_dgemm = t_Info(3);
%%    flp_cnt_trsm = t_Info(4);
%%    flp_cnt_dger = t_Info(5);
%%    hardware_flp_cnt = flp_cnt_dgemm + flp_cnt_trsm + flp_cnt_dger;



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

    umfpnnz= nnz(l)+nnz(u) - m;
    mynnz = nnz(LU); %+ nnz(paddingZ);



    myflop = luflop(L,U);

    umfflop = luflop(l,u);

    %%matlab format
    fprintf(ff,'%d %d %g %g %g', id, nnz(A), myErr, umfErr, ...
    log10(myErr/umfErr));
    fprintf(ff,' %g %g %g', myElaps, umfElaps, myElaps/umfElaps);
    fprintf(ff,' %g %g %g', mynnz , umfpnnz, mynnz/umfpnnz );
    fprintf(ff,' %g %g %g', myflop, umfflop, myflop/umfflop);

    %% fprintf(ff,' %g ', hardware_flp_cnt);
    fprintf(ff,' \n');

    %%csv format
    %    fprintf(ff,'%d %d %g %g %g', id, nnz(A), myErr, umfErr, ...
    %    fprintf(ff,'%d, %d, %g, %g, %g,', id, nnz(A), myErr, umfErr, ...
    %        log10(myErr/umfErr));
    %    fprintf(ff,' %g, %g, %g,', myElaps, umfElaps, myElaps/umfElaps);
    %    fprintf(ff,' %g, %g, %g,', mynnz , umfpnnz, mynnz/umfpnnz );
    %    fprintf(ff,' %g, %g, %g\n', myflop, umfflop, myflop/umfflop);



    % cleaning the files because of the memory problem
    str = ['rm  ' path LU_name];    system(str);
    str = ['rm  ' path col_name];    system(str);
    str = ['rm  ' path row_name];    system(str);
    str = ['rm  ' path info_name];    system(str);
    if (s == 1)
        str = ['rm  ' path s_name];    system(str);
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

fprintf(ff,'subplot (1,3,2) ;\n');
fprintf(ff,'loglog (intensity,  nzratio, ''o'', ''MarkerSize'', 10) ;\n');

fprintf(ff,'subplot (1,3,3) ;\n');
fprintf(ff,'loglog (intensity,  flratio, ''o'', ''MarkerSize'', 10) ;\n');

fclose (ff) ;
