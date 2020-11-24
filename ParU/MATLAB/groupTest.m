ids = [ 313 2015];
ids = [1867 2034 2056 2842 2843 2844 2845];
ff = fopen ('groupRes.m', 'w') ;
err = 1e-5;

%%matlab format
fprintf(ff,'%% id nnzA myErr umfErr logratio' );
fprintf(ff,' myElaps umfElaps ratio');
fprintf(ff,' mynnz umfnnz ratio');
fprintf(ff,' myflop umfflop ratio\n results = [');


for k = 1:size(ids,2)

    id = ids(k) 
    Prob = ssget(id);
    Aorig = Prob.A;
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
    A = spdiags (1./max (abs(A),[], 2), 0, size(A,1), size(A,2)) * A ;
    mmwrite('../Matrix/ParUTst/tmp.mtx', A);
    intel = sprintf('. /home/grads/a/aznaveh/intel/bin/compilervars.sh intel64;');
    intel = sprintf('. /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64;');
    str = sprintf ('../Demo/umfout %d < ../Matrix/ParUTst/tmp.mtx', id );
    %str = strcat(intel, str);
    system(str);



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