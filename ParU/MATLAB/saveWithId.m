function A = saveWithId(id)
index = ssget ;

Prob = ssget(id);
A = Prob.A;

[dp,dq,dr,ds,dcc,drr] = dmperm(A);

[m n] = size (A);
if (size(dr) ~= 2 )
    if (norm(diff(dr)-diff(ds)) ~= 0 )
        sprintf('Unexpected')
    end
    B = A(dp,dq);
    [M,I] = max(diff(dr));
    A = B(dr(I):dr(I+1)-1, dr(I):dr(I+1)-1 );

    [m n] = size (A);
    if ( m== 1)
        sprintf('not worth trying');
    end

end


mmwrite('../Matrix/ParUTst/tmp.mtx',A)
end
