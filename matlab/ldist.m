function d=ldist( Q1, Q2, P )
    d = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
