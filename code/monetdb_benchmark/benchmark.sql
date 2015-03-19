-----------------
-- Covariance ---
-----------------
DROP FUNCTION covar;
-- TODO: use non-reshape2 implementation here
CREATE FUNCTION covar(i1 int, i2 int, i3 real) RETURNS TABLE(row_num int, col_num int, val float)  LANGUAGE R {
    library(reshape2)
    x <- data.frame(i1,i2,i3)
    A <- acast(x, list(names(x)[1], names(x)[2]))
    S <- cov(A)
    melt(S)
};

SELECT * FROM covar((SELECT g.patientid, g.geneid, g.expr_value FROM geo g, patients p WHERE g.patientid = p.id AND p.disease = 5 )); 

-- see cstore benchmark for moar