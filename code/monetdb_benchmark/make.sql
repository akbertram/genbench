-- Create the geo table to hold the expression data
DROP TABLE geo;
CREATE TABLE geo (
       geneid integer,     -- row number in the matrix
       patientid integer,  -- column number in the matrix
       expr_value real     -- expression values
);

-- Create the go_matrix table to hold gene ontology data
DROP TABLE go_matrix;
CREATE TABLE go_matrix (
       geneid integer, -- geneid from the geo table
       goid integer,   -- gene ontology id
       belongs integer -- whether the gene belongs to the go category (either 0 or 1)
);

-- Create the genes table with metadata about
-- each gene
DROP TABLE genes;
CREATE TABLE genes (
       geneid integer, -- geneid from the geo table
       target integer, -- gene target
       pos bigint,     -- position in the genome
       len integer,    -- length of the gene       
       func integer    -- function of the gene
);

-- Create the patients table with metadata about
-- each patient
DROP TABLE patients;
CREATE TABLE patients (
       id integer,      -- patientid from the geo table
       age integer,     -- patient age (in years)
       gender integer,  -- patient gender (0 or 1)
       zipcode integer, -- patient's zip code
       disease integer, -- patient's disease
       response float   -- drug response
);
