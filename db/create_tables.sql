CREATE TABLE meta (
    meta_id INT(6),
    timestamp varchar(30),
    benchmark varchar(30),
    benchmark_group varchar(30),
    sys_name varchar(30),
    sys_release varchar(30),
    lang varchar(30),
    lang_major varchar(30),
    lang_minor varchar(30)
    
);

ALTER TABLE meta 
  ADD CONSTRAINT meta_id_pk
    PRIMARY KEY (meta_id);
ALTER TABLE meta
  MODIFY meta_id INT(6) AUTO_INCREMENT;

ALTER TABLE meta
  ALTER timestamp SET DEFAULT 'not_set'
  
  ;

CREATE TABLE timings (
    meta_id INT,
    block varchar(30) not null,
    user_self NUMERIC(8,3),
    sys_self NUMERIC(8,3),
    elapsed NUMERIC(8,3),
    user_child NUMERIC(8,3),
    sys_child NUMERIC(8,3)
);

ALTER TABLE timings 
  ADD CONSTRAINT meta_id_fk
  FOREIGN KEY (meta_id) 
    REFERENCES meta(meta_id)
    ON DELETE CASCADE;


