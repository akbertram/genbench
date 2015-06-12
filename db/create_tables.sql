CREATE TABLE meta (
    meta_id INT(6),
    insert_ts TIMESTAMP,
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

ALTER TABLE meta ALTER timestamp SET DEFAULT 'not_set' ;
ALTER TABLE meta ALTER benchmark SET DEFAULT 'not_set';
ALTER TABLE meta ALTER benchmark_group SET DEFAULT 'not_set';
ALTER TABLE meta ALTER sys_name SET DEFAULT 'not_set';
ALTER TABLE meta ALTER sys_release SET DEFAULT 'not_set';
ALTER TABLE meta ALTER lang SET DEFAULT 'not_set';
ALTER TABLE meta ALTER lang_major SET DEFAULT 'not_set';
ALTER TABLE meta ALTER lang_minor SET DEFAULT 'not_set';
ALTER TABLE meta MODIFY COLUMN insert_ts TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP;

CREATE TABLE timings (
    meta_id INT,
    block varchar(30) not null,
    variable varchar(30) not null,
    value NUMERIC(8,3)
);

ALTER TABLE timings 
  ADD CONSTRAINT meta_id_fk
  FOREIGN KEY (meta_id) 
    REFERENCES meta(meta_id)
    ON DELETE CASCADE;

CREATE TABLE extra_meta (
    meta_id INT,
    variable varchar(30) not null,
    value varchar(30)
);

ALTER TABLE extra_meta 
  ADD CONSTRAINT meta_id_fk
  FOREIGN KEY (meta_id) 
    REFERENCES meta(meta_id)
    ON DELETE CASCADE;


