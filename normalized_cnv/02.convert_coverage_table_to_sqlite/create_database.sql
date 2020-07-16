CREATE TABLE "lines" (
	"line_id"	INTEGER PRIMARY KEY,
	"line"	TEXT
);

CREATE TABLE "regions" (
	"chrom"	TEXT,
	"chromStart"	INTEGER,
	"chromEnd"	INTEGER,
	"reg_id"	INTEGER PRIMARY KEY,
	"sd"	REAL
);

CREATE TABLE "norm_covs" (
	"reg_id"	INTEGER,
	"line_id"	INTEGER,
	"norm_cov"	REAL,
	PRIMARY KEY (reg_id, line_id),
	FOREIGN KEY (reg_id) REFERENCES regions (reg_id) 
            ON DELETE CASCADE ON UPDATE NO ACTION,
	FOREIGN KEY (line_id) REFERENCES lines (line_id) 
            ON DELETE CASCADE ON UPDATE NO ACTION
);


CREATE TABLE "covs" (
	"reg_id"	INTEGER,
	"line_id"	INTEGER,
	"cov"	INTEGER,
	PRIMARY KEY (reg_id, line_id),
	FOREIGN KEY (reg_id) REFERENCES regions (reg_id) 
            ON DELETE CASCADE ON UPDATE NO ACTION,
	FOREIGN KEY (line_id) REFERENCES lines (line_id) 
            ON DELETE CASCADE ON UPDATE NO ACTION
);